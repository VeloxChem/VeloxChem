#include "ElectronRepulsionGeom1010ContrRecIDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_idxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_idxx,
                                              const size_t idx_geom_0010_hdxx,
                                              const size_t idx_geom_1010_hdxx,
                                              const size_t idx_geom_1010_hfxx,
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

            const auto hd_geom_0010_off = idx_geom_0010_hdxx + i * dcomps + j;

            auto g_0_0_x_0_xxxxx_xx = cbuffer.data(hd_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxx_xy = cbuffer.data(hd_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxx_xz = cbuffer.data(hd_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxx_yy = cbuffer.data(hd_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxx_yz = cbuffer.data(hd_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxx_zz = cbuffer.data(hd_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_xx = cbuffer.data(hd_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_xy = cbuffer.data(hd_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_xz = cbuffer.data(hd_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_yy = cbuffer.data(hd_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_yz = cbuffer.data(hd_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_zz = cbuffer.data(hd_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_xx = cbuffer.data(hd_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_xy = cbuffer.data(hd_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_xz = cbuffer.data(hd_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_yy = cbuffer.data(hd_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_yz = cbuffer.data(hd_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_zz = cbuffer.data(hd_geom_0010_off + 17 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_xx = cbuffer.data(hd_geom_0010_off + 18 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_xy = cbuffer.data(hd_geom_0010_off + 19 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_xz = cbuffer.data(hd_geom_0010_off + 20 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_yy = cbuffer.data(hd_geom_0010_off + 21 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_yz = cbuffer.data(hd_geom_0010_off + 22 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_zz = cbuffer.data(hd_geom_0010_off + 23 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_xx = cbuffer.data(hd_geom_0010_off + 24 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_xy = cbuffer.data(hd_geom_0010_off + 25 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_xz = cbuffer.data(hd_geom_0010_off + 26 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_yy = cbuffer.data(hd_geom_0010_off + 27 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_yz = cbuffer.data(hd_geom_0010_off + 28 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_zz = cbuffer.data(hd_geom_0010_off + 29 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_xx = cbuffer.data(hd_geom_0010_off + 30 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_xy = cbuffer.data(hd_geom_0010_off + 31 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_xz = cbuffer.data(hd_geom_0010_off + 32 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_yy = cbuffer.data(hd_geom_0010_off + 33 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_yz = cbuffer.data(hd_geom_0010_off + 34 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_zz = cbuffer.data(hd_geom_0010_off + 35 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_xx = cbuffer.data(hd_geom_0010_off + 36 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_xy = cbuffer.data(hd_geom_0010_off + 37 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_xz = cbuffer.data(hd_geom_0010_off + 38 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_yy = cbuffer.data(hd_geom_0010_off + 39 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_yz = cbuffer.data(hd_geom_0010_off + 40 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_zz = cbuffer.data(hd_geom_0010_off + 41 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_xx = cbuffer.data(hd_geom_0010_off + 42 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_xy = cbuffer.data(hd_geom_0010_off + 43 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_xz = cbuffer.data(hd_geom_0010_off + 44 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_yy = cbuffer.data(hd_geom_0010_off + 45 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_yz = cbuffer.data(hd_geom_0010_off + 46 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_zz = cbuffer.data(hd_geom_0010_off + 47 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_xx = cbuffer.data(hd_geom_0010_off + 48 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_xy = cbuffer.data(hd_geom_0010_off + 49 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_xz = cbuffer.data(hd_geom_0010_off + 50 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_yy = cbuffer.data(hd_geom_0010_off + 51 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_yz = cbuffer.data(hd_geom_0010_off + 52 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_zz = cbuffer.data(hd_geom_0010_off + 53 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_xx = cbuffer.data(hd_geom_0010_off + 54 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_xy = cbuffer.data(hd_geom_0010_off + 55 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_xz = cbuffer.data(hd_geom_0010_off + 56 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_yy = cbuffer.data(hd_geom_0010_off + 57 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_yz = cbuffer.data(hd_geom_0010_off + 58 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_zz = cbuffer.data(hd_geom_0010_off + 59 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_xx = cbuffer.data(hd_geom_0010_off + 60 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_xy = cbuffer.data(hd_geom_0010_off + 61 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_xz = cbuffer.data(hd_geom_0010_off + 62 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_yy = cbuffer.data(hd_geom_0010_off + 63 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_yz = cbuffer.data(hd_geom_0010_off + 64 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_zz = cbuffer.data(hd_geom_0010_off + 65 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_xx = cbuffer.data(hd_geom_0010_off + 66 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_xy = cbuffer.data(hd_geom_0010_off + 67 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_xz = cbuffer.data(hd_geom_0010_off + 68 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_yy = cbuffer.data(hd_geom_0010_off + 69 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_yz = cbuffer.data(hd_geom_0010_off + 70 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_zz = cbuffer.data(hd_geom_0010_off + 71 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_xx = cbuffer.data(hd_geom_0010_off + 72 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_xy = cbuffer.data(hd_geom_0010_off + 73 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_xz = cbuffer.data(hd_geom_0010_off + 74 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_yy = cbuffer.data(hd_geom_0010_off + 75 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_yz = cbuffer.data(hd_geom_0010_off + 76 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_zz = cbuffer.data(hd_geom_0010_off + 77 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_xx = cbuffer.data(hd_geom_0010_off + 78 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_xy = cbuffer.data(hd_geom_0010_off + 79 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_xz = cbuffer.data(hd_geom_0010_off + 80 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_yy = cbuffer.data(hd_geom_0010_off + 81 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_yz = cbuffer.data(hd_geom_0010_off + 82 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_zz = cbuffer.data(hd_geom_0010_off + 83 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_xx = cbuffer.data(hd_geom_0010_off + 84 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_xy = cbuffer.data(hd_geom_0010_off + 85 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_xz = cbuffer.data(hd_geom_0010_off + 86 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_yy = cbuffer.data(hd_geom_0010_off + 87 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_yz = cbuffer.data(hd_geom_0010_off + 88 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_zz = cbuffer.data(hd_geom_0010_off + 89 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_xx = cbuffer.data(hd_geom_0010_off + 90 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_xy = cbuffer.data(hd_geom_0010_off + 91 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_xz = cbuffer.data(hd_geom_0010_off + 92 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_yy = cbuffer.data(hd_geom_0010_off + 93 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_yz = cbuffer.data(hd_geom_0010_off + 94 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_zz = cbuffer.data(hd_geom_0010_off + 95 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_xx = cbuffer.data(hd_geom_0010_off + 96 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_xy = cbuffer.data(hd_geom_0010_off + 97 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_xz = cbuffer.data(hd_geom_0010_off + 98 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_yy = cbuffer.data(hd_geom_0010_off + 99 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_yz = cbuffer.data(hd_geom_0010_off + 100 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_zz = cbuffer.data(hd_geom_0010_off + 101 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_xx = cbuffer.data(hd_geom_0010_off + 102 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_xy = cbuffer.data(hd_geom_0010_off + 103 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_xz = cbuffer.data(hd_geom_0010_off + 104 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_yy = cbuffer.data(hd_geom_0010_off + 105 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_yz = cbuffer.data(hd_geom_0010_off + 106 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_zz = cbuffer.data(hd_geom_0010_off + 107 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_xx = cbuffer.data(hd_geom_0010_off + 108 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_xy = cbuffer.data(hd_geom_0010_off + 109 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_xz = cbuffer.data(hd_geom_0010_off + 110 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_yy = cbuffer.data(hd_geom_0010_off + 111 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_yz = cbuffer.data(hd_geom_0010_off + 112 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_zz = cbuffer.data(hd_geom_0010_off + 113 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_xx = cbuffer.data(hd_geom_0010_off + 114 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_xy = cbuffer.data(hd_geom_0010_off + 115 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_xz = cbuffer.data(hd_geom_0010_off + 116 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_yy = cbuffer.data(hd_geom_0010_off + 117 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_yz = cbuffer.data(hd_geom_0010_off + 118 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_zz = cbuffer.data(hd_geom_0010_off + 119 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_xx = cbuffer.data(hd_geom_0010_off + 120 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_xy = cbuffer.data(hd_geom_0010_off + 121 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_xz = cbuffer.data(hd_geom_0010_off + 122 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_yy = cbuffer.data(hd_geom_0010_off + 123 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_yz = cbuffer.data(hd_geom_0010_off + 124 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_zz = cbuffer.data(hd_geom_0010_off + 125 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_xx = cbuffer.data(hd_geom_0010_off + 126 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_xy = cbuffer.data(hd_geom_0010_off + 127 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_xz = cbuffer.data(hd_geom_0010_off + 128 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_yy = cbuffer.data(hd_geom_0010_off + 129 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_yz = cbuffer.data(hd_geom_0010_off + 130 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_zz = cbuffer.data(hd_geom_0010_off + 131 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_xx = cbuffer.data(hd_geom_0010_off + 132 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_xy = cbuffer.data(hd_geom_0010_off + 133 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_xz = cbuffer.data(hd_geom_0010_off + 134 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_yy = cbuffer.data(hd_geom_0010_off + 135 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_yz = cbuffer.data(hd_geom_0010_off + 136 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_zz = cbuffer.data(hd_geom_0010_off + 137 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_xx = cbuffer.data(hd_geom_0010_off + 138 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_xy = cbuffer.data(hd_geom_0010_off + 139 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_xz = cbuffer.data(hd_geom_0010_off + 140 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_yy = cbuffer.data(hd_geom_0010_off + 141 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_yz = cbuffer.data(hd_geom_0010_off + 142 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_zz = cbuffer.data(hd_geom_0010_off + 143 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_xx = cbuffer.data(hd_geom_0010_off + 144 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_xy = cbuffer.data(hd_geom_0010_off + 145 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_xz = cbuffer.data(hd_geom_0010_off + 146 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_yy = cbuffer.data(hd_geom_0010_off + 147 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_yz = cbuffer.data(hd_geom_0010_off + 148 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_zz = cbuffer.data(hd_geom_0010_off + 149 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_xx = cbuffer.data(hd_geom_0010_off + 150 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_xy = cbuffer.data(hd_geom_0010_off + 151 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_xz = cbuffer.data(hd_geom_0010_off + 152 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_yy = cbuffer.data(hd_geom_0010_off + 153 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_yz = cbuffer.data(hd_geom_0010_off + 154 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_zz = cbuffer.data(hd_geom_0010_off + 155 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_xx = cbuffer.data(hd_geom_0010_off + 156 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_xy = cbuffer.data(hd_geom_0010_off + 157 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_xz = cbuffer.data(hd_geom_0010_off + 158 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_yy = cbuffer.data(hd_geom_0010_off + 159 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_yz = cbuffer.data(hd_geom_0010_off + 160 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_zz = cbuffer.data(hd_geom_0010_off + 161 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_xx = cbuffer.data(hd_geom_0010_off + 162 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_xy = cbuffer.data(hd_geom_0010_off + 163 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_xz = cbuffer.data(hd_geom_0010_off + 164 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_yy = cbuffer.data(hd_geom_0010_off + 165 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_yz = cbuffer.data(hd_geom_0010_off + 166 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_zz = cbuffer.data(hd_geom_0010_off + 167 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_xx = cbuffer.data(hd_geom_0010_off + 168 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_xy = cbuffer.data(hd_geom_0010_off + 169 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_xz = cbuffer.data(hd_geom_0010_off + 170 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_yy = cbuffer.data(hd_geom_0010_off + 171 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_yz = cbuffer.data(hd_geom_0010_off + 172 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_zz = cbuffer.data(hd_geom_0010_off + 173 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_xx = cbuffer.data(hd_geom_0010_off + 174 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_xy = cbuffer.data(hd_geom_0010_off + 175 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_xz = cbuffer.data(hd_geom_0010_off + 176 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_yy = cbuffer.data(hd_geom_0010_off + 177 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_yz = cbuffer.data(hd_geom_0010_off + 178 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_zz = cbuffer.data(hd_geom_0010_off + 179 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_xx = cbuffer.data(hd_geom_0010_off + 180 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_xy = cbuffer.data(hd_geom_0010_off + 181 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_xz = cbuffer.data(hd_geom_0010_off + 182 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_yy = cbuffer.data(hd_geom_0010_off + 183 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_yz = cbuffer.data(hd_geom_0010_off + 184 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_zz = cbuffer.data(hd_geom_0010_off + 185 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_xx = cbuffer.data(hd_geom_0010_off + 186 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_xy = cbuffer.data(hd_geom_0010_off + 187 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_xz = cbuffer.data(hd_geom_0010_off + 188 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_yy = cbuffer.data(hd_geom_0010_off + 189 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_yz = cbuffer.data(hd_geom_0010_off + 190 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_zz = cbuffer.data(hd_geom_0010_off + 191 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_xx = cbuffer.data(hd_geom_0010_off + 192 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_xy = cbuffer.data(hd_geom_0010_off + 193 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_xz = cbuffer.data(hd_geom_0010_off + 194 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_yy = cbuffer.data(hd_geom_0010_off + 195 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_yz = cbuffer.data(hd_geom_0010_off + 196 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_zz = cbuffer.data(hd_geom_0010_off + 197 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_xx = cbuffer.data(hd_geom_0010_off + 198 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_xy = cbuffer.data(hd_geom_0010_off + 199 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_xz = cbuffer.data(hd_geom_0010_off + 200 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_yy = cbuffer.data(hd_geom_0010_off + 201 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_yz = cbuffer.data(hd_geom_0010_off + 202 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_zz = cbuffer.data(hd_geom_0010_off + 203 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_xx = cbuffer.data(hd_geom_0010_off + 204 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_xy = cbuffer.data(hd_geom_0010_off + 205 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_xz = cbuffer.data(hd_geom_0010_off + 206 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_yy = cbuffer.data(hd_geom_0010_off + 207 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_yz = cbuffer.data(hd_geom_0010_off + 208 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_zz = cbuffer.data(hd_geom_0010_off + 209 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_xx = cbuffer.data(hd_geom_0010_off + 210 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_xy = cbuffer.data(hd_geom_0010_off + 211 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_xz = cbuffer.data(hd_geom_0010_off + 212 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_yy = cbuffer.data(hd_geom_0010_off + 213 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_yz = cbuffer.data(hd_geom_0010_off + 214 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_zz = cbuffer.data(hd_geom_0010_off + 215 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_xx = cbuffer.data(hd_geom_0010_off + 216 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_xy = cbuffer.data(hd_geom_0010_off + 217 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_xz = cbuffer.data(hd_geom_0010_off + 218 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_yy = cbuffer.data(hd_geom_0010_off + 219 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_yz = cbuffer.data(hd_geom_0010_off + 220 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_zz = cbuffer.data(hd_geom_0010_off + 221 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_xx = cbuffer.data(hd_geom_0010_off + 222 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_xy = cbuffer.data(hd_geom_0010_off + 223 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_xz = cbuffer.data(hd_geom_0010_off + 224 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_yy = cbuffer.data(hd_geom_0010_off + 225 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_yz = cbuffer.data(hd_geom_0010_off + 226 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_zz = cbuffer.data(hd_geom_0010_off + 227 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_xx = cbuffer.data(hd_geom_0010_off + 228 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_xy = cbuffer.data(hd_geom_0010_off + 229 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_xz = cbuffer.data(hd_geom_0010_off + 230 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_yy = cbuffer.data(hd_geom_0010_off + 231 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_yz = cbuffer.data(hd_geom_0010_off + 232 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_zz = cbuffer.data(hd_geom_0010_off + 233 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_xx = cbuffer.data(hd_geom_0010_off + 234 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_xy = cbuffer.data(hd_geom_0010_off + 235 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_xz = cbuffer.data(hd_geom_0010_off + 236 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_yy = cbuffer.data(hd_geom_0010_off + 237 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_yz = cbuffer.data(hd_geom_0010_off + 238 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_zz = cbuffer.data(hd_geom_0010_off + 239 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_xx = cbuffer.data(hd_geom_0010_off + 240 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_xy = cbuffer.data(hd_geom_0010_off + 241 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_xz = cbuffer.data(hd_geom_0010_off + 242 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_yy = cbuffer.data(hd_geom_0010_off + 243 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_yz = cbuffer.data(hd_geom_0010_off + 244 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_zz = cbuffer.data(hd_geom_0010_off + 245 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_xx = cbuffer.data(hd_geom_0010_off + 246 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_xy = cbuffer.data(hd_geom_0010_off + 247 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_xz = cbuffer.data(hd_geom_0010_off + 248 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_yy = cbuffer.data(hd_geom_0010_off + 249 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_yz = cbuffer.data(hd_geom_0010_off + 250 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_zz = cbuffer.data(hd_geom_0010_off + 251 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_xx = cbuffer.data(hd_geom_0010_off + 252 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_xy = cbuffer.data(hd_geom_0010_off + 253 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_xz = cbuffer.data(hd_geom_0010_off + 254 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_yy = cbuffer.data(hd_geom_0010_off + 255 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_yz = cbuffer.data(hd_geom_0010_off + 256 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_zz = cbuffer.data(hd_geom_0010_off + 257 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_xx = cbuffer.data(hd_geom_0010_off + 258 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_xy = cbuffer.data(hd_geom_0010_off + 259 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_xz = cbuffer.data(hd_geom_0010_off + 260 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_yy = cbuffer.data(hd_geom_0010_off + 261 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_yz = cbuffer.data(hd_geom_0010_off + 262 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_zz = cbuffer.data(hd_geom_0010_off + 263 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_xx = cbuffer.data(hd_geom_0010_off + 264 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_xy = cbuffer.data(hd_geom_0010_off + 265 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_xz = cbuffer.data(hd_geom_0010_off + 266 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_yy = cbuffer.data(hd_geom_0010_off + 267 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_yz = cbuffer.data(hd_geom_0010_off + 268 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_zz = cbuffer.data(hd_geom_0010_off + 269 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_xx = cbuffer.data(hd_geom_0010_off + 270 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_xy = cbuffer.data(hd_geom_0010_off + 271 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_xz = cbuffer.data(hd_geom_0010_off + 272 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_yy = cbuffer.data(hd_geom_0010_off + 273 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_yz = cbuffer.data(hd_geom_0010_off + 274 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_zz = cbuffer.data(hd_geom_0010_off + 275 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_xx = cbuffer.data(hd_geom_0010_off + 276 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_xy = cbuffer.data(hd_geom_0010_off + 277 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_xz = cbuffer.data(hd_geom_0010_off + 278 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_yy = cbuffer.data(hd_geom_0010_off + 279 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_yz = cbuffer.data(hd_geom_0010_off + 280 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_zz = cbuffer.data(hd_geom_0010_off + 281 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_xx = cbuffer.data(hd_geom_0010_off + 282 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_xy = cbuffer.data(hd_geom_0010_off + 283 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_xz = cbuffer.data(hd_geom_0010_off + 284 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_yy = cbuffer.data(hd_geom_0010_off + 285 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_yz = cbuffer.data(hd_geom_0010_off + 286 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_zz = cbuffer.data(hd_geom_0010_off + 287 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_xx = cbuffer.data(hd_geom_0010_off + 288 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_xy = cbuffer.data(hd_geom_0010_off + 289 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_xz = cbuffer.data(hd_geom_0010_off + 290 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_yy = cbuffer.data(hd_geom_0010_off + 291 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_yz = cbuffer.data(hd_geom_0010_off + 292 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_zz = cbuffer.data(hd_geom_0010_off + 293 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_xx = cbuffer.data(hd_geom_0010_off + 294 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_xy = cbuffer.data(hd_geom_0010_off + 295 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_xz = cbuffer.data(hd_geom_0010_off + 296 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_yy = cbuffer.data(hd_geom_0010_off + 297 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_yz = cbuffer.data(hd_geom_0010_off + 298 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_zz = cbuffer.data(hd_geom_0010_off + 299 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_xx = cbuffer.data(hd_geom_0010_off + 300 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_xy = cbuffer.data(hd_geom_0010_off + 301 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_xz = cbuffer.data(hd_geom_0010_off + 302 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_yy = cbuffer.data(hd_geom_0010_off + 303 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_yz = cbuffer.data(hd_geom_0010_off + 304 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_zz = cbuffer.data(hd_geom_0010_off + 305 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_xx = cbuffer.data(hd_geom_0010_off + 306 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_xy = cbuffer.data(hd_geom_0010_off + 307 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_xz = cbuffer.data(hd_geom_0010_off + 308 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_yy = cbuffer.data(hd_geom_0010_off + 309 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_yz = cbuffer.data(hd_geom_0010_off + 310 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_zz = cbuffer.data(hd_geom_0010_off + 311 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_xx = cbuffer.data(hd_geom_0010_off + 312 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_xy = cbuffer.data(hd_geom_0010_off + 313 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_xz = cbuffer.data(hd_geom_0010_off + 314 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_yy = cbuffer.data(hd_geom_0010_off + 315 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_yz = cbuffer.data(hd_geom_0010_off + 316 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_zz = cbuffer.data(hd_geom_0010_off + 317 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_xx = cbuffer.data(hd_geom_0010_off + 318 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_xy = cbuffer.data(hd_geom_0010_off + 319 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_xz = cbuffer.data(hd_geom_0010_off + 320 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_yy = cbuffer.data(hd_geom_0010_off + 321 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_yz = cbuffer.data(hd_geom_0010_off + 322 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_zz = cbuffer.data(hd_geom_0010_off + 323 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_xx = cbuffer.data(hd_geom_0010_off + 324 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_xy = cbuffer.data(hd_geom_0010_off + 325 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_xz = cbuffer.data(hd_geom_0010_off + 326 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_yy = cbuffer.data(hd_geom_0010_off + 327 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_yz = cbuffer.data(hd_geom_0010_off + 328 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_zz = cbuffer.data(hd_geom_0010_off + 329 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_xx = cbuffer.data(hd_geom_0010_off + 330 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_xy = cbuffer.data(hd_geom_0010_off + 331 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_xz = cbuffer.data(hd_geom_0010_off + 332 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_yy = cbuffer.data(hd_geom_0010_off + 333 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_yz = cbuffer.data(hd_geom_0010_off + 334 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_zz = cbuffer.data(hd_geom_0010_off + 335 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_xx = cbuffer.data(hd_geom_0010_off + 336 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_xy = cbuffer.data(hd_geom_0010_off + 337 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_xz = cbuffer.data(hd_geom_0010_off + 338 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_yy = cbuffer.data(hd_geom_0010_off + 339 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_yz = cbuffer.data(hd_geom_0010_off + 340 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_zz = cbuffer.data(hd_geom_0010_off + 341 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_xx = cbuffer.data(hd_geom_0010_off + 342 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_xy = cbuffer.data(hd_geom_0010_off + 343 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_xz = cbuffer.data(hd_geom_0010_off + 344 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_yy = cbuffer.data(hd_geom_0010_off + 345 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_yz = cbuffer.data(hd_geom_0010_off + 346 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_zz = cbuffer.data(hd_geom_0010_off + 347 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_xx = cbuffer.data(hd_geom_0010_off + 348 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_xy = cbuffer.data(hd_geom_0010_off + 349 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_xz = cbuffer.data(hd_geom_0010_off + 350 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_yy = cbuffer.data(hd_geom_0010_off + 351 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_yz = cbuffer.data(hd_geom_0010_off + 352 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_zz = cbuffer.data(hd_geom_0010_off + 353 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_xx = cbuffer.data(hd_geom_0010_off + 354 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_xy = cbuffer.data(hd_geom_0010_off + 355 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_xz = cbuffer.data(hd_geom_0010_off + 356 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_yy = cbuffer.data(hd_geom_0010_off + 357 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_yz = cbuffer.data(hd_geom_0010_off + 358 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_zz = cbuffer.data(hd_geom_0010_off + 359 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_xx = cbuffer.data(hd_geom_0010_off + 360 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_xy = cbuffer.data(hd_geom_0010_off + 361 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_xz = cbuffer.data(hd_geom_0010_off + 362 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_yy = cbuffer.data(hd_geom_0010_off + 363 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_yz = cbuffer.data(hd_geom_0010_off + 364 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_zz = cbuffer.data(hd_geom_0010_off + 365 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_xx = cbuffer.data(hd_geom_0010_off + 366 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_xy = cbuffer.data(hd_geom_0010_off + 367 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_xz = cbuffer.data(hd_geom_0010_off + 368 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_yy = cbuffer.data(hd_geom_0010_off + 369 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_yz = cbuffer.data(hd_geom_0010_off + 370 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_zz = cbuffer.data(hd_geom_0010_off + 371 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_xx = cbuffer.data(hd_geom_0010_off + 372 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_xy = cbuffer.data(hd_geom_0010_off + 373 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_xz = cbuffer.data(hd_geom_0010_off + 374 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_yy = cbuffer.data(hd_geom_0010_off + 375 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_yz = cbuffer.data(hd_geom_0010_off + 376 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_zz = cbuffer.data(hd_geom_0010_off + 377 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : HFSS

            const auto hf_geom_1010_off = idx_geom_1010_hfxx + i * dcomps + j;

            auto g_x_0_x_0_xxxxx_xxx = cbuffer.data(hf_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_xxy = cbuffer.data(hf_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_xxz = cbuffer.data(hf_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_xyy = cbuffer.data(hf_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_xyz = cbuffer.data(hf_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_xzz = cbuffer.data(hf_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_yyy = cbuffer.data(hf_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_yyz = cbuffer.data(hf_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_yzz = cbuffer.data(hf_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_zzz = cbuffer.data(hf_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_xxx = cbuffer.data(hf_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_xxy = cbuffer.data(hf_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_xxz = cbuffer.data(hf_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_xyy = cbuffer.data(hf_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_xyz = cbuffer.data(hf_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_xzz = cbuffer.data(hf_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_yyy = cbuffer.data(hf_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_yyz = cbuffer.data(hf_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_yzz = cbuffer.data(hf_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_zzz = cbuffer.data(hf_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_xxx = cbuffer.data(hf_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_xxy = cbuffer.data(hf_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_xxz = cbuffer.data(hf_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_xyy = cbuffer.data(hf_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_xyz = cbuffer.data(hf_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_xzz = cbuffer.data(hf_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_yyy = cbuffer.data(hf_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_yyz = cbuffer.data(hf_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_yzz = cbuffer.data(hf_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_zzz = cbuffer.data(hf_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_xxx = cbuffer.data(hf_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_xxy = cbuffer.data(hf_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_xxz = cbuffer.data(hf_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_xyy = cbuffer.data(hf_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_xyz = cbuffer.data(hf_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_xzz = cbuffer.data(hf_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_yyy = cbuffer.data(hf_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_yyz = cbuffer.data(hf_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_yzz = cbuffer.data(hf_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_zzz = cbuffer.data(hf_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_xxx = cbuffer.data(hf_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_xxy = cbuffer.data(hf_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_xxz = cbuffer.data(hf_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_xyy = cbuffer.data(hf_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_xyz = cbuffer.data(hf_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_xzz = cbuffer.data(hf_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_yyy = cbuffer.data(hf_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_yyz = cbuffer.data(hf_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_yzz = cbuffer.data(hf_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_zzz = cbuffer.data(hf_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_xxx = cbuffer.data(hf_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_xxy = cbuffer.data(hf_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_xxz = cbuffer.data(hf_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_xyy = cbuffer.data(hf_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_xyz = cbuffer.data(hf_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_xzz = cbuffer.data(hf_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_yyy = cbuffer.data(hf_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_yyz = cbuffer.data(hf_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_yzz = cbuffer.data(hf_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_zzz = cbuffer.data(hf_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_xxx = cbuffer.data(hf_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_xxy = cbuffer.data(hf_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_xxz = cbuffer.data(hf_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_xyy = cbuffer.data(hf_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_xyz = cbuffer.data(hf_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_xzz = cbuffer.data(hf_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_yyy = cbuffer.data(hf_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_yyz = cbuffer.data(hf_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_yzz = cbuffer.data(hf_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_zzz = cbuffer.data(hf_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_xxx = cbuffer.data(hf_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_xxy = cbuffer.data(hf_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_xxz = cbuffer.data(hf_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_xyy = cbuffer.data(hf_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_xyz = cbuffer.data(hf_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_xzz = cbuffer.data(hf_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_yyy = cbuffer.data(hf_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_yyz = cbuffer.data(hf_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_yzz = cbuffer.data(hf_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_zzz = cbuffer.data(hf_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_xxx = cbuffer.data(hf_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_xxy = cbuffer.data(hf_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_xxz = cbuffer.data(hf_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_xyy = cbuffer.data(hf_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_xyz = cbuffer.data(hf_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_xzz = cbuffer.data(hf_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_yyy = cbuffer.data(hf_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_yyz = cbuffer.data(hf_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_yzz = cbuffer.data(hf_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_zzz = cbuffer.data(hf_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_xxx = cbuffer.data(hf_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_xxy = cbuffer.data(hf_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_xxz = cbuffer.data(hf_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_xyy = cbuffer.data(hf_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_xyz = cbuffer.data(hf_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_xzz = cbuffer.data(hf_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_yyy = cbuffer.data(hf_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_yyz = cbuffer.data(hf_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_yzz = cbuffer.data(hf_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_zzz = cbuffer.data(hf_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_xxx = cbuffer.data(hf_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_xxy = cbuffer.data(hf_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_xxz = cbuffer.data(hf_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_xyy = cbuffer.data(hf_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_xyz = cbuffer.data(hf_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_xzz = cbuffer.data(hf_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_yyy = cbuffer.data(hf_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_yyz = cbuffer.data(hf_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_yzz = cbuffer.data(hf_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_zzz = cbuffer.data(hf_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_xxx = cbuffer.data(hf_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_xxy = cbuffer.data(hf_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_xxz = cbuffer.data(hf_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_xyy = cbuffer.data(hf_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_xyz = cbuffer.data(hf_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_xzz = cbuffer.data(hf_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_yyy = cbuffer.data(hf_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_yyz = cbuffer.data(hf_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_yzz = cbuffer.data(hf_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_zzz = cbuffer.data(hf_geom_1010_off + 119 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_xxx = cbuffer.data(hf_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_xxy = cbuffer.data(hf_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_xxz = cbuffer.data(hf_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_xyy = cbuffer.data(hf_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_xyz = cbuffer.data(hf_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_xzz = cbuffer.data(hf_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_yyy = cbuffer.data(hf_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_yyz = cbuffer.data(hf_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_yzz = cbuffer.data(hf_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_zzz = cbuffer.data(hf_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_xxx = cbuffer.data(hf_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_xxy = cbuffer.data(hf_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_xxz = cbuffer.data(hf_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_xyy = cbuffer.data(hf_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_xyz = cbuffer.data(hf_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_xzz = cbuffer.data(hf_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_yyy = cbuffer.data(hf_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_yyz = cbuffer.data(hf_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_yzz = cbuffer.data(hf_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_zzz = cbuffer.data(hf_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_xxx = cbuffer.data(hf_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_xxy = cbuffer.data(hf_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_xxz = cbuffer.data(hf_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_xyy = cbuffer.data(hf_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_xyz = cbuffer.data(hf_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_xzz = cbuffer.data(hf_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_yyy = cbuffer.data(hf_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_yyz = cbuffer.data(hf_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_yzz = cbuffer.data(hf_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_zzz = cbuffer.data(hf_geom_1010_off + 149 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_xxx = cbuffer.data(hf_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_xxy = cbuffer.data(hf_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_xxz = cbuffer.data(hf_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_xyy = cbuffer.data(hf_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_xyz = cbuffer.data(hf_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_xzz = cbuffer.data(hf_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_yyy = cbuffer.data(hf_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_yyz = cbuffer.data(hf_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_yzz = cbuffer.data(hf_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_zzz = cbuffer.data(hf_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_xxx = cbuffer.data(hf_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_xxy = cbuffer.data(hf_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_xxz = cbuffer.data(hf_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_xyy = cbuffer.data(hf_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_xyz = cbuffer.data(hf_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_xzz = cbuffer.data(hf_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_yyy = cbuffer.data(hf_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_yyz = cbuffer.data(hf_geom_1010_off + 167 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_yzz = cbuffer.data(hf_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_zzz = cbuffer.data(hf_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_xxx = cbuffer.data(hf_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_xxy = cbuffer.data(hf_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_xxz = cbuffer.data(hf_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_xyy = cbuffer.data(hf_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_xyz = cbuffer.data(hf_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_xzz = cbuffer.data(hf_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_yyy = cbuffer.data(hf_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_yyz = cbuffer.data(hf_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_yzz = cbuffer.data(hf_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_zzz = cbuffer.data(hf_geom_1010_off + 179 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_xxx = cbuffer.data(hf_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_xxy = cbuffer.data(hf_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_xxz = cbuffer.data(hf_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_xyy = cbuffer.data(hf_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_xyz = cbuffer.data(hf_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_xzz = cbuffer.data(hf_geom_1010_off + 185 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_yyy = cbuffer.data(hf_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_yyz = cbuffer.data(hf_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_yzz = cbuffer.data(hf_geom_1010_off + 188 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_zzz = cbuffer.data(hf_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_xxx = cbuffer.data(hf_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_xxy = cbuffer.data(hf_geom_1010_off + 191 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_xxz = cbuffer.data(hf_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_xyy = cbuffer.data(hf_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_xyz = cbuffer.data(hf_geom_1010_off + 194 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_xzz = cbuffer.data(hf_geom_1010_off + 195 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_yyy = cbuffer.data(hf_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_yyz = cbuffer.data(hf_geom_1010_off + 197 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_yzz = cbuffer.data(hf_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_zzz = cbuffer.data(hf_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_xxx = cbuffer.data(hf_geom_1010_off + 200 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_xxy = cbuffer.data(hf_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_xxz = cbuffer.data(hf_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_xyy = cbuffer.data(hf_geom_1010_off + 203 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_xyz = cbuffer.data(hf_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_xzz = cbuffer.data(hf_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_yyy = cbuffer.data(hf_geom_1010_off + 206 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_yyz = cbuffer.data(hf_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_yzz = cbuffer.data(hf_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_zzz = cbuffer.data(hf_geom_1010_off + 209 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_xxx = cbuffer.data(hf_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_xxy = cbuffer.data(hf_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_xxz = cbuffer.data(hf_geom_1010_off + 212 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_xyy = cbuffer.data(hf_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_xyz = cbuffer.data(hf_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_xzz = cbuffer.data(hf_geom_1010_off + 215 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_yyy = cbuffer.data(hf_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_yyz = cbuffer.data(hf_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_yzz = cbuffer.data(hf_geom_1010_off + 218 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_zzz = cbuffer.data(hf_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_xxx = cbuffer.data(hf_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_xxy = cbuffer.data(hf_geom_1010_off + 221 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_xxz = cbuffer.data(hf_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_xyy = cbuffer.data(hf_geom_1010_off + 223 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_xyz = cbuffer.data(hf_geom_1010_off + 224 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_xzz = cbuffer.data(hf_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_yyy = cbuffer.data(hf_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_yyz = cbuffer.data(hf_geom_1010_off + 227 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_yzz = cbuffer.data(hf_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_zzz = cbuffer.data(hf_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_xxx = cbuffer.data(hf_geom_1010_off + 230 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_xxy = cbuffer.data(hf_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_xxz = cbuffer.data(hf_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_xyy = cbuffer.data(hf_geom_1010_off + 233 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_xyz = cbuffer.data(hf_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_xzz = cbuffer.data(hf_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_yyy = cbuffer.data(hf_geom_1010_off + 236 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_yyz = cbuffer.data(hf_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_yzz = cbuffer.data(hf_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_zzz = cbuffer.data(hf_geom_1010_off + 239 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_xxx = cbuffer.data(hf_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_xxy = cbuffer.data(hf_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_xxz = cbuffer.data(hf_geom_1010_off + 242 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_xyy = cbuffer.data(hf_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_xyz = cbuffer.data(hf_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_xzz = cbuffer.data(hf_geom_1010_off + 245 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_yyy = cbuffer.data(hf_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_yyz = cbuffer.data(hf_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_yzz = cbuffer.data(hf_geom_1010_off + 248 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_zzz = cbuffer.data(hf_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_xxx = cbuffer.data(hf_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_xxy = cbuffer.data(hf_geom_1010_off + 251 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_xxz = cbuffer.data(hf_geom_1010_off + 252 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_xyy = cbuffer.data(hf_geom_1010_off + 253 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_xyz = cbuffer.data(hf_geom_1010_off + 254 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_xzz = cbuffer.data(hf_geom_1010_off + 255 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_yyy = cbuffer.data(hf_geom_1010_off + 256 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_yyz = cbuffer.data(hf_geom_1010_off + 257 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_yzz = cbuffer.data(hf_geom_1010_off + 258 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_zzz = cbuffer.data(hf_geom_1010_off + 259 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_xxx = cbuffer.data(hf_geom_1010_off + 260 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_xxy = cbuffer.data(hf_geom_1010_off + 261 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_xxz = cbuffer.data(hf_geom_1010_off + 262 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_xyy = cbuffer.data(hf_geom_1010_off + 263 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_xyz = cbuffer.data(hf_geom_1010_off + 264 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_xzz = cbuffer.data(hf_geom_1010_off + 265 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_yyy = cbuffer.data(hf_geom_1010_off + 266 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_yyz = cbuffer.data(hf_geom_1010_off + 267 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_yzz = cbuffer.data(hf_geom_1010_off + 268 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_zzz = cbuffer.data(hf_geom_1010_off + 269 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_xxx = cbuffer.data(hf_geom_1010_off + 270 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_xxy = cbuffer.data(hf_geom_1010_off + 271 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_xxz = cbuffer.data(hf_geom_1010_off + 272 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_xyy = cbuffer.data(hf_geom_1010_off + 273 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_xyz = cbuffer.data(hf_geom_1010_off + 274 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_xzz = cbuffer.data(hf_geom_1010_off + 275 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_yyy = cbuffer.data(hf_geom_1010_off + 276 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_yyz = cbuffer.data(hf_geom_1010_off + 277 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_yzz = cbuffer.data(hf_geom_1010_off + 278 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_zzz = cbuffer.data(hf_geom_1010_off + 279 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_xxx = cbuffer.data(hf_geom_1010_off + 280 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_xxy = cbuffer.data(hf_geom_1010_off + 281 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_xxz = cbuffer.data(hf_geom_1010_off + 282 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_xyy = cbuffer.data(hf_geom_1010_off + 283 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_xyz = cbuffer.data(hf_geom_1010_off + 284 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_xzz = cbuffer.data(hf_geom_1010_off + 285 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_yyy = cbuffer.data(hf_geom_1010_off + 286 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_yyz = cbuffer.data(hf_geom_1010_off + 287 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_yzz = cbuffer.data(hf_geom_1010_off + 288 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_zzz = cbuffer.data(hf_geom_1010_off + 289 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_xxx = cbuffer.data(hf_geom_1010_off + 290 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_xxy = cbuffer.data(hf_geom_1010_off + 291 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_xxz = cbuffer.data(hf_geom_1010_off + 292 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_xyy = cbuffer.data(hf_geom_1010_off + 293 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_xyz = cbuffer.data(hf_geom_1010_off + 294 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_xzz = cbuffer.data(hf_geom_1010_off + 295 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_yyy = cbuffer.data(hf_geom_1010_off + 296 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_yyz = cbuffer.data(hf_geom_1010_off + 297 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_yzz = cbuffer.data(hf_geom_1010_off + 298 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_zzz = cbuffer.data(hf_geom_1010_off + 299 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_xxx = cbuffer.data(hf_geom_1010_off + 300 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_xxy = cbuffer.data(hf_geom_1010_off + 301 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_xxz = cbuffer.data(hf_geom_1010_off + 302 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_xyy = cbuffer.data(hf_geom_1010_off + 303 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_xyz = cbuffer.data(hf_geom_1010_off + 304 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_xzz = cbuffer.data(hf_geom_1010_off + 305 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_yyy = cbuffer.data(hf_geom_1010_off + 306 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_yyz = cbuffer.data(hf_geom_1010_off + 307 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_yzz = cbuffer.data(hf_geom_1010_off + 308 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_zzz = cbuffer.data(hf_geom_1010_off + 309 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_xxx = cbuffer.data(hf_geom_1010_off + 310 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_xxy = cbuffer.data(hf_geom_1010_off + 311 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_xxz = cbuffer.data(hf_geom_1010_off + 312 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_xyy = cbuffer.data(hf_geom_1010_off + 313 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_xyz = cbuffer.data(hf_geom_1010_off + 314 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_xzz = cbuffer.data(hf_geom_1010_off + 315 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_yyy = cbuffer.data(hf_geom_1010_off + 316 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_yyz = cbuffer.data(hf_geom_1010_off + 317 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_yzz = cbuffer.data(hf_geom_1010_off + 318 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_zzz = cbuffer.data(hf_geom_1010_off + 319 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_xxx = cbuffer.data(hf_geom_1010_off + 320 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_xxy = cbuffer.data(hf_geom_1010_off + 321 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_xxz = cbuffer.data(hf_geom_1010_off + 322 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_xyy = cbuffer.data(hf_geom_1010_off + 323 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_xyz = cbuffer.data(hf_geom_1010_off + 324 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_xzz = cbuffer.data(hf_geom_1010_off + 325 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_yyy = cbuffer.data(hf_geom_1010_off + 326 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_yyz = cbuffer.data(hf_geom_1010_off + 327 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_yzz = cbuffer.data(hf_geom_1010_off + 328 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_zzz = cbuffer.data(hf_geom_1010_off + 329 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_xxx = cbuffer.data(hf_geom_1010_off + 330 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_xxy = cbuffer.data(hf_geom_1010_off + 331 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_xxz = cbuffer.data(hf_geom_1010_off + 332 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_xyy = cbuffer.data(hf_geom_1010_off + 333 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_xyz = cbuffer.data(hf_geom_1010_off + 334 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_xzz = cbuffer.data(hf_geom_1010_off + 335 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_yyy = cbuffer.data(hf_geom_1010_off + 336 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_yyz = cbuffer.data(hf_geom_1010_off + 337 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_yzz = cbuffer.data(hf_geom_1010_off + 338 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_zzz = cbuffer.data(hf_geom_1010_off + 339 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_xxx = cbuffer.data(hf_geom_1010_off + 340 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_xxy = cbuffer.data(hf_geom_1010_off + 341 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_xxz = cbuffer.data(hf_geom_1010_off + 342 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_xyy = cbuffer.data(hf_geom_1010_off + 343 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_xyz = cbuffer.data(hf_geom_1010_off + 344 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_xzz = cbuffer.data(hf_geom_1010_off + 345 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_yyy = cbuffer.data(hf_geom_1010_off + 346 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_yyz = cbuffer.data(hf_geom_1010_off + 347 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_yzz = cbuffer.data(hf_geom_1010_off + 348 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_zzz = cbuffer.data(hf_geom_1010_off + 349 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_xxx = cbuffer.data(hf_geom_1010_off + 350 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_xxy = cbuffer.data(hf_geom_1010_off + 351 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_xxz = cbuffer.data(hf_geom_1010_off + 352 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_xyy = cbuffer.data(hf_geom_1010_off + 353 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_xyz = cbuffer.data(hf_geom_1010_off + 354 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_xzz = cbuffer.data(hf_geom_1010_off + 355 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_yyy = cbuffer.data(hf_geom_1010_off + 356 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_yyz = cbuffer.data(hf_geom_1010_off + 357 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_yzz = cbuffer.data(hf_geom_1010_off + 358 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_zzz = cbuffer.data(hf_geom_1010_off + 359 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_xxx = cbuffer.data(hf_geom_1010_off + 360 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_xxy = cbuffer.data(hf_geom_1010_off + 361 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_xxz = cbuffer.data(hf_geom_1010_off + 362 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_xyy = cbuffer.data(hf_geom_1010_off + 363 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_xyz = cbuffer.data(hf_geom_1010_off + 364 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_xzz = cbuffer.data(hf_geom_1010_off + 365 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_yyy = cbuffer.data(hf_geom_1010_off + 366 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_yyz = cbuffer.data(hf_geom_1010_off + 367 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_yzz = cbuffer.data(hf_geom_1010_off + 368 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_zzz = cbuffer.data(hf_geom_1010_off + 369 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_xxx = cbuffer.data(hf_geom_1010_off + 370 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_xxy = cbuffer.data(hf_geom_1010_off + 371 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_xxz = cbuffer.data(hf_geom_1010_off + 372 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_xyy = cbuffer.data(hf_geom_1010_off + 373 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_xyz = cbuffer.data(hf_geom_1010_off + 374 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_xzz = cbuffer.data(hf_geom_1010_off + 375 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_yyy = cbuffer.data(hf_geom_1010_off + 376 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_yyz = cbuffer.data(hf_geom_1010_off + 377 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_yzz = cbuffer.data(hf_geom_1010_off + 378 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_zzz = cbuffer.data(hf_geom_1010_off + 379 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_xxx = cbuffer.data(hf_geom_1010_off + 380 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_xxy = cbuffer.data(hf_geom_1010_off + 381 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_xxz = cbuffer.data(hf_geom_1010_off + 382 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_xyy = cbuffer.data(hf_geom_1010_off + 383 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_xyz = cbuffer.data(hf_geom_1010_off + 384 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_xzz = cbuffer.data(hf_geom_1010_off + 385 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_yyy = cbuffer.data(hf_geom_1010_off + 386 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_yyz = cbuffer.data(hf_geom_1010_off + 387 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_yzz = cbuffer.data(hf_geom_1010_off + 388 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_zzz = cbuffer.data(hf_geom_1010_off + 389 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_xxx = cbuffer.data(hf_geom_1010_off + 390 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_xxy = cbuffer.data(hf_geom_1010_off + 391 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_xxz = cbuffer.data(hf_geom_1010_off + 392 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_xyy = cbuffer.data(hf_geom_1010_off + 393 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_xyz = cbuffer.data(hf_geom_1010_off + 394 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_xzz = cbuffer.data(hf_geom_1010_off + 395 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_yyy = cbuffer.data(hf_geom_1010_off + 396 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_yyz = cbuffer.data(hf_geom_1010_off + 397 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_yzz = cbuffer.data(hf_geom_1010_off + 398 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_zzz = cbuffer.data(hf_geom_1010_off + 399 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_xxx = cbuffer.data(hf_geom_1010_off + 400 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_xxy = cbuffer.data(hf_geom_1010_off + 401 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_xxz = cbuffer.data(hf_geom_1010_off + 402 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_xyy = cbuffer.data(hf_geom_1010_off + 403 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_xyz = cbuffer.data(hf_geom_1010_off + 404 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_xzz = cbuffer.data(hf_geom_1010_off + 405 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_yyy = cbuffer.data(hf_geom_1010_off + 406 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_yyz = cbuffer.data(hf_geom_1010_off + 407 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_yzz = cbuffer.data(hf_geom_1010_off + 408 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_zzz = cbuffer.data(hf_geom_1010_off + 409 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_xxx = cbuffer.data(hf_geom_1010_off + 410 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_xxy = cbuffer.data(hf_geom_1010_off + 411 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_xxz = cbuffer.data(hf_geom_1010_off + 412 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_xyy = cbuffer.data(hf_geom_1010_off + 413 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_xyz = cbuffer.data(hf_geom_1010_off + 414 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_xzz = cbuffer.data(hf_geom_1010_off + 415 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_yyy = cbuffer.data(hf_geom_1010_off + 416 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_yyz = cbuffer.data(hf_geom_1010_off + 417 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_yzz = cbuffer.data(hf_geom_1010_off + 418 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_zzz = cbuffer.data(hf_geom_1010_off + 419 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_xxx = cbuffer.data(hf_geom_1010_off + 420 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_xxy = cbuffer.data(hf_geom_1010_off + 421 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_xxz = cbuffer.data(hf_geom_1010_off + 422 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_xyy = cbuffer.data(hf_geom_1010_off + 423 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_xyz = cbuffer.data(hf_geom_1010_off + 424 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_xzz = cbuffer.data(hf_geom_1010_off + 425 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_yyy = cbuffer.data(hf_geom_1010_off + 426 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_yyz = cbuffer.data(hf_geom_1010_off + 427 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_yzz = cbuffer.data(hf_geom_1010_off + 428 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_zzz = cbuffer.data(hf_geom_1010_off + 429 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_xxx = cbuffer.data(hf_geom_1010_off + 430 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_xxy = cbuffer.data(hf_geom_1010_off + 431 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_xxz = cbuffer.data(hf_geom_1010_off + 432 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_xyy = cbuffer.data(hf_geom_1010_off + 433 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_xyz = cbuffer.data(hf_geom_1010_off + 434 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_xzz = cbuffer.data(hf_geom_1010_off + 435 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_yyy = cbuffer.data(hf_geom_1010_off + 436 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_yyz = cbuffer.data(hf_geom_1010_off + 437 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_yzz = cbuffer.data(hf_geom_1010_off + 438 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_zzz = cbuffer.data(hf_geom_1010_off + 439 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_xxx = cbuffer.data(hf_geom_1010_off + 440 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_xxy = cbuffer.data(hf_geom_1010_off + 441 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_xxz = cbuffer.data(hf_geom_1010_off + 442 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_xyy = cbuffer.data(hf_geom_1010_off + 443 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_xyz = cbuffer.data(hf_geom_1010_off + 444 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_xzz = cbuffer.data(hf_geom_1010_off + 445 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_yyy = cbuffer.data(hf_geom_1010_off + 446 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_yyz = cbuffer.data(hf_geom_1010_off + 447 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_yzz = cbuffer.data(hf_geom_1010_off + 448 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_zzz = cbuffer.data(hf_geom_1010_off + 449 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_xxx = cbuffer.data(hf_geom_1010_off + 450 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_xxy = cbuffer.data(hf_geom_1010_off + 451 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_xxz = cbuffer.data(hf_geom_1010_off + 452 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_xyy = cbuffer.data(hf_geom_1010_off + 453 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_xyz = cbuffer.data(hf_geom_1010_off + 454 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_xzz = cbuffer.data(hf_geom_1010_off + 455 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_yyy = cbuffer.data(hf_geom_1010_off + 456 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_yyz = cbuffer.data(hf_geom_1010_off + 457 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_yzz = cbuffer.data(hf_geom_1010_off + 458 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_zzz = cbuffer.data(hf_geom_1010_off + 459 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_xxx = cbuffer.data(hf_geom_1010_off + 460 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_xxy = cbuffer.data(hf_geom_1010_off + 461 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_xxz = cbuffer.data(hf_geom_1010_off + 462 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_xyy = cbuffer.data(hf_geom_1010_off + 463 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_xyz = cbuffer.data(hf_geom_1010_off + 464 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_xzz = cbuffer.data(hf_geom_1010_off + 465 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_yyy = cbuffer.data(hf_geom_1010_off + 466 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_yyz = cbuffer.data(hf_geom_1010_off + 467 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_yzz = cbuffer.data(hf_geom_1010_off + 468 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_zzz = cbuffer.data(hf_geom_1010_off + 469 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_xxx = cbuffer.data(hf_geom_1010_off + 470 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_xxy = cbuffer.data(hf_geom_1010_off + 471 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_xxz = cbuffer.data(hf_geom_1010_off + 472 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_xyy = cbuffer.data(hf_geom_1010_off + 473 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_xyz = cbuffer.data(hf_geom_1010_off + 474 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_xzz = cbuffer.data(hf_geom_1010_off + 475 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_yyy = cbuffer.data(hf_geom_1010_off + 476 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_yyz = cbuffer.data(hf_geom_1010_off + 477 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_yzz = cbuffer.data(hf_geom_1010_off + 478 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_zzz = cbuffer.data(hf_geom_1010_off + 479 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_xxx = cbuffer.data(hf_geom_1010_off + 480 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_xxy = cbuffer.data(hf_geom_1010_off + 481 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_xxz = cbuffer.data(hf_geom_1010_off + 482 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_xyy = cbuffer.data(hf_geom_1010_off + 483 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_xyz = cbuffer.data(hf_geom_1010_off + 484 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_xzz = cbuffer.data(hf_geom_1010_off + 485 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_yyy = cbuffer.data(hf_geom_1010_off + 486 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_yyz = cbuffer.data(hf_geom_1010_off + 487 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_yzz = cbuffer.data(hf_geom_1010_off + 488 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_zzz = cbuffer.data(hf_geom_1010_off + 489 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_xxx = cbuffer.data(hf_geom_1010_off + 490 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_xxy = cbuffer.data(hf_geom_1010_off + 491 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_xxz = cbuffer.data(hf_geom_1010_off + 492 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_xyy = cbuffer.data(hf_geom_1010_off + 493 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_xyz = cbuffer.data(hf_geom_1010_off + 494 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_xzz = cbuffer.data(hf_geom_1010_off + 495 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_yyy = cbuffer.data(hf_geom_1010_off + 496 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_yyz = cbuffer.data(hf_geom_1010_off + 497 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_yzz = cbuffer.data(hf_geom_1010_off + 498 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_zzz = cbuffer.data(hf_geom_1010_off + 499 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_xxx = cbuffer.data(hf_geom_1010_off + 500 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_xxy = cbuffer.data(hf_geom_1010_off + 501 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_xxz = cbuffer.data(hf_geom_1010_off + 502 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_xyy = cbuffer.data(hf_geom_1010_off + 503 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_xyz = cbuffer.data(hf_geom_1010_off + 504 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_xzz = cbuffer.data(hf_geom_1010_off + 505 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_yyy = cbuffer.data(hf_geom_1010_off + 506 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_yyz = cbuffer.data(hf_geom_1010_off + 507 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_yzz = cbuffer.data(hf_geom_1010_off + 508 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_zzz = cbuffer.data(hf_geom_1010_off + 509 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_xxx = cbuffer.data(hf_geom_1010_off + 510 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_xxy = cbuffer.data(hf_geom_1010_off + 511 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_xxz = cbuffer.data(hf_geom_1010_off + 512 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_xyy = cbuffer.data(hf_geom_1010_off + 513 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_xyz = cbuffer.data(hf_geom_1010_off + 514 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_xzz = cbuffer.data(hf_geom_1010_off + 515 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_yyy = cbuffer.data(hf_geom_1010_off + 516 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_yyz = cbuffer.data(hf_geom_1010_off + 517 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_yzz = cbuffer.data(hf_geom_1010_off + 518 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_zzz = cbuffer.data(hf_geom_1010_off + 519 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_xxx = cbuffer.data(hf_geom_1010_off + 520 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_xxy = cbuffer.data(hf_geom_1010_off + 521 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_xxz = cbuffer.data(hf_geom_1010_off + 522 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_xyy = cbuffer.data(hf_geom_1010_off + 523 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_xyz = cbuffer.data(hf_geom_1010_off + 524 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_xzz = cbuffer.data(hf_geom_1010_off + 525 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_yyy = cbuffer.data(hf_geom_1010_off + 526 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_yyz = cbuffer.data(hf_geom_1010_off + 527 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_yzz = cbuffer.data(hf_geom_1010_off + 528 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_zzz = cbuffer.data(hf_geom_1010_off + 529 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_xxx = cbuffer.data(hf_geom_1010_off + 530 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_xxy = cbuffer.data(hf_geom_1010_off + 531 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_xxz = cbuffer.data(hf_geom_1010_off + 532 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_xyy = cbuffer.data(hf_geom_1010_off + 533 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_xyz = cbuffer.data(hf_geom_1010_off + 534 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_xzz = cbuffer.data(hf_geom_1010_off + 535 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_yyy = cbuffer.data(hf_geom_1010_off + 536 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_yyz = cbuffer.data(hf_geom_1010_off + 537 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_yzz = cbuffer.data(hf_geom_1010_off + 538 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_zzz = cbuffer.data(hf_geom_1010_off + 539 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_xxx = cbuffer.data(hf_geom_1010_off + 540 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_xxy = cbuffer.data(hf_geom_1010_off + 541 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_xxz = cbuffer.data(hf_geom_1010_off + 542 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_xyy = cbuffer.data(hf_geom_1010_off + 543 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_xyz = cbuffer.data(hf_geom_1010_off + 544 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_xzz = cbuffer.data(hf_geom_1010_off + 545 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_yyy = cbuffer.data(hf_geom_1010_off + 546 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_yyz = cbuffer.data(hf_geom_1010_off + 547 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_yzz = cbuffer.data(hf_geom_1010_off + 548 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_zzz = cbuffer.data(hf_geom_1010_off + 549 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_xxx = cbuffer.data(hf_geom_1010_off + 550 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_xxy = cbuffer.data(hf_geom_1010_off + 551 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_xxz = cbuffer.data(hf_geom_1010_off + 552 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_xyy = cbuffer.data(hf_geom_1010_off + 553 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_xyz = cbuffer.data(hf_geom_1010_off + 554 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_xzz = cbuffer.data(hf_geom_1010_off + 555 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_yyy = cbuffer.data(hf_geom_1010_off + 556 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_yyz = cbuffer.data(hf_geom_1010_off + 557 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_yzz = cbuffer.data(hf_geom_1010_off + 558 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_zzz = cbuffer.data(hf_geom_1010_off + 559 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_xxx = cbuffer.data(hf_geom_1010_off + 560 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_xxy = cbuffer.data(hf_geom_1010_off + 561 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_xxz = cbuffer.data(hf_geom_1010_off + 562 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_xyy = cbuffer.data(hf_geom_1010_off + 563 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_xyz = cbuffer.data(hf_geom_1010_off + 564 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_xzz = cbuffer.data(hf_geom_1010_off + 565 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_yyy = cbuffer.data(hf_geom_1010_off + 566 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_yyz = cbuffer.data(hf_geom_1010_off + 567 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_yzz = cbuffer.data(hf_geom_1010_off + 568 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_zzz = cbuffer.data(hf_geom_1010_off + 569 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_xxx = cbuffer.data(hf_geom_1010_off + 570 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_xxy = cbuffer.data(hf_geom_1010_off + 571 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_xxz = cbuffer.data(hf_geom_1010_off + 572 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_xyy = cbuffer.data(hf_geom_1010_off + 573 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_xyz = cbuffer.data(hf_geom_1010_off + 574 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_xzz = cbuffer.data(hf_geom_1010_off + 575 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_yyy = cbuffer.data(hf_geom_1010_off + 576 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_yyz = cbuffer.data(hf_geom_1010_off + 577 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_yzz = cbuffer.data(hf_geom_1010_off + 578 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_zzz = cbuffer.data(hf_geom_1010_off + 579 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_xxx = cbuffer.data(hf_geom_1010_off + 580 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_xxy = cbuffer.data(hf_geom_1010_off + 581 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_xxz = cbuffer.data(hf_geom_1010_off + 582 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_xyy = cbuffer.data(hf_geom_1010_off + 583 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_xyz = cbuffer.data(hf_geom_1010_off + 584 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_xzz = cbuffer.data(hf_geom_1010_off + 585 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_yyy = cbuffer.data(hf_geom_1010_off + 586 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_yyz = cbuffer.data(hf_geom_1010_off + 587 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_yzz = cbuffer.data(hf_geom_1010_off + 588 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_zzz = cbuffer.data(hf_geom_1010_off + 589 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_xxx = cbuffer.data(hf_geom_1010_off + 590 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_xxy = cbuffer.data(hf_geom_1010_off + 591 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_xxz = cbuffer.data(hf_geom_1010_off + 592 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_xyy = cbuffer.data(hf_geom_1010_off + 593 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_xyz = cbuffer.data(hf_geom_1010_off + 594 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_xzz = cbuffer.data(hf_geom_1010_off + 595 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_yyy = cbuffer.data(hf_geom_1010_off + 596 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_yyz = cbuffer.data(hf_geom_1010_off + 597 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_yzz = cbuffer.data(hf_geom_1010_off + 598 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_zzz = cbuffer.data(hf_geom_1010_off + 599 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_xxx = cbuffer.data(hf_geom_1010_off + 600 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_xxy = cbuffer.data(hf_geom_1010_off + 601 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_xxz = cbuffer.data(hf_geom_1010_off + 602 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_xyy = cbuffer.data(hf_geom_1010_off + 603 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_xyz = cbuffer.data(hf_geom_1010_off + 604 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_xzz = cbuffer.data(hf_geom_1010_off + 605 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_yyy = cbuffer.data(hf_geom_1010_off + 606 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_yyz = cbuffer.data(hf_geom_1010_off + 607 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_yzz = cbuffer.data(hf_geom_1010_off + 608 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_zzz = cbuffer.data(hf_geom_1010_off + 609 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_xxx = cbuffer.data(hf_geom_1010_off + 610 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_xxy = cbuffer.data(hf_geom_1010_off + 611 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_xxz = cbuffer.data(hf_geom_1010_off + 612 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_xyy = cbuffer.data(hf_geom_1010_off + 613 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_xyz = cbuffer.data(hf_geom_1010_off + 614 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_xzz = cbuffer.data(hf_geom_1010_off + 615 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_yyy = cbuffer.data(hf_geom_1010_off + 616 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_yyz = cbuffer.data(hf_geom_1010_off + 617 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_yzz = cbuffer.data(hf_geom_1010_off + 618 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_zzz = cbuffer.data(hf_geom_1010_off + 619 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_xxx = cbuffer.data(hf_geom_1010_off + 620 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_xxy = cbuffer.data(hf_geom_1010_off + 621 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_xxz = cbuffer.data(hf_geom_1010_off + 622 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_xyy = cbuffer.data(hf_geom_1010_off + 623 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_xyz = cbuffer.data(hf_geom_1010_off + 624 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_xzz = cbuffer.data(hf_geom_1010_off + 625 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_yyy = cbuffer.data(hf_geom_1010_off + 626 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_yyz = cbuffer.data(hf_geom_1010_off + 627 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_yzz = cbuffer.data(hf_geom_1010_off + 628 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_zzz = cbuffer.data(hf_geom_1010_off + 629 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_xxx = cbuffer.data(hf_geom_1010_off + 630 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_xxy = cbuffer.data(hf_geom_1010_off + 631 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_xxz = cbuffer.data(hf_geom_1010_off + 632 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_xyy = cbuffer.data(hf_geom_1010_off + 633 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_xyz = cbuffer.data(hf_geom_1010_off + 634 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_xzz = cbuffer.data(hf_geom_1010_off + 635 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_yyy = cbuffer.data(hf_geom_1010_off + 636 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_yyz = cbuffer.data(hf_geom_1010_off + 637 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_yzz = cbuffer.data(hf_geom_1010_off + 638 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_zzz = cbuffer.data(hf_geom_1010_off + 639 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_xxx = cbuffer.data(hf_geom_1010_off + 640 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_xxy = cbuffer.data(hf_geom_1010_off + 641 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_xxz = cbuffer.data(hf_geom_1010_off + 642 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_xyy = cbuffer.data(hf_geom_1010_off + 643 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_xyz = cbuffer.data(hf_geom_1010_off + 644 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_xzz = cbuffer.data(hf_geom_1010_off + 645 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_yyy = cbuffer.data(hf_geom_1010_off + 646 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_yyz = cbuffer.data(hf_geom_1010_off + 647 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_yzz = cbuffer.data(hf_geom_1010_off + 648 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_zzz = cbuffer.data(hf_geom_1010_off + 649 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_xxx = cbuffer.data(hf_geom_1010_off + 650 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_xxy = cbuffer.data(hf_geom_1010_off + 651 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_xxz = cbuffer.data(hf_geom_1010_off + 652 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_xyy = cbuffer.data(hf_geom_1010_off + 653 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_xyz = cbuffer.data(hf_geom_1010_off + 654 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_xzz = cbuffer.data(hf_geom_1010_off + 655 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_yyy = cbuffer.data(hf_geom_1010_off + 656 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_yyz = cbuffer.data(hf_geom_1010_off + 657 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_yzz = cbuffer.data(hf_geom_1010_off + 658 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_zzz = cbuffer.data(hf_geom_1010_off + 659 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_xxx = cbuffer.data(hf_geom_1010_off + 660 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_xxy = cbuffer.data(hf_geom_1010_off + 661 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_xxz = cbuffer.data(hf_geom_1010_off + 662 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_xyy = cbuffer.data(hf_geom_1010_off + 663 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_xyz = cbuffer.data(hf_geom_1010_off + 664 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_xzz = cbuffer.data(hf_geom_1010_off + 665 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_yyy = cbuffer.data(hf_geom_1010_off + 666 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_yyz = cbuffer.data(hf_geom_1010_off + 667 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_yzz = cbuffer.data(hf_geom_1010_off + 668 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_zzz = cbuffer.data(hf_geom_1010_off + 669 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_xxx = cbuffer.data(hf_geom_1010_off + 670 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_xxy = cbuffer.data(hf_geom_1010_off + 671 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_xxz = cbuffer.data(hf_geom_1010_off + 672 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_xyy = cbuffer.data(hf_geom_1010_off + 673 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_xyz = cbuffer.data(hf_geom_1010_off + 674 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_xzz = cbuffer.data(hf_geom_1010_off + 675 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_yyy = cbuffer.data(hf_geom_1010_off + 676 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_yyz = cbuffer.data(hf_geom_1010_off + 677 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_yzz = cbuffer.data(hf_geom_1010_off + 678 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_zzz = cbuffer.data(hf_geom_1010_off + 679 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_xxx = cbuffer.data(hf_geom_1010_off + 680 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_xxy = cbuffer.data(hf_geom_1010_off + 681 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_xxz = cbuffer.data(hf_geom_1010_off + 682 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_xyy = cbuffer.data(hf_geom_1010_off + 683 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_xyz = cbuffer.data(hf_geom_1010_off + 684 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_xzz = cbuffer.data(hf_geom_1010_off + 685 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_yyy = cbuffer.data(hf_geom_1010_off + 686 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_yyz = cbuffer.data(hf_geom_1010_off + 687 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_yzz = cbuffer.data(hf_geom_1010_off + 688 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_zzz = cbuffer.data(hf_geom_1010_off + 689 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_xxx = cbuffer.data(hf_geom_1010_off + 690 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_xxy = cbuffer.data(hf_geom_1010_off + 691 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_xxz = cbuffer.data(hf_geom_1010_off + 692 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_xyy = cbuffer.data(hf_geom_1010_off + 693 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_xyz = cbuffer.data(hf_geom_1010_off + 694 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_xzz = cbuffer.data(hf_geom_1010_off + 695 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_yyy = cbuffer.data(hf_geom_1010_off + 696 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_yyz = cbuffer.data(hf_geom_1010_off + 697 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_yzz = cbuffer.data(hf_geom_1010_off + 698 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_zzz = cbuffer.data(hf_geom_1010_off + 699 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_xxx = cbuffer.data(hf_geom_1010_off + 700 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_xxy = cbuffer.data(hf_geom_1010_off + 701 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_xxz = cbuffer.data(hf_geom_1010_off + 702 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_xyy = cbuffer.data(hf_geom_1010_off + 703 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_xyz = cbuffer.data(hf_geom_1010_off + 704 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_xzz = cbuffer.data(hf_geom_1010_off + 705 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_yyy = cbuffer.data(hf_geom_1010_off + 706 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_yyz = cbuffer.data(hf_geom_1010_off + 707 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_yzz = cbuffer.data(hf_geom_1010_off + 708 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_zzz = cbuffer.data(hf_geom_1010_off + 709 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_xxx = cbuffer.data(hf_geom_1010_off + 710 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_xxy = cbuffer.data(hf_geom_1010_off + 711 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_xxz = cbuffer.data(hf_geom_1010_off + 712 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_xyy = cbuffer.data(hf_geom_1010_off + 713 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_xyz = cbuffer.data(hf_geom_1010_off + 714 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_xzz = cbuffer.data(hf_geom_1010_off + 715 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_yyy = cbuffer.data(hf_geom_1010_off + 716 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_yyz = cbuffer.data(hf_geom_1010_off + 717 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_yzz = cbuffer.data(hf_geom_1010_off + 718 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_zzz = cbuffer.data(hf_geom_1010_off + 719 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_xxx = cbuffer.data(hf_geom_1010_off + 720 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_xxy = cbuffer.data(hf_geom_1010_off + 721 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_xxz = cbuffer.data(hf_geom_1010_off + 722 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_xyy = cbuffer.data(hf_geom_1010_off + 723 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_xyz = cbuffer.data(hf_geom_1010_off + 724 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_xzz = cbuffer.data(hf_geom_1010_off + 725 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_yyy = cbuffer.data(hf_geom_1010_off + 726 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_yyz = cbuffer.data(hf_geom_1010_off + 727 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_yzz = cbuffer.data(hf_geom_1010_off + 728 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_zzz = cbuffer.data(hf_geom_1010_off + 729 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_xxx = cbuffer.data(hf_geom_1010_off + 730 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_xxy = cbuffer.data(hf_geom_1010_off + 731 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_xxz = cbuffer.data(hf_geom_1010_off + 732 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_xyy = cbuffer.data(hf_geom_1010_off + 733 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_xyz = cbuffer.data(hf_geom_1010_off + 734 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_xzz = cbuffer.data(hf_geom_1010_off + 735 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_yyy = cbuffer.data(hf_geom_1010_off + 736 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_yyz = cbuffer.data(hf_geom_1010_off + 737 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_yzz = cbuffer.data(hf_geom_1010_off + 738 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_zzz = cbuffer.data(hf_geom_1010_off + 739 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_xxx = cbuffer.data(hf_geom_1010_off + 740 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_xxy = cbuffer.data(hf_geom_1010_off + 741 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_xxz = cbuffer.data(hf_geom_1010_off + 742 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_xyy = cbuffer.data(hf_geom_1010_off + 743 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_xyz = cbuffer.data(hf_geom_1010_off + 744 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_xzz = cbuffer.data(hf_geom_1010_off + 745 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_yyy = cbuffer.data(hf_geom_1010_off + 746 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_yyz = cbuffer.data(hf_geom_1010_off + 747 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_yzz = cbuffer.data(hf_geom_1010_off + 748 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_zzz = cbuffer.data(hf_geom_1010_off + 749 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_xxx = cbuffer.data(hf_geom_1010_off + 750 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_xxy = cbuffer.data(hf_geom_1010_off + 751 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_xxz = cbuffer.data(hf_geom_1010_off + 752 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_xyy = cbuffer.data(hf_geom_1010_off + 753 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_xyz = cbuffer.data(hf_geom_1010_off + 754 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_xzz = cbuffer.data(hf_geom_1010_off + 755 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_yyy = cbuffer.data(hf_geom_1010_off + 756 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_yyz = cbuffer.data(hf_geom_1010_off + 757 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_yzz = cbuffer.data(hf_geom_1010_off + 758 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_zzz = cbuffer.data(hf_geom_1010_off + 759 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_xxx = cbuffer.data(hf_geom_1010_off + 760 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_xxy = cbuffer.data(hf_geom_1010_off + 761 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_xxz = cbuffer.data(hf_geom_1010_off + 762 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_xyy = cbuffer.data(hf_geom_1010_off + 763 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_xyz = cbuffer.data(hf_geom_1010_off + 764 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_xzz = cbuffer.data(hf_geom_1010_off + 765 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_yyy = cbuffer.data(hf_geom_1010_off + 766 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_yyz = cbuffer.data(hf_geom_1010_off + 767 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_yzz = cbuffer.data(hf_geom_1010_off + 768 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_zzz = cbuffer.data(hf_geom_1010_off + 769 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_xxx = cbuffer.data(hf_geom_1010_off + 770 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_xxy = cbuffer.data(hf_geom_1010_off + 771 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_xxz = cbuffer.data(hf_geom_1010_off + 772 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_xyy = cbuffer.data(hf_geom_1010_off + 773 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_xyz = cbuffer.data(hf_geom_1010_off + 774 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_xzz = cbuffer.data(hf_geom_1010_off + 775 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_yyy = cbuffer.data(hf_geom_1010_off + 776 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_yyz = cbuffer.data(hf_geom_1010_off + 777 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_yzz = cbuffer.data(hf_geom_1010_off + 778 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_zzz = cbuffer.data(hf_geom_1010_off + 779 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_xxx = cbuffer.data(hf_geom_1010_off + 780 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_xxy = cbuffer.data(hf_geom_1010_off + 781 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_xxz = cbuffer.data(hf_geom_1010_off + 782 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_xyy = cbuffer.data(hf_geom_1010_off + 783 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_xyz = cbuffer.data(hf_geom_1010_off + 784 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_xzz = cbuffer.data(hf_geom_1010_off + 785 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_yyy = cbuffer.data(hf_geom_1010_off + 786 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_yyz = cbuffer.data(hf_geom_1010_off + 787 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_yzz = cbuffer.data(hf_geom_1010_off + 788 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_zzz = cbuffer.data(hf_geom_1010_off + 789 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_xxx = cbuffer.data(hf_geom_1010_off + 790 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_xxy = cbuffer.data(hf_geom_1010_off + 791 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_xxz = cbuffer.data(hf_geom_1010_off + 792 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_xyy = cbuffer.data(hf_geom_1010_off + 793 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_xyz = cbuffer.data(hf_geom_1010_off + 794 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_xzz = cbuffer.data(hf_geom_1010_off + 795 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_yyy = cbuffer.data(hf_geom_1010_off + 796 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_yyz = cbuffer.data(hf_geom_1010_off + 797 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_yzz = cbuffer.data(hf_geom_1010_off + 798 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_zzz = cbuffer.data(hf_geom_1010_off + 799 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_xxx = cbuffer.data(hf_geom_1010_off + 800 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_xxy = cbuffer.data(hf_geom_1010_off + 801 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_xxz = cbuffer.data(hf_geom_1010_off + 802 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_xyy = cbuffer.data(hf_geom_1010_off + 803 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_xyz = cbuffer.data(hf_geom_1010_off + 804 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_xzz = cbuffer.data(hf_geom_1010_off + 805 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_yyy = cbuffer.data(hf_geom_1010_off + 806 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_yyz = cbuffer.data(hf_geom_1010_off + 807 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_yzz = cbuffer.data(hf_geom_1010_off + 808 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_zzz = cbuffer.data(hf_geom_1010_off + 809 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_xxx = cbuffer.data(hf_geom_1010_off + 810 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_xxy = cbuffer.data(hf_geom_1010_off + 811 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_xxz = cbuffer.data(hf_geom_1010_off + 812 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_xyy = cbuffer.data(hf_geom_1010_off + 813 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_xyz = cbuffer.data(hf_geom_1010_off + 814 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_xzz = cbuffer.data(hf_geom_1010_off + 815 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_yyy = cbuffer.data(hf_geom_1010_off + 816 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_yyz = cbuffer.data(hf_geom_1010_off + 817 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_yzz = cbuffer.data(hf_geom_1010_off + 818 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_zzz = cbuffer.data(hf_geom_1010_off + 819 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_xxx = cbuffer.data(hf_geom_1010_off + 820 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_xxy = cbuffer.data(hf_geom_1010_off + 821 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_xxz = cbuffer.data(hf_geom_1010_off + 822 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_xyy = cbuffer.data(hf_geom_1010_off + 823 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_xyz = cbuffer.data(hf_geom_1010_off + 824 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_xzz = cbuffer.data(hf_geom_1010_off + 825 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_yyy = cbuffer.data(hf_geom_1010_off + 826 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_yyz = cbuffer.data(hf_geom_1010_off + 827 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_yzz = cbuffer.data(hf_geom_1010_off + 828 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_zzz = cbuffer.data(hf_geom_1010_off + 829 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_xxx = cbuffer.data(hf_geom_1010_off + 830 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_xxy = cbuffer.data(hf_geom_1010_off + 831 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_xxz = cbuffer.data(hf_geom_1010_off + 832 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_xyy = cbuffer.data(hf_geom_1010_off + 833 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_xyz = cbuffer.data(hf_geom_1010_off + 834 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_xzz = cbuffer.data(hf_geom_1010_off + 835 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_yyy = cbuffer.data(hf_geom_1010_off + 836 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_yyz = cbuffer.data(hf_geom_1010_off + 837 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_yzz = cbuffer.data(hf_geom_1010_off + 838 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_zzz = cbuffer.data(hf_geom_1010_off + 839 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_xxx = cbuffer.data(hf_geom_1010_off + 840 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_xxy = cbuffer.data(hf_geom_1010_off + 841 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_xxz = cbuffer.data(hf_geom_1010_off + 842 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_xyy = cbuffer.data(hf_geom_1010_off + 843 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_xyz = cbuffer.data(hf_geom_1010_off + 844 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_xzz = cbuffer.data(hf_geom_1010_off + 845 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_yyy = cbuffer.data(hf_geom_1010_off + 846 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_yyz = cbuffer.data(hf_geom_1010_off + 847 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_yzz = cbuffer.data(hf_geom_1010_off + 848 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_zzz = cbuffer.data(hf_geom_1010_off + 849 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_xxx = cbuffer.data(hf_geom_1010_off + 850 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_xxy = cbuffer.data(hf_geom_1010_off + 851 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_xxz = cbuffer.data(hf_geom_1010_off + 852 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_xyy = cbuffer.data(hf_geom_1010_off + 853 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_xyz = cbuffer.data(hf_geom_1010_off + 854 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_xzz = cbuffer.data(hf_geom_1010_off + 855 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_yyy = cbuffer.data(hf_geom_1010_off + 856 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_yyz = cbuffer.data(hf_geom_1010_off + 857 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_yzz = cbuffer.data(hf_geom_1010_off + 858 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_zzz = cbuffer.data(hf_geom_1010_off + 859 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_xxx = cbuffer.data(hf_geom_1010_off + 860 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_xxy = cbuffer.data(hf_geom_1010_off + 861 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_xxz = cbuffer.data(hf_geom_1010_off + 862 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_xyy = cbuffer.data(hf_geom_1010_off + 863 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_xyz = cbuffer.data(hf_geom_1010_off + 864 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_xzz = cbuffer.data(hf_geom_1010_off + 865 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_yyy = cbuffer.data(hf_geom_1010_off + 866 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_yyz = cbuffer.data(hf_geom_1010_off + 867 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_yzz = cbuffer.data(hf_geom_1010_off + 868 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_zzz = cbuffer.data(hf_geom_1010_off + 869 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_xxx = cbuffer.data(hf_geom_1010_off + 870 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_xxy = cbuffer.data(hf_geom_1010_off + 871 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_xxz = cbuffer.data(hf_geom_1010_off + 872 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_xyy = cbuffer.data(hf_geom_1010_off + 873 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_xyz = cbuffer.data(hf_geom_1010_off + 874 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_xzz = cbuffer.data(hf_geom_1010_off + 875 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_yyy = cbuffer.data(hf_geom_1010_off + 876 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_yyz = cbuffer.data(hf_geom_1010_off + 877 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_yzz = cbuffer.data(hf_geom_1010_off + 878 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_zzz = cbuffer.data(hf_geom_1010_off + 879 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_xxx = cbuffer.data(hf_geom_1010_off + 880 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_xxy = cbuffer.data(hf_geom_1010_off + 881 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_xxz = cbuffer.data(hf_geom_1010_off + 882 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_xyy = cbuffer.data(hf_geom_1010_off + 883 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_xyz = cbuffer.data(hf_geom_1010_off + 884 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_xzz = cbuffer.data(hf_geom_1010_off + 885 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_yyy = cbuffer.data(hf_geom_1010_off + 886 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_yyz = cbuffer.data(hf_geom_1010_off + 887 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_yzz = cbuffer.data(hf_geom_1010_off + 888 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_zzz = cbuffer.data(hf_geom_1010_off + 889 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_xxx = cbuffer.data(hf_geom_1010_off + 890 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_xxy = cbuffer.data(hf_geom_1010_off + 891 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_xxz = cbuffer.data(hf_geom_1010_off + 892 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_xyy = cbuffer.data(hf_geom_1010_off + 893 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_xyz = cbuffer.data(hf_geom_1010_off + 894 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_xzz = cbuffer.data(hf_geom_1010_off + 895 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_yyy = cbuffer.data(hf_geom_1010_off + 896 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_yyz = cbuffer.data(hf_geom_1010_off + 897 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_yzz = cbuffer.data(hf_geom_1010_off + 898 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_zzz = cbuffer.data(hf_geom_1010_off + 899 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_xxx = cbuffer.data(hf_geom_1010_off + 900 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_xxy = cbuffer.data(hf_geom_1010_off + 901 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_xxz = cbuffer.data(hf_geom_1010_off + 902 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_xyy = cbuffer.data(hf_geom_1010_off + 903 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_xyz = cbuffer.data(hf_geom_1010_off + 904 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_xzz = cbuffer.data(hf_geom_1010_off + 905 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_yyy = cbuffer.data(hf_geom_1010_off + 906 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_yyz = cbuffer.data(hf_geom_1010_off + 907 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_yzz = cbuffer.data(hf_geom_1010_off + 908 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_zzz = cbuffer.data(hf_geom_1010_off + 909 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_xxx = cbuffer.data(hf_geom_1010_off + 910 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_xxy = cbuffer.data(hf_geom_1010_off + 911 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_xxz = cbuffer.data(hf_geom_1010_off + 912 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_xyy = cbuffer.data(hf_geom_1010_off + 913 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_xyz = cbuffer.data(hf_geom_1010_off + 914 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_xzz = cbuffer.data(hf_geom_1010_off + 915 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_yyy = cbuffer.data(hf_geom_1010_off + 916 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_yyz = cbuffer.data(hf_geom_1010_off + 917 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_yzz = cbuffer.data(hf_geom_1010_off + 918 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_zzz = cbuffer.data(hf_geom_1010_off + 919 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_xxx = cbuffer.data(hf_geom_1010_off + 920 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_xxy = cbuffer.data(hf_geom_1010_off + 921 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_xxz = cbuffer.data(hf_geom_1010_off + 922 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_xyy = cbuffer.data(hf_geom_1010_off + 923 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_xyz = cbuffer.data(hf_geom_1010_off + 924 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_xzz = cbuffer.data(hf_geom_1010_off + 925 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_yyy = cbuffer.data(hf_geom_1010_off + 926 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_yyz = cbuffer.data(hf_geom_1010_off + 927 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_yzz = cbuffer.data(hf_geom_1010_off + 928 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_zzz = cbuffer.data(hf_geom_1010_off + 929 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_xxx = cbuffer.data(hf_geom_1010_off + 930 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_xxy = cbuffer.data(hf_geom_1010_off + 931 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_xxz = cbuffer.data(hf_geom_1010_off + 932 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_xyy = cbuffer.data(hf_geom_1010_off + 933 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_xyz = cbuffer.data(hf_geom_1010_off + 934 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_xzz = cbuffer.data(hf_geom_1010_off + 935 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_yyy = cbuffer.data(hf_geom_1010_off + 936 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_yyz = cbuffer.data(hf_geom_1010_off + 937 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_yzz = cbuffer.data(hf_geom_1010_off + 938 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_zzz = cbuffer.data(hf_geom_1010_off + 939 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_xxx = cbuffer.data(hf_geom_1010_off + 940 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_xxy = cbuffer.data(hf_geom_1010_off + 941 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_xxz = cbuffer.data(hf_geom_1010_off + 942 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_xyy = cbuffer.data(hf_geom_1010_off + 943 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_xyz = cbuffer.data(hf_geom_1010_off + 944 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_xzz = cbuffer.data(hf_geom_1010_off + 945 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_yyy = cbuffer.data(hf_geom_1010_off + 946 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_yyz = cbuffer.data(hf_geom_1010_off + 947 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_yzz = cbuffer.data(hf_geom_1010_off + 948 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_zzz = cbuffer.data(hf_geom_1010_off + 949 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_xxx = cbuffer.data(hf_geom_1010_off + 950 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_xxy = cbuffer.data(hf_geom_1010_off + 951 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_xxz = cbuffer.data(hf_geom_1010_off + 952 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_xyy = cbuffer.data(hf_geom_1010_off + 953 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_xyz = cbuffer.data(hf_geom_1010_off + 954 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_xzz = cbuffer.data(hf_geom_1010_off + 955 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_yyy = cbuffer.data(hf_geom_1010_off + 956 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_yyz = cbuffer.data(hf_geom_1010_off + 957 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_yzz = cbuffer.data(hf_geom_1010_off + 958 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_zzz = cbuffer.data(hf_geom_1010_off + 959 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_xxx = cbuffer.data(hf_geom_1010_off + 960 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_xxy = cbuffer.data(hf_geom_1010_off + 961 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_xxz = cbuffer.data(hf_geom_1010_off + 962 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_xyy = cbuffer.data(hf_geom_1010_off + 963 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_xyz = cbuffer.data(hf_geom_1010_off + 964 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_xzz = cbuffer.data(hf_geom_1010_off + 965 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_yyy = cbuffer.data(hf_geom_1010_off + 966 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_yyz = cbuffer.data(hf_geom_1010_off + 967 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_yzz = cbuffer.data(hf_geom_1010_off + 968 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_zzz = cbuffer.data(hf_geom_1010_off + 969 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_xxx = cbuffer.data(hf_geom_1010_off + 970 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_xxy = cbuffer.data(hf_geom_1010_off + 971 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_xxz = cbuffer.data(hf_geom_1010_off + 972 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_xyy = cbuffer.data(hf_geom_1010_off + 973 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_xyz = cbuffer.data(hf_geom_1010_off + 974 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_xzz = cbuffer.data(hf_geom_1010_off + 975 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_yyy = cbuffer.data(hf_geom_1010_off + 976 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_yyz = cbuffer.data(hf_geom_1010_off + 977 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_yzz = cbuffer.data(hf_geom_1010_off + 978 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_zzz = cbuffer.data(hf_geom_1010_off + 979 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_xxx = cbuffer.data(hf_geom_1010_off + 980 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_xxy = cbuffer.data(hf_geom_1010_off + 981 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_xxz = cbuffer.data(hf_geom_1010_off + 982 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_xyy = cbuffer.data(hf_geom_1010_off + 983 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_xyz = cbuffer.data(hf_geom_1010_off + 984 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_xzz = cbuffer.data(hf_geom_1010_off + 985 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_yyy = cbuffer.data(hf_geom_1010_off + 986 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_yyz = cbuffer.data(hf_geom_1010_off + 987 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_yzz = cbuffer.data(hf_geom_1010_off + 988 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_zzz = cbuffer.data(hf_geom_1010_off + 989 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_xxx = cbuffer.data(hf_geom_1010_off + 990 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_xxy = cbuffer.data(hf_geom_1010_off + 991 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_xxz = cbuffer.data(hf_geom_1010_off + 992 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_xyy = cbuffer.data(hf_geom_1010_off + 993 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_xyz = cbuffer.data(hf_geom_1010_off + 994 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_xzz = cbuffer.data(hf_geom_1010_off + 995 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_yyy = cbuffer.data(hf_geom_1010_off + 996 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_yyz = cbuffer.data(hf_geom_1010_off + 997 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_yzz = cbuffer.data(hf_geom_1010_off + 998 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_zzz = cbuffer.data(hf_geom_1010_off + 999 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_xxx = cbuffer.data(hf_geom_1010_off + 1000 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_xxy = cbuffer.data(hf_geom_1010_off + 1001 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_xxz = cbuffer.data(hf_geom_1010_off + 1002 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_xyy = cbuffer.data(hf_geom_1010_off + 1003 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_xyz = cbuffer.data(hf_geom_1010_off + 1004 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_xzz = cbuffer.data(hf_geom_1010_off + 1005 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_yyy = cbuffer.data(hf_geom_1010_off + 1006 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_yyz = cbuffer.data(hf_geom_1010_off + 1007 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_yzz = cbuffer.data(hf_geom_1010_off + 1008 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_zzz = cbuffer.data(hf_geom_1010_off + 1009 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_xxx = cbuffer.data(hf_geom_1010_off + 1010 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_xxy = cbuffer.data(hf_geom_1010_off + 1011 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_xxz = cbuffer.data(hf_geom_1010_off + 1012 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_xyy = cbuffer.data(hf_geom_1010_off + 1013 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_xyz = cbuffer.data(hf_geom_1010_off + 1014 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_xzz = cbuffer.data(hf_geom_1010_off + 1015 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_yyy = cbuffer.data(hf_geom_1010_off + 1016 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_yyz = cbuffer.data(hf_geom_1010_off + 1017 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_yzz = cbuffer.data(hf_geom_1010_off + 1018 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_zzz = cbuffer.data(hf_geom_1010_off + 1019 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_xxx = cbuffer.data(hf_geom_1010_off + 1020 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_xxy = cbuffer.data(hf_geom_1010_off + 1021 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_xxz = cbuffer.data(hf_geom_1010_off + 1022 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_xyy = cbuffer.data(hf_geom_1010_off + 1023 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_xyz = cbuffer.data(hf_geom_1010_off + 1024 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_xzz = cbuffer.data(hf_geom_1010_off + 1025 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_yyy = cbuffer.data(hf_geom_1010_off + 1026 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_yyz = cbuffer.data(hf_geom_1010_off + 1027 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_yzz = cbuffer.data(hf_geom_1010_off + 1028 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_zzz = cbuffer.data(hf_geom_1010_off + 1029 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1030 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1031 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1032 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1033 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1034 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1035 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1036 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1037 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1038 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1039 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1040 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1041 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1042 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1043 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1044 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1045 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1046 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1047 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1048 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1049 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_xxx = cbuffer.data(hf_geom_1010_off + 1050 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_xxy = cbuffer.data(hf_geom_1010_off + 1051 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_xxz = cbuffer.data(hf_geom_1010_off + 1052 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_xyy = cbuffer.data(hf_geom_1010_off + 1053 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_xyz = cbuffer.data(hf_geom_1010_off + 1054 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_xzz = cbuffer.data(hf_geom_1010_off + 1055 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_yyy = cbuffer.data(hf_geom_1010_off + 1056 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_yyz = cbuffer.data(hf_geom_1010_off + 1057 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_yzz = cbuffer.data(hf_geom_1010_off + 1058 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_zzz = cbuffer.data(hf_geom_1010_off + 1059 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_xxx = cbuffer.data(hf_geom_1010_off + 1060 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_xxy = cbuffer.data(hf_geom_1010_off + 1061 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_xxz = cbuffer.data(hf_geom_1010_off + 1062 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_xyy = cbuffer.data(hf_geom_1010_off + 1063 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_xyz = cbuffer.data(hf_geom_1010_off + 1064 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_xzz = cbuffer.data(hf_geom_1010_off + 1065 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_yyy = cbuffer.data(hf_geom_1010_off + 1066 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_yyz = cbuffer.data(hf_geom_1010_off + 1067 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_yzz = cbuffer.data(hf_geom_1010_off + 1068 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_zzz = cbuffer.data(hf_geom_1010_off + 1069 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_xxx = cbuffer.data(hf_geom_1010_off + 1070 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_xxy = cbuffer.data(hf_geom_1010_off + 1071 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_xxz = cbuffer.data(hf_geom_1010_off + 1072 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_xyy = cbuffer.data(hf_geom_1010_off + 1073 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_xyz = cbuffer.data(hf_geom_1010_off + 1074 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_xzz = cbuffer.data(hf_geom_1010_off + 1075 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_yyy = cbuffer.data(hf_geom_1010_off + 1076 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_yyz = cbuffer.data(hf_geom_1010_off + 1077 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_yzz = cbuffer.data(hf_geom_1010_off + 1078 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_zzz = cbuffer.data(hf_geom_1010_off + 1079 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_xxx = cbuffer.data(hf_geom_1010_off + 1080 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_xxy = cbuffer.data(hf_geom_1010_off + 1081 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_xxz = cbuffer.data(hf_geom_1010_off + 1082 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_xyy = cbuffer.data(hf_geom_1010_off + 1083 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_xyz = cbuffer.data(hf_geom_1010_off + 1084 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_xzz = cbuffer.data(hf_geom_1010_off + 1085 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_yyy = cbuffer.data(hf_geom_1010_off + 1086 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_yyz = cbuffer.data(hf_geom_1010_off + 1087 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_yzz = cbuffer.data(hf_geom_1010_off + 1088 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_zzz = cbuffer.data(hf_geom_1010_off + 1089 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_xxx = cbuffer.data(hf_geom_1010_off + 1090 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_xxy = cbuffer.data(hf_geom_1010_off + 1091 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_xxz = cbuffer.data(hf_geom_1010_off + 1092 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_xyy = cbuffer.data(hf_geom_1010_off + 1093 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_xyz = cbuffer.data(hf_geom_1010_off + 1094 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_xzz = cbuffer.data(hf_geom_1010_off + 1095 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_yyy = cbuffer.data(hf_geom_1010_off + 1096 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_yyz = cbuffer.data(hf_geom_1010_off + 1097 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_yzz = cbuffer.data(hf_geom_1010_off + 1098 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_zzz = cbuffer.data(hf_geom_1010_off + 1099 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_xxx = cbuffer.data(hf_geom_1010_off + 1100 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_xxy = cbuffer.data(hf_geom_1010_off + 1101 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_xxz = cbuffer.data(hf_geom_1010_off + 1102 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_xyy = cbuffer.data(hf_geom_1010_off + 1103 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_xyz = cbuffer.data(hf_geom_1010_off + 1104 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_xzz = cbuffer.data(hf_geom_1010_off + 1105 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_yyy = cbuffer.data(hf_geom_1010_off + 1106 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_yyz = cbuffer.data(hf_geom_1010_off + 1107 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_yzz = cbuffer.data(hf_geom_1010_off + 1108 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_zzz = cbuffer.data(hf_geom_1010_off + 1109 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_xxx = cbuffer.data(hf_geom_1010_off + 1110 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_xxy = cbuffer.data(hf_geom_1010_off + 1111 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_xxz = cbuffer.data(hf_geom_1010_off + 1112 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_xyy = cbuffer.data(hf_geom_1010_off + 1113 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_xyz = cbuffer.data(hf_geom_1010_off + 1114 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_xzz = cbuffer.data(hf_geom_1010_off + 1115 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_yyy = cbuffer.data(hf_geom_1010_off + 1116 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_yyz = cbuffer.data(hf_geom_1010_off + 1117 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_yzz = cbuffer.data(hf_geom_1010_off + 1118 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_zzz = cbuffer.data(hf_geom_1010_off + 1119 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_xxx = cbuffer.data(hf_geom_1010_off + 1120 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_xxy = cbuffer.data(hf_geom_1010_off + 1121 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_xxz = cbuffer.data(hf_geom_1010_off + 1122 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_xyy = cbuffer.data(hf_geom_1010_off + 1123 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_xyz = cbuffer.data(hf_geom_1010_off + 1124 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_xzz = cbuffer.data(hf_geom_1010_off + 1125 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_yyy = cbuffer.data(hf_geom_1010_off + 1126 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_yyz = cbuffer.data(hf_geom_1010_off + 1127 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_yzz = cbuffer.data(hf_geom_1010_off + 1128 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_zzz = cbuffer.data(hf_geom_1010_off + 1129 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_xxx = cbuffer.data(hf_geom_1010_off + 1130 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_xxy = cbuffer.data(hf_geom_1010_off + 1131 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_xxz = cbuffer.data(hf_geom_1010_off + 1132 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_xyy = cbuffer.data(hf_geom_1010_off + 1133 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_xyz = cbuffer.data(hf_geom_1010_off + 1134 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_xzz = cbuffer.data(hf_geom_1010_off + 1135 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_yyy = cbuffer.data(hf_geom_1010_off + 1136 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_yyz = cbuffer.data(hf_geom_1010_off + 1137 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_yzz = cbuffer.data(hf_geom_1010_off + 1138 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_zzz = cbuffer.data(hf_geom_1010_off + 1139 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_xxx = cbuffer.data(hf_geom_1010_off + 1140 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_xxy = cbuffer.data(hf_geom_1010_off + 1141 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_xxz = cbuffer.data(hf_geom_1010_off + 1142 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_xyy = cbuffer.data(hf_geom_1010_off + 1143 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_xyz = cbuffer.data(hf_geom_1010_off + 1144 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_xzz = cbuffer.data(hf_geom_1010_off + 1145 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_yyy = cbuffer.data(hf_geom_1010_off + 1146 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_yyz = cbuffer.data(hf_geom_1010_off + 1147 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_yzz = cbuffer.data(hf_geom_1010_off + 1148 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_zzz = cbuffer.data(hf_geom_1010_off + 1149 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_xxx = cbuffer.data(hf_geom_1010_off + 1150 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_xxy = cbuffer.data(hf_geom_1010_off + 1151 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_xxz = cbuffer.data(hf_geom_1010_off + 1152 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_xyy = cbuffer.data(hf_geom_1010_off + 1153 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_xyz = cbuffer.data(hf_geom_1010_off + 1154 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_xzz = cbuffer.data(hf_geom_1010_off + 1155 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_yyy = cbuffer.data(hf_geom_1010_off + 1156 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_yyz = cbuffer.data(hf_geom_1010_off + 1157 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_yzz = cbuffer.data(hf_geom_1010_off + 1158 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_zzz = cbuffer.data(hf_geom_1010_off + 1159 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_xxx = cbuffer.data(hf_geom_1010_off + 1160 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_xxy = cbuffer.data(hf_geom_1010_off + 1161 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_xxz = cbuffer.data(hf_geom_1010_off + 1162 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_xyy = cbuffer.data(hf_geom_1010_off + 1163 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_xyz = cbuffer.data(hf_geom_1010_off + 1164 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_xzz = cbuffer.data(hf_geom_1010_off + 1165 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_yyy = cbuffer.data(hf_geom_1010_off + 1166 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_yyz = cbuffer.data(hf_geom_1010_off + 1167 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_yzz = cbuffer.data(hf_geom_1010_off + 1168 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_zzz = cbuffer.data(hf_geom_1010_off + 1169 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_xxx = cbuffer.data(hf_geom_1010_off + 1170 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_xxy = cbuffer.data(hf_geom_1010_off + 1171 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_xxz = cbuffer.data(hf_geom_1010_off + 1172 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_xyy = cbuffer.data(hf_geom_1010_off + 1173 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_xyz = cbuffer.data(hf_geom_1010_off + 1174 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_xzz = cbuffer.data(hf_geom_1010_off + 1175 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_yyy = cbuffer.data(hf_geom_1010_off + 1176 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_yyz = cbuffer.data(hf_geom_1010_off + 1177 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_yzz = cbuffer.data(hf_geom_1010_off + 1178 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_zzz = cbuffer.data(hf_geom_1010_off + 1179 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_xxx = cbuffer.data(hf_geom_1010_off + 1180 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_xxy = cbuffer.data(hf_geom_1010_off + 1181 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_xxz = cbuffer.data(hf_geom_1010_off + 1182 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_xyy = cbuffer.data(hf_geom_1010_off + 1183 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_xyz = cbuffer.data(hf_geom_1010_off + 1184 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_xzz = cbuffer.data(hf_geom_1010_off + 1185 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_yyy = cbuffer.data(hf_geom_1010_off + 1186 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_yyz = cbuffer.data(hf_geom_1010_off + 1187 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_yzz = cbuffer.data(hf_geom_1010_off + 1188 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_zzz = cbuffer.data(hf_geom_1010_off + 1189 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1190 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1191 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1192 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1193 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1194 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1195 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1196 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1197 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1198 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1199 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_xxx = cbuffer.data(hf_geom_1010_off + 1200 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_xxy = cbuffer.data(hf_geom_1010_off + 1201 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_xxz = cbuffer.data(hf_geom_1010_off + 1202 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_xyy = cbuffer.data(hf_geom_1010_off + 1203 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_xyz = cbuffer.data(hf_geom_1010_off + 1204 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_xzz = cbuffer.data(hf_geom_1010_off + 1205 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_yyy = cbuffer.data(hf_geom_1010_off + 1206 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_yyz = cbuffer.data(hf_geom_1010_off + 1207 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_yzz = cbuffer.data(hf_geom_1010_off + 1208 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_zzz = cbuffer.data(hf_geom_1010_off + 1209 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_xxx = cbuffer.data(hf_geom_1010_off + 1210 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_xxy = cbuffer.data(hf_geom_1010_off + 1211 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_xxz = cbuffer.data(hf_geom_1010_off + 1212 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_xyy = cbuffer.data(hf_geom_1010_off + 1213 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_xyz = cbuffer.data(hf_geom_1010_off + 1214 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_xzz = cbuffer.data(hf_geom_1010_off + 1215 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_yyy = cbuffer.data(hf_geom_1010_off + 1216 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_yyz = cbuffer.data(hf_geom_1010_off + 1217 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_yzz = cbuffer.data(hf_geom_1010_off + 1218 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_zzz = cbuffer.data(hf_geom_1010_off + 1219 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_xxx = cbuffer.data(hf_geom_1010_off + 1220 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_xxy = cbuffer.data(hf_geom_1010_off + 1221 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_xxz = cbuffer.data(hf_geom_1010_off + 1222 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_xyy = cbuffer.data(hf_geom_1010_off + 1223 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_xyz = cbuffer.data(hf_geom_1010_off + 1224 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_xzz = cbuffer.data(hf_geom_1010_off + 1225 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_yyy = cbuffer.data(hf_geom_1010_off + 1226 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_yyz = cbuffer.data(hf_geom_1010_off + 1227 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_yzz = cbuffer.data(hf_geom_1010_off + 1228 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_zzz = cbuffer.data(hf_geom_1010_off + 1229 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_xxx = cbuffer.data(hf_geom_1010_off + 1230 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_xxy = cbuffer.data(hf_geom_1010_off + 1231 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_xxz = cbuffer.data(hf_geom_1010_off + 1232 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_xyy = cbuffer.data(hf_geom_1010_off + 1233 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_xyz = cbuffer.data(hf_geom_1010_off + 1234 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_xzz = cbuffer.data(hf_geom_1010_off + 1235 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_yyy = cbuffer.data(hf_geom_1010_off + 1236 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_yyz = cbuffer.data(hf_geom_1010_off + 1237 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_yzz = cbuffer.data(hf_geom_1010_off + 1238 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_zzz = cbuffer.data(hf_geom_1010_off + 1239 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1240 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1241 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1242 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1243 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1244 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1245 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1246 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1247 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1248 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1249 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1250 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1251 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1252 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1253 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1254 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1255 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1256 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1257 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1258 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1259 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_xxx = cbuffer.data(hf_geom_1010_off + 1260 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_xxy = cbuffer.data(hf_geom_1010_off + 1261 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_xxz = cbuffer.data(hf_geom_1010_off + 1262 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_xyy = cbuffer.data(hf_geom_1010_off + 1263 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_xyz = cbuffer.data(hf_geom_1010_off + 1264 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_xzz = cbuffer.data(hf_geom_1010_off + 1265 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_yyy = cbuffer.data(hf_geom_1010_off + 1266 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_yyz = cbuffer.data(hf_geom_1010_off + 1267 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_yzz = cbuffer.data(hf_geom_1010_off + 1268 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_zzz = cbuffer.data(hf_geom_1010_off + 1269 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_xxx = cbuffer.data(hf_geom_1010_off + 1270 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_xxy = cbuffer.data(hf_geom_1010_off + 1271 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_xxz = cbuffer.data(hf_geom_1010_off + 1272 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_xyy = cbuffer.data(hf_geom_1010_off + 1273 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_xyz = cbuffer.data(hf_geom_1010_off + 1274 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_xzz = cbuffer.data(hf_geom_1010_off + 1275 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_yyy = cbuffer.data(hf_geom_1010_off + 1276 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_yyz = cbuffer.data(hf_geom_1010_off + 1277 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_yzz = cbuffer.data(hf_geom_1010_off + 1278 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_zzz = cbuffer.data(hf_geom_1010_off + 1279 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_xxx = cbuffer.data(hf_geom_1010_off + 1280 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_xxy = cbuffer.data(hf_geom_1010_off + 1281 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_xxz = cbuffer.data(hf_geom_1010_off + 1282 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_xyy = cbuffer.data(hf_geom_1010_off + 1283 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_xyz = cbuffer.data(hf_geom_1010_off + 1284 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_xzz = cbuffer.data(hf_geom_1010_off + 1285 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_yyy = cbuffer.data(hf_geom_1010_off + 1286 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_yyz = cbuffer.data(hf_geom_1010_off + 1287 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_yzz = cbuffer.data(hf_geom_1010_off + 1288 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_zzz = cbuffer.data(hf_geom_1010_off + 1289 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_xxx = cbuffer.data(hf_geom_1010_off + 1290 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_xxy = cbuffer.data(hf_geom_1010_off + 1291 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_xxz = cbuffer.data(hf_geom_1010_off + 1292 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_xyy = cbuffer.data(hf_geom_1010_off + 1293 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_xyz = cbuffer.data(hf_geom_1010_off + 1294 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_xzz = cbuffer.data(hf_geom_1010_off + 1295 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_yyy = cbuffer.data(hf_geom_1010_off + 1296 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_yyz = cbuffer.data(hf_geom_1010_off + 1297 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_yzz = cbuffer.data(hf_geom_1010_off + 1298 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_zzz = cbuffer.data(hf_geom_1010_off + 1299 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_xxx = cbuffer.data(hf_geom_1010_off + 1300 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_xxy = cbuffer.data(hf_geom_1010_off + 1301 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_xxz = cbuffer.data(hf_geom_1010_off + 1302 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_xyy = cbuffer.data(hf_geom_1010_off + 1303 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_xyz = cbuffer.data(hf_geom_1010_off + 1304 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_xzz = cbuffer.data(hf_geom_1010_off + 1305 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_yyy = cbuffer.data(hf_geom_1010_off + 1306 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_yyz = cbuffer.data(hf_geom_1010_off + 1307 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_yzz = cbuffer.data(hf_geom_1010_off + 1308 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_zzz = cbuffer.data(hf_geom_1010_off + 1309 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_xxx = cbuffer.data(hf_geom_1010_off + 1310 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_xxy = cbuffer.data(hf_geom_1010_off + 1311 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_xxz = cbuffer.data(hf_geom_1010_off + 1312 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_xyy = cbuffer.data(hf_geom_1010_off + 1313 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_xyz = cbuffer.data(hf_geom_1010_off + 1314 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_xzz = cbuffer.data(hf_geom_1010_off + 1315 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_yyy = cbuffer.data(hf_geom_1010_off + 1316 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_yyz = cbuffer.data(hf_geom_1010_off + 1317 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_yzz = cbuffer.data(hf_geom_1010_off + 1318 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_zzz = cbuffer.data(hf_geom_1010_off + 1319 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_xxx = cbuffer.data(hf_geom_1010_off + 1320 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_xxy = cbuffer.data(hf_geom_1010_off + 1321 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_xxz = cbuffer.data(hf_geom_1010_off + 1322 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_xyy = cbuffer.data(hf_geom_1010_off + 1323 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_xyz = cbuffer.data(hf_geom_1010_off + 1324 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_xzz = cbuffer.data(hf_geom_1010_off + 1325 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_yyy = cbuffer.data(hf_geom_1010_off + 1326 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_yyz = cbuffer.data(hf_geom_1010_off + 1327 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_yzz = cbuffer.data(hf_geom_1010_off + 1328 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_zzz = cbuffer.data(hf_geom_1010_off + 1329 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_xxx = cbuffer.data(hf_geom_1010_off + 1330 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_xxy = cbuffer.data(hf_geom_1010_off + 1331 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_xxz = cbuffer.data(hf_geom_1010_off + 1332 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_xyy = cbuffer.data(hf_geom_1010_off + 1333 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_xyz = cbuffer.data(hf_geom_1010_off + 1334 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_xzz = cbuffer.data(hf_geom_1010_off + 1335 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_yyy = cbuffer.data(hf_geom_1010_off + 1336 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_yyz = cbuffer.data(hf_geom_1010_off + 1337 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_yzz = cbuffer.data(hf_geom_1010_off + 1338 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_zzz = cbuffer.data(hf_geom_1010_off + 1339 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_xxx = cbuffer.data(hf_geom_1010_off + 1340 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_xxy = cbuffer.data(hf_geom_1010_off + 1341 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_xxz = cbuffer.data(hf_geom_1010_off + 1342 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_xyy = cbuffer.data(hf_geom_1010_off + 1343 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_xyz = cbuffer.data(hf_geom_1010_off + 1344 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_xzz = cbuffer.data(hf_geom_1010_off + 1345 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_yyy = cbuffer.data(hf_geom_1010_off + 1346 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_yyz = cbuffer.data(hf_geom_1010_off + 1347 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_yzz = cbuffer.data(hf_geom_1010_off + 1348 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_zzz = cbuffer.data(hf_geom_1010_off + 1349 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_xxx = cbuffer.data(hf_geom_1010_off + 1350 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_xxy = cbuffer.data(hf_geom_1010_off + 1351 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_xxz = cbuffer.data(hf_geom_1010_off + 1352 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_xyy = cbuffer.data(hf_geom_1010_off + 1353 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_xyz = cbuffer.data(hf_geom_1010_off + 1354 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_xzz = cbuffer.data(hf_geom_1010_off + 1355 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_yyy = cbuffer.data(hf_geom_1010_off + 1356 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_yyz = cbuffer.data(hf_geom_1010_off + 1357 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_yzz = cbuffer.data(hf_geom_1010_off + 1358 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_zzz = cbuffer.data(hf_geom_1010_off + 1359 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_xxx = cbuffer.data(hf_geom_1010_off + 1360 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_xxy = cbuffer.data(hf_geom_1010_off + 1361 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_xxz = cbuffer.data(hf_geom_1010_off + 1362 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_xyy = cbuffer.data(hf_geom_1010_off + 1363 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_xyz = cbuffer.data(hf_geom_1010_off + 1364 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_xzz = cbuffer.data(hf_geom_1010_off + 1365 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_yyy = cbuffer.data(hf_geom_1010_off + 1366 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_yyz = cbuffer.data(hf_geom_1010_off + 1367 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_yzz = cbuffer.data(hf_geom_1010_off + 1368 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_zzz = cbuffer.data(hf_geom_1010_off + 1369 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_xxx = cbuffer.data(hf_geom_1010_off + 1370 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_xxy = cbuffer.data(hf_geom_1010_off + 1371 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_xxz = cbuffer.data(hf_geom_1010_off + 1372 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_xyy = cbuffer.data(hf_geom_1010_off + 1373 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_xyz = cbuffer.data(hf_geom_1010_off + 1374 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_xzz = cbuffer.data(hf_geom_1010_off + 1375 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_yyy = cbuffer.data(hf_geom_1010_off + 1376 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_yyz = cbuffer.data(hf_geom_1010_off + 1377 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_yzz = cbuffer.data(hf_geom_1010_off + 1378 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_zzz = cbuffer.data(hf_geom_1010_off + 1379 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_xxx = cbuffer.data(hf_geom_1010_off + 1380 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_xxy = cbuffer.data(hf_geom_1010_off + 1381 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_xxz = cbuffer.data(hf_geom_1010_off + 1382 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_xyy = cbuffer.data(hf_geom_1010_off + 1383 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_xyz = cbuffer.data(hf_geom_1010_off + 1384 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_xzz = cbuffer.data(hf_geom_1010_off + 1385 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_yyy = cbuffer.data(hf_geom_1010_off + 1386 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_yyz = cbuffer.data(hf_geom_1010_off + 1387 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_yzz = cbuffer.data(hf_geom_1010_off + 1388 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_zzz = cbuffer.data(hf_geom_1010_off + 1389 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_xxx = cbuffer.data(hf_geom_1010_off + 1390 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_xxy = cbuffer.data(hf_geom_1010_off + 1391 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_xxz = cbuffer.data(hf_geom_1010_off + 1392 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_xyy = cbuffer.data(hf_geom_1010_off + 1393 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_xyz = cbuffer.data(hf_geom_1010_off + 1394 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_xzz = cbuffer.data(hf_geom_1010_off + 1395 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_yyy = cbuffer.data(hf_geom_1010_off + 1396 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_yyz = cbuffer.data(hf_geom_1010_off + 1397 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_yzz = cbuffer.data(hf_geom_1010_off + 1398 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_zzz = cbuffer.data(hf_geom_1010_off + 1399 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1400 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1401 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1402 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1403 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1404 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1405 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1406 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1407 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1408 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1409 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_xxx = cbuffer.data(hf_geom_1010_off + 1410 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_xxy = cbuffer.data(hf_geom_1010_off + 1411 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_xxz = cbuffer.data(hf_geom_1010_off + 1412 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_xyy = cbuffer.data(hf_geom_1010_off + 1413 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_xyz = cbuffer.data(hf_geom_1010_off + 1414 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_xzz = cbuffer.data(hf_geom_1010_off + 1415 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_yyy = cbuffer.data(hf_geom_1010_off + 1416 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_yyz = cbuffer.data(hf_geom_1010_off + 1417 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_yzz = cbuffer.data(hf_geom_1010_off + 1418 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_zzz = cbuffer.data(hf_geom_1010_off + 1419 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_xxx = cbuffer.data(hf_geom_1010_off + 1420 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_xxy = cbuffer.data(hf_geom_1010_off + 1421 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_xxz = cbuffer.data(hf_geom_1010_off + 1422 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_xyy = cbuffer.data(hf_geom_1010_off + 1423 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_xyz = cbuffer.data(hf_geom_1010_off + 1424 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_xzz = cbuffer.data(hf_geom_1010_off + 1425 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_yyy = cbuffer.data(hf_geom_1010_off + 1426 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_yyz = cbuffer.data(hf_geom_1010_off + 1427 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_yzz = cbuffer.data(hf_geom_1010_off + 1428 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_zzz = cbuffer.data(hf_geom_1010_off + 1429 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_xxx = cbuffer.data(hf_geom_1010_off + 1430 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_xxy = cbuffer.data(hf_geom_1010_off + 1431 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_xxz = cbuffer.data(hf_geom_1010_off + 1432 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_xyy = cbuffer.data(hf_geom_1010_off + 1433 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_xyz = cbuffer.data(hf_geom_1010_off + 1434 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_xzz = cbuffer.data(hf_geom_1010_off + 1435 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_yyy = cbuffer.data(hf_geom_1010_off + 1436 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_yyz = cbuffer.data(hf_geom_1010_off + 1437 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_yzz = cbuffer.data(hf_geom_1010_off + 1438 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_zzz = cbuffer.data(hf_geom_1010_off + 1439 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_xxx = cbuffer.data(hf_geom_1010_off + 1440 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_xxy = cbuffer.data(hf_geom_1010_off + 1441 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_xxz = cbuffer.data(hf_geom_1010_off + 1442 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_xyy = cbuffer.data(hf_geom_1010_off + 1443 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_xyz = cbuffer.data(hf_geom_1010_off + 1444 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_xzz = cbuffer.data(hf_geom_1010_off + 1445 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_yyy = cbuffer.data(hf_geom_1010_off + 1446 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_yyz = cbuffer.data(hf_geom_1010_off + 1447 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_yzz = cbuffer.data(hf_geom_1010_off + 1448 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_zzz = cbuffer.data(hf_geom_1010_off + 1449 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1450 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1451 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1452 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1453 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1454 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1455 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1456 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1457 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1458 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1459 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1460 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1461 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1462 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1463 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1464 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1465 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1466 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1467 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1468 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1469 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_xxx = cbuffer.data(hf_geom_1010_off + 1470 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_xxy = cbuffer.data(hf_geom_1010_off + 1471 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_xxz = cbuffer.data(hf_geom_1010_off + 1472 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_xyy = cbuffer.data(hf_geom_1010_off + 1473 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_xyz = cbuffer.data(hf_geom_1010_off + 1474 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_xzz = cbuffer.data(hf_geom_1010_off + 1475 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_yyy = cbuffer.data(hf_geom_1010_off + 1476 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_yyz = cbuffer.data(hf_geom_1010_off + 1477 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_yzz = cbuffer.data(hf_geom_1010_off + 1478 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_zzz = cbuffer.data(hf_geom_1010_off + 1479 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_xxx = cbuffer.data(hf_geom_1010_off + 1480 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_xxy = cbuffer.data(hf_geom_1010_off + 1481 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_xxz = cbuffer.data(hf_geom_1010_off + 1482 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_xyy = cbuffer.data(hf_geom_1010_off + 1483 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_xyz = cbuffer.data(hf_geom_1010_off + 1484 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_xzz = cbuffer.data(hf_geom_1010_off + 1485 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_yyy = cbuffer.data(hf_geom_1010_off + 1486 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_yyz = cbuffer.data(hf_geom_1010_off + 1487 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_yzz = cbuffer.data(hf_geom_1010_off + 1488 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_zzz = cbuffer.data(hf_geom_1010_off + 1489 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_xxx = cbuffer.data(hf_geom_1010_off + 1490 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_xxy = cbuffer.data(hf_geom_1010_off + 1491 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_xxz = cbuffer.data(hf_geom_1010_off + 1492 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_xyy = cbuffer.data(hf_geom_1010_off + 1493 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_xyz = cbuffer.data(hf_geom_1010_off + 1494 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_xzz = cbuffer.data(hf_geom_1010_off + 1495 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_yyy = cbuffer.data(hf_geom_1010_off + 1496 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_yyz = cbuffer.data(hf_geom_1010_off + 1497 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_yzz = cbuffer.data(hf_geom_1010_off + 1498 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_zzz = cbuffer.data(hf_geom_1010_off + 1499 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_xxx = cbuffer.data(hf_geom_1010_off + 1500 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_xxy = cbuffer.data(hf_geom_1010_off + 1501 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_xxz = cbuffer.data(hf_geom_1010_off + 1502 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_xyy = cbuffer.data(hf_geom_1010_off + 1503 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_xyz = cbuffer.data(hf_geom_1010_off + 1504 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_xzz = cbuffer.data(hf_geom_1010_off + 1505 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_yyy = cbuffer.data(hf_geom_1010_off + 1506 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_yyz = cbuffer.data(hf_geom_1010_off + 1507 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_yzz = cbuffer.data(hf_geom_1010_off + 1508 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_zzz = cbuffer.data(hf_geom_1010_off + 1509 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_xxx = cbuffer.data(hf_geom_1010_off + 1510 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_xxy = cbuffer.data(hf_geom_1010_off + 1511 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_xxz = cbuffer.data(hf_geom_1010_off + 1512 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_xyy = cbuffer.data(hf_geom_1010_off + 1513 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_xyz = cbuffer.data(hf_geom_1010_off + 1514 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_xzz = cbuffer.data(hf_geom_1010_off + 1515 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_yyy = cbuffer.data(hf_geom_1010_off + 1516 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_yyz = cbuffer.data(hf_geom_1010_off + 1517 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_yzz = cbuffer.data(hf_geom_1010_off + 1518 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_zzz = cbuffer.data(hf_geom_1010_off + 1519 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_xxx = cbuffer.data(hf_geom_1010_off + 1520 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_xxy = cbuffer.data(hf_geom_1010_off + 1521 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_xxz = cbuffer.data(hf_geom_1010_off + 1522 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_xyy = cbuffer.data(hf_geom_1010_off + 1523 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_xyz = cbuffer.data(hf_geom_1010_off + 1524 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_xzz = cbuffer.data(hf_geom_1010_off + 1525 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_yyy = cbuffer.data(hf_geom_1010_off + 1526 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_yyz = cbuffer.data(hf_geom_1010_off + 1527 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_yzz = cbuffer.data(hf_geom_1010_off + 1528 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_zzz = cbuffer.data(hf_geom_1010_off + 1529 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_xxx = cbuffer.data(hf_geom_1010_off + 1530 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_xxy = cbuffer.data(hf_geom_1010_off + 1531 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_xxz = cbuffer.data(hf_geom_1010_off + 1532 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_xyy = cbuffer.data(hf_geom_1010_off + 1533 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_xyz = cbuffer.data(hf_geom_1010_off + 1534 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_xzz = cbuffer.data(hf_geom_1010_off + 1535 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_yyy = cbuffer.data(hf_geom_1010_off + 1536 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_yyz = cbuffer.data(hf_geom_1010_off + 1537 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_yzz = cbuffer.data(hf_geom_1010_off + 1538 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_zzz = cbuffer.data(hf_geom_1010_off + 1539 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_xxx = cbuffer.data(hf_geom_1010_off + 1540 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_xxy = cbuffer.data(hf_geom_1010_off + 1541 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_xxz = cbuffer.data(hf_geom_1010_off + 1542 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_xyy = cbuffer.data(hf_geom_1010_off + 1543 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_xyz = cbuffer.data(hf_geom_1010_off + 1544 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_xzz = cbuffer.data(hf_geom_1010_off + 1545 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_yyy = cbuffer.data(hf_geom_1010_off + 1546 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_yyz = cbuffer.data(hf_geom_1010_off + 1547 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_yzz = cbuffer.data(hf_geom_1010_off + 1548 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_zzz = cbuffer.data(hf_geom_1010_off + 1549 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_xxx = cbuffer.data(hf_geom_1010_off + 1550 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_xxy = cbuffer.data(hf_geom_1010_off + 1551 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_xxz = cbuffer.data(hf_geom_1010_off + 1552 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_xyy = cbuffer.data(hf_geom_1010_off + 1553 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_xyz = cbuffer.data(hf_geom_1010_off + 1554 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_xzz = cbuffer.data(hf_geom_1010_off + 1555 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_yyy = cbuffer.data(hf_geom_1010_off + 1556 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_yyz = cbuffer.data(hf_geom_1010_off + 1557 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_yzz = cbuffer.data(hf_geom_1010_off + 1558 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_zzz = cbuffer.data(hf_geom_1010_off + 1559 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_xxx = cbuffer.data(hf_geom_1010_off + 1560 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_xxy = cbuffer.data(hf_geom_1010_off + 1561 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_xxz = cbuffer.data(hf_geom_1010_off + 1562 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_xyy = cbuffer.data(hf_geom_1010_off + 1563 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_xyz = cbuffer.data(hf_geom_1010_off + 1564 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_xzz = cbuffer.data(hf_geom_1010_off + 1565 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_yyy = cbuffer.data(hf_geom_1010_off + 1566 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_yyz = cbuffer.data(hf_geom_1010_off + 1567 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_yzz = cbuffer.data(hf_geom_1010_off + 1568 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_zzz = cbuffer.data(hf_geom_1010_off + 1569 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_xxx = cbuffer.data(hf_geom_1010_off + 1570 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_xxy = cbuffer.data(hf_geom_1010_off + 1571 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_xxz = cbuffer.data(hf_geom_1010_off + 1572 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_xyy = cbuffer.data(hf_geom_1010_off + 1573 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_xyz = cbuffer.data(hf_geom_1010_off + 1574 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_xzz = cbuffer.data(hf_geom_1010_off + 1575 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_yyy = cbuffer.data(hf_geom_1010_off + 1576 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_yyz = cbuffer.data(hf_geom_1010_off + 1577 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_yzz = cbuffer.data(hf_geom_1010_off + 1578 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_zzz = cbuffer.data(hf_geom_1010_off + 1579 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_xxx = cbuffer.data(hf_geom_1010_off + 1580 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_xxy = cbuffer.data(hf_geom_1010_off + 1581 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_xxz = cbuffer.data(hf_geom_1010_off + 1582 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_xyy = cbuffer.data(hf_geom_1010_off + 1583 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_xyz = cbuffer.data(hf_geom_1010_off + 1584 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_xzz = cbuffer.data(hf_geom_1010_off + 1585 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_yyy = cbuffer.data(hf_geom_1010_off + 1586 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_yyz = cbuffer.data(hf_geom_1010_off + 1587 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_yzz = cbuffer.data(hf_geom_1010_off + 1588 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_zzz = cbuffer.data(hf_geom_1010_off + 1589 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_xxx = cbuffer.data(hf_geom_1010_off + 1590 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_xxy = cbuffer.data(hf_geom_1010_off + 1591 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_xxz = cbuffer.data(hf_geom_1010_off + 1592 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_xyy = cbuffer.data(hf_geom_1010_off + 1593 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_xyz = cbuffer.data(hf_geom_1010_off + 1594 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_xzz = cbuffer.data(hf_geom_1010_off + 1595 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_yyy = cbuffer.data(hf_geom_1010_off + 1596 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_yyz = cbuffer.data(hf_geom_1010_off + 1597 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_yzz = cbuffer.data(hf_geom_1010_off + 1598 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_zzz = cbuffer.data(hf_geom_1010_off + 1599 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_xxx = cbuffer.data(hf_geom_1010_off + 1600 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_xxy = cbuffer.data(hf_geom_1010_off + 1601 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_xxz = cbuffer.data(hf_geom_1010_off + 1602 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_xyy = cbuffer.data(hf_geom_1010_off + 1603 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_xyz = cbuffer.data(hf_geom_1010_off + 1604 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_xzz = cbuffer.data(hf_geom_1010_off + 1605 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_yyy = cbuffer.data(hf_geom_1010_off + 1606 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_yyz = cbuffer.data(hf_geom_1010_off + 1607 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_yzz = cbuffer.data(hf_geom_1010_off + 1608 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_zzz = cbuffer.data(hf_geom_1010_off + 1609 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1610 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1611 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1612 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1613 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1614 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1615 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1616 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1617 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1618 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1619 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_xxx = cbuffer.data(hf_geom_1010_off + 1620 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_xxy = cbuffer.data(hf_geom_1010_off + 1621 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_xxz = cbuffer.data(hf_geom_1010_off + 1622 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_xyy = cbuffer.data(hf_geom_1010_off + 1623 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_xyz = cbuffer.data(hf_geom_1010_off + 1624 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_xzz = cbuffer.data(hf_geom_1010_off + 1625 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_yyy = cbuffer.data(hf_geom_1010_off + 1626 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_yyz = cbuffer.data(hf_geom_1010_off + 1627 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_yzz = cbuffer.data(hf_geom_1010_off + 1628 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_zzz = cbuffer.data(hf_geom_1010_off + 1629 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_xxx = cbuffer.data(hf_geom_1010_off + 1630 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_xxy = cbuffer.data(hf_geom_1010_off + 1631 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_xxz = cbuffer.data(hf_geom_1010_off + 1632 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_xyy = cbuffer.data(hf_geom_1010_off + 1633 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_xyz = cbuffer.data(hf_geom_1010_off + 1634 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_xzz = cbuffer.data(hf_geom_1010_off + 1635 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_yyy = cbuffer.data(hf_geom_1010_off + 1636 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_yyz = cbuffer.data(hf_geom_1010_off + 1637 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_yzz = cbuffer.data(hf_geom_1010_off + 1638 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_zzz = cbuffer.data(hf_geom_1010_off + 1639 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_xxx = cbuffer.data(hf_geom_1010_off + 1640 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_xxy = cbuffer.data(hf_geom_1010_off + 1641 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_xxz = cbuffer.data(hf_geom_1010_off + 1642 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_xyy = cbuffer.data(hf_geom_1010_off + 1643 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_xyz = cbuffer.data(hf_geom_1010_off + 1644 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_xzz = cbuffer.data(hf_geom_1010_off + 1645 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_yyy = cbuffer.data(hf_geom_1010_off + 1646 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_yyz = cbuffer.data(hf_geom_1010_off + 1647 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_yzz = cbuffer.data(hf_geom_1010_off + 1648 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_zzz = cbuffer.data(hf_geom_1010_off + 1649 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_xxx = cbuffer.data(hf_geom_1010_off + 1650 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_xxy = cbuffer.data(hf_geom_1010_off + 1651 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_xxz = cbuffer.data(hf_geom_1010_off + 1652 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_xyy = cbuffer.data(hf_geom_1010_off + 1653 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_xyz = cbuffer.data(hf_geom_1010_off + 1654 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_xzz = cbuffer.data(hf_geom_1010_off + 1655 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_yyy = cbuffer.data(hf_geom_1010_off + 1656 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_yyz = cbuffer.data(hf_geom_1010_off + 1657 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_yzz = cbuffer.data(hf_geom_1010_off + 1658 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_zzz = cbuffer.data(hf_geom_1010_off + 1659 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1660 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1661 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1662 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1663 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1664 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1665 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1666 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1667 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1668 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1669 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1670 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1671 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1672 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1673 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1674 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1675 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1676 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1677 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1678 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1679 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_xxx = cbuffer.data(hf_geom_1010_off + 1680 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_xxy = cbuffer.data(hf_geom_1010_off + 1681 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_xxz = cbuffer.data(hf_geom_1010_off + 1682 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_xyy = cbuffer.data(hf_geom_1010_off + 1683 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_xyz = cbuffer.data(hf_geom_1010_off + 1684 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_xzz = cbuffer.data(hf_geom_1010_off + 1685 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_yyy = cbuffer.data(hf_geom_1010_off + 1686 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_yyz = cbuffer.data(hf_geom_1010_off + 1687 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_yzz = cbuffer.data(hf_geom_1010_off + 1688 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_zzz = cbuffer.data(hf_geom_1010_off + 1689 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_xxx = cbuffer.data(hf_geom_1010_off + 1690 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_xxy = cbuffer.data(hf_geom_1010_off + 1691 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_xxz = cbuffer.data(hf_geom_1010_off + 1692 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_xyy = cbuffer.data(hf_geom_1010_off + 1693 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_xyz = cbuffer.data(hf_geom_1010_off + 1694 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_xzz = cbuffer.data(hf_geom_1010_off + 1695 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_yyy = cbuffer.data(hf_geom_1010_off + 1696 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_yyz = cbuffer.data(hf_geom_1010_off + 1697 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_yzz = cbuffer.data(hf_geom_1010_off + 1698 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_zzz = cbuffer.data(hf_geom_1010_off + 1699 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_xxx = cbuffer.data(hf_geom_1010_off + 1700 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_xxy = cbuffer.data(hf_geom_1010_off + 1701 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_xxz = cbuffer.data(hf_geom_1010_off + 1702 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_xyy = cbuffer.data(hf_geom_1010_off + 1703 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_xyz = cbuffer.data(hf_geom_1010_off + 1704 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_xzz = cbuffer.data(hf_geom_1010_off + 1705 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_yyy = cbuffer.data(hf_geom_1010_off + 1706 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_yyz = cbuffer.data(hf_geom_1010_off + 1707 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_yzz = cbuffer.data(hf_geom_1010_off + 1708 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_zzz = cbuffer.data(hf_geom_1010_off + 1709 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_xxx = cbuffer.data(hf_geom_1010_off + 1710 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_xxy = cbuffer.data(hf_geom_1010_off + 1711 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_xxz = cbuffer.data(hf_geom_1010_off + 1712 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_xyy = cbuffer.data(hf_geom_1010_off + 1713 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_xyz = cbuffer.data(hf_geom_1010_off + 1714 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_xzz = cbuffer.data(hf_geom_1010_off + 1715 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_yyy = cbuffer.data(hf_geom_1010_off + 1716 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_yyz = cbuffer.data(hf_geom_1010_off + 1717 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_yzz = cbuffer.data(hf_geom_1010_off + 1718 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_zzz = cbuffer.data(hf_geom_1010_off + 1719 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_xxx = cbuffer.data(hf_geom_1010_off + 1720 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_xxy = cbuffer.data(hf_geom_1010_off + 1721 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_xxz = cbuffer.data(hf_geom_1010_off + 1722 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_xyy = cbuffer.data(hf_geom_1010_off + 1723 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_xyz = cbuffer.data(hf_geom_1010_off + 1724 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_xzz = cbuffer.data(hf_geom_1010_off + 1725 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_yyy = cbuffer.data(hf_geom_1010_off + 1726 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_yyz = cbuffer.data(hf_geom_1010_off + 1727 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_yzz = cbuffer.data(hf_geom_1010_off + 1728 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_zzz = cbuffer.data(hf_geom_1010_off + 1729 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_xxx = cbuffer.data(hf_geom_1010_off + 1730 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_xxy = cbuffer.data(hf_geom_1010_off + 1731 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_xxz = cbuffer.data(hf_geom_1010_off + 1732 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_xyy = cbuffer.data(hf_geom_1010_off + 1733 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_xyz = cbuffer.data(hf_geom_1010_off + 1734 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_xzz = cbuffer.data(hf_geom_1010_off + 1735 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_yyy = cbuffer.data(hf_geom_1010_off + 1736 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_yyz = cbuffer.data(hf_geom_1010_off + 1737 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_yzz = cbuffer.data(hf_geom_1010_off + 1738 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_zzz = cbuffer.data(hf_geom_1010_off + 1739 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_xxx = cbuffer.data(hf_geom_1010_off + 1740 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_xxy = cbuffer.data(hf_geom_1010_off + 1741 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_xxz = cbuffer.data(hf_geom_1010_off + 1742 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_xyy = cbuffer.data(hf_geom_1010_off + 1743 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_xyz = cbuffer.data(hf_geom_1010_off + 1744 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_xzz = cbuffer.data(hf_geom_1010_off + 1745 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_yyy = cbuffer.data(hf_geom_1010_off + 1746 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_yyz = cbuffer.data(hf_geom_1010_off + 1747 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_yzz = cbuffer.data(hf_geom_1010_off + 1748 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_zzz = cbuffer.data(hf_geom_1010_off + 1749 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_xxx = cbuffer.data(hf_geom_1010_off + 1750 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_xxy = cbuffer.data(hf_geom_1010_off + 1751 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_xxz = cbuffer.data(hf_geom_1010_off + 1752 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_xyy = cbuffer.data(hf_geom_1010_off + 1753 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_xyz = cbuffer.data(hf_geom_1010_off + 1754 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_xzz = cbuffer.data(hf_geom_1010_off + 1755 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_yyy = cbuffer.data(hf_geom_1010_off + 1756 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_yyz = cbuffer.data(hf_geom_1010_off + 1757 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_yzz = cbuffer.data(hf_geom_1010_off + 1758 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_zzz = cbuffer.data(hf_geom_1010_off + 1759 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_xxx = cbuffer.data(hf_geom_1010_off + 1760 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_xxy = cbuffer.data(hf_geom_1010_off + 1761 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_xxz = cbuffer.data(hf_geom_1010_off + 1762 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_xyy = cbuffer.data(hf_geom_1010_off + 1763 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_xyz = cbuffer.data(hf_geom_1010_off + 1764 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_xzz = cbuffer.data(hf_geom_1010_off + 1765 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_yyy = cbuffer.data(hf_geom_1010_off + 1766 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_yyz = cbuffer.data(hf_geom_1010_off + 1767 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_yzz = cbuffer.data(hf_geom_1010_off + 1768 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_zzz = cbuffer.data(hf_geom_1010_off + 1769 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_xxx = cbuffer.data(hf_geom_1010_off + 1770 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_xxy = cbuffer.data(hf_geom_1010_off + 1771 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_xxz = cbuffer.data(hf_geom_1010_off + 1772 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_xyy = cbuffer.data(hf_geom_1010_off + 1773 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_xyz = cbuffer.data(hf_geom_1010_off + 1774 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_xzz = cbuffer.data(hf_geom_1010_off + 1775 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_yyy = cbuffer.data(hf_geom_1010_off + 1776 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_yyz = cbuffer.data(hf_geom_1010_off + 1777 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_yzz = cbuffer.data(hf_geom_1010_off + 1778 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_zzz = cbuffer.data(hf_geom_1010_off + 1779 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_xxx = cbuffer.data(hf_geom_1010_off + 1780 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_xxy = cbuffer.data(hf_geom_1010_off + 1781 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_xxz = cbuffer.data(hf_geom_1010_off + 1782 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_xyy = cbuffer.data(hf_geom_1010_off + 1783 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_xyz = cbuffer.data(hf_geom_1010_off + 1784 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_xzz = cbuffer.data(hf_geom_1010_off + 1785 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_yyy = cbuffer.data(hf_geom_1010_off + 1786 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_yyz = cbuffer.data(hf_geom_1010_off + 1787 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_yzz = cbuffer.data(hf_geom_1010_off + 1788 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_zzz = cbuffer.data(hf_geom_1010_off + 1789 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_xxx = cbuffer.data(hf_geom_1010_off + 1790 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_xxy = cbuffer.data(hf_geom_1010_off + 1791 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_xxz = cbuffer.data(hf_geom_1010_off + 1792 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_xyy = cbuffer.data(hf_geom_1010_off + 1793 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_xyz = cbuffer.data(hf_geom_1010_off + 1794 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_xzz = cbuffer.data(hf_geom_1010_off + 1795 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_yyy = cbuffer.data(hf_geom_1010_off + 1796 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_yyz = cbuffer.data(hf_geom_1010_off + 1797 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_yzz = cbuffer.data(hf_geom_1010_off + 1798 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_zzz = cbuffer.data(hf_geom_1010_off + 1799 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_xxx = cbuffer.data(hf_geom_1010_off + 1800 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_xxy = cbuffer.data(hf_geom_1010_off + 1801 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_xxz = cbuffer.data(hf_geom_1010_off + 1802 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_xyy = cbuffer.data(hf_geom_1010_off + 1803 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_xyz = cbuffer.data(hf_geom_1010_off + 1804 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_xzz = cbuffer.data(hf_geom_1010_off + 1805 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_yyy = cbuffer.data(hf_geom_1010_off + 1806 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_yyz = cbuffer.data(hf_geom_1010_off + 1807 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_yzz = cbuffer.data(hf_geom_1010_off + 1808 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_zzz = cbuffer.data(hf_geom_1010_off + 1809 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_xxx = cbuffer.data(hf_geom_1010_off + 1810 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_xxy = cbuffer.data(hf_geom_1010_off + 1811 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_xxz = cbuffer.data(hf_geom_1010_off + 1812 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_xyy = cbuffer.data(hf_geom_1010_off + 1813 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_xyz = cbuffer.data(hf_geom_1010_off + 1814 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_xzz = cbuffer.data(hf_geom_1010_off + 1815 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_yyy = cbuffer.data(hf_geom_1010_off + 1816 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_yyz = cbuffer.data(hf_geom_1010_off + 1817 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_yzz = cbuffer.data(hf_geom_1010_off + 1818 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_zzz = cbuffer.data(hf_geom_1010_off + 1819 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1820 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1821 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1822 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1823 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1824 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1825 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1826 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1827 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1828 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1829 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_xxx = cbuffer.data(hf_geom_1010_off + 1830 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_xxy = cbuffer.data(hf_geom_1010_off + 1831 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_xxz = cbuffer.data(hf_geom_1010_off + 1832 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_xyy = cbuffer.data(hf_geom_1010_off + 1833 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_xyz = cbuffer.data(hf_geom_1010_off + 1834 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_xzz = cbuffer.data(hf_geom_1010_off + 1835 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_yyy = cbuffer.data(hf_geom_1010_off + 1836 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_yyz = cbuffer.data(hf_geom_1010_off + 1837 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_yzz = cbuffer.data(hf_geom_1010_off + 1838 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_zzz = cbuffer.data(hf_geom_1010_off + 1839 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_xxx = cbuffer.data(hf_geom_1010_off + 1840 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_xxy = cbuffer.data(hf_geom_1010_off + 1841 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_xxz = cbuffer.data(hf_geom_1010_off + 1842 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_xyy = cbuffer.data(hf_geom_1010_off + 1843 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_xyz = cbuffer.data(hf_geom_1010_off + 1844 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_xzz = cbuffer.data(hf_geom_1010_off + 1845 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_yyy = cbuffer.data(hf_geom_1010_off + 1846 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_yyz = cbuffer.data(hf_geom_1010_off + 1847 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_yzz = cbuffer.data(hf_geom_1010_off + 1848 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_zzz = cbuffer.data(hf_geom_1010_off + 1849 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_xxx = cbuffer.data(hf_geom_1010_off + 1850 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_xxy = cbuffer.data(hf_geom_1010_off + 1851 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_xxz = cbuffer.data(hf_geom_1010_off + 1852 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_xyy = cbuffer.data(hf_geom_1010_off + 1853 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_xyz = cbuffer.data(hf_geom_1010_off + 1854 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_xzz = cbuffer.data(hf_geom_1010_off + 1855 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_yyy = cbuffer.data(hf_geom_1010_off + 1856 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_yyz = cbuffer.data(hf_geom_1010_off + 1857 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_yzz = cbuffer.data(hf_geom_1010_off + 1858 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_zzz = cbuffer.data(hf_geom_1010_off + 1859 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_xxx = cbuffer.data(hf_geom_1010_off + 1860 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_xxy = cbuffer.data(hf_geom_1010_off + 1861 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_xxz = cbuffer.data(hf_geom_1010_off + 1862 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_xyy = cbuffer.data(hf_geom_1010_off + 1863 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_xyz = cbuffer.data(hf_geom_1010_off + 1864 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_xzz = cbuffer.data(hf_geom_1010_off + 1865 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_yyy = cbuffer.data(hf_geom_1010_off + 1866 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_yyz = cbuffer.data(hf_geom_1010_off + 1867 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_yzz = cbuffer.data(hf_geom_1010_off + 1868 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_zzz = cbuffer.data(hf_geom_1010_off + 1869 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1870 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1871 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1872 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1873 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1874 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1875 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1876 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1877 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1878 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1879 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_xxx = cbuffer.data(hf_geom_1010_off + 1880 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_xxy = cbuffer.data(hf_geom_1010_off + 1881 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_xxz = cbuffer.data(hf_geom_1010_off + 1882 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_xyy = cbuffer.data(hf_geom_1010_off + 1883 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_xyz = cbuffer.data(hf_geom_1010_off + 1884 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_xzz = cbuffer.data(hf_geom_1010_off + 1885 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_yyy = cbuffer.data(hf_geom_1010_off + 1886 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_yyz = cbuffer.data(hf_geom_1010_off + 1887 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_yzz = cbuffer.data(hf_geom_1010_off + 1888 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_zzz = cbuffer.data(hf_geom_1010_off + 1889 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_idxx

            const auto id_geom_1010_off = idx_geom_1010_idxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxxx_xx = cbuffer.data(id_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxx_xy = cbuffer.data(id_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxx_xz = cbuffer.data(id_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxx_yy = cbuffer.data(id_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxx_yz = cbuffer.data(id_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxx_zz = cbuffer.data(id_geom_1010_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_xxxxx_xx, g_0_0_x_0_xxxxx_xy, g_0_0_x_0_xxxxx_xz, g_0_0_x_0_xxxxx_yy, g_0_0_x_0_xxxxx_yz, g_0_0_x_0_xxxxx_zz, g_x_0_x_0_xxxxx_xx, g_x_0_x_0_xxxxx_xxx, g_x_0_x_0_xxxxx_xxy, g_x_0_x_0_xxxxx_xxz, g_x_0_x_0_xxxxx_xy, g_x_0_x_0_xxxxx_xyy, g_x_0_x_0_xxxxx_xyz, g_x_0_x_0_xxxxx_xz, g_x_0_x_0_xxxxx_xzz, g_x_0_x_0_xxxxx_yy, g_x_0_x_0_xxxxx_yz, g_x_0_x_0_xxxxx_zz, g_x_0_x_0_xxxxxx_xx, g_x_0_x_0_xxxxxx_xy, g_x_0_x_0_xxxxxx_xz, g_x_0_x_0_xxxxxx_yy, g_x_0_x_0_xxxxxx_yz, g_x_0_x_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxxx_xx[k] = -g_0_0_x_0_xxxxx_xx[k] - g_x_0_x_0_xxxxx_xx[k] * ab_x + g_x_0_x_0_xxxxx_xxx[k];

                g_x_0_x_0_xxxxxx_xy[k] = -g_0_0_x_0_xxxxx_xy[k] - g_x_0_x_0_xxxxx_xy[k] * ab_x + g_x_0_x_0_xxxxx_xxy[k];

                g_x_0_x_0_xxxxxx_xz[k] = -g_0_0_x_0_xxxxx_xz[k] - g_x_0_x_0_xxxxx_xz[k] * ab_x + g_x_0_x_0_xxxxx_xxz[k];

                g_x_0_x_0_xxxxxx_yy[k] = -g_0_0_x_0_xxxxx_yy[k] - g_x_0_x_0_xxxxx_yy[k] * ab_x + g_x_0_x_0_xxxxx_xyy[k];

                g_x_0_x_0_xxxxxx_yz[k] = -g_0_0_x_0_xxxxx_yz[k] - g_x_0_x_0_xxxxx_yz[k] * ab_x + g_x_0_x_0_xxxxx_xyz[k];

                g_x_0_x_0_xxxxxx_zz[k] = -g_0_0_x_0_xxxxx_zz[k] - g_x_0_x_0_xxxxx_zz[k] * ab_x + g_x_0_x_0_xxxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxxy_xx = cbuffer.data(id_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxy_xy = cbuffer.data(id_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxy_xz = cbuffer.data(id_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxy_yy = cbuffer.data(id_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxy_yz = cbuffer.data(id_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxy_zz = cbuffer.data(id_geom_1010_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxx_xx, g_x_0_x_0_xxxxx_xxy, g_x_0_x_0_xxxxx_xy, g_x_0_x_0_xxxxx_xyy, g_x_0_x_0_xxxxx_xyz, g_x_0_x_0_xxxxx_xz, g_x_0_x_0_xxxxx_yy, g_x_0_x_0_xxxxx_yyy, g_x_0_x_0_xxxxx_yyz, g_x_0_x_0_xxxxx_yz, g_x_0_x_0_xxxxx_yzz, g_x_0_x_0_xxxxx_zz, g_x_0_x_0_xxxxxy_xx, g_x_0_x_0_xxxxxy_xy, g_x_0_x_0_xxxxxy_xz, g_x_0_x_0_xxxxxy_yy, g_x_0_x_0_xxxxxy_yz, g_x_0_x_0_xxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxxy_xx[k] = -g_x_0_x_0_xxxxx_xx[k] * ab_y + g_x_0_x_0_xxxxx_xxy[k];

                g_x_0_x_0_xxxxxy_xy[k] = -g_x_0_x_0_xxxxx_xy[k] * ab_y + g_x_0_x_0_xxxxx_xyy[k];

                g_x_0_x_0_xxxxxy_xz[k] = -g_x_0_x_0_xxxxx_xz[k] * ab_y + g_x_0_x_0_xxxxx_xyz[k];

                g_x_0_x_0_xxxxxy_yy[k] = -g_x_0_x_0_xxxxx_yy[k] * ab_y + g_x_0_x_0_xxxxx_yyy[k];

                g_x_0_x_0_xxxxxy_yz[k] = -g_x_0_x_0_xxxxx_yz[k] * ab_y + g_x_0_x_0_xxxxx_yyz[k];

                g_x_0_x_0_xxxxxy_zz[k] = -g_x_0_x_0_xxxxx_zz[k] * ab_y + g_x_0_x_0_xxxxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxxz_xx = cbuffer.data(id_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxz_xy = cbuffer.data(id_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxz_xz = cbuffer.data(id_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxz_yy = cbuffer.data(id_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxz_yz = cbuffer.data(id_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxz_zz = cbuffer.data(id_geom_1010_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxx_xx, g_x_0_x_0_xxxxx_xxz, g_x_0_x_0_xxxxx_xy, g_x_0_x_0_xxxxx_xyz, g_x_0_x_0_xxxxx_xz, g_x_0_x_0_xxxxx_xzz, g_x_0_x_0_xxxxx_yy, g_x_0_x_0_xxxxx_yyz, g_x_0_x_0_xxxxx_yz, g_x_0_x_0_xxxxx_yzz, g_x_0_x_0_xxxxx_zz, g_x_0_x_0_xxxxx_zzz, g_x_0_x_0_xxxxxz_xx, g_x_0_x_0_xxxxxz_xy, g_x_0_x_0_xxxxxz_xz, g_x_0_x_0_xxxxxz_yy, g_x_0_x_0_xxxxxz_yz, g_x_0_x_0_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxxz_xx[k] = -g_x_0_x_0_xxxxx_xx[k] * ab_z + g_x_0_x_0_xxxxx_xxz[k];

                g_x_0_x_0_xxxxxz_xy[k] = -g_x_0_x_0_xxxxx_xy[k] * ab_z + g_x_0_x_0_xxxxx_xyz[k];

                g_x_0_x_0_xxxxxz_xz[k] = -g_x_0_x_0_xxxxx_xz[k] * ab_z + g_x_0_x_0_xxxxx_xzz[k];

                g_x_0_x_0_xxxxxz_yy[k] = -g_x_0_x_0_xxxxx_yy[k] * ab_z + g_x_0_x_0_xxxxx_yyz[k];

                g_x_0_x_0_xxxxxz_yz[k] = -g_x_0_x_0_xxxxx_yz[k] * ab_z + g_x_0_x_0_xxxxx_yzz[k];

                g_x_0_x_0_xxxxxz_zz[k] = -g_x_0_x_0_xxxxx_zz[k] * ab_z + g_x_0_x_0_xxxxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxyy_xx = cbuffer.data(id_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyy_xy = cbuffer.data(id_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyy_xz = cbuffer.data(id_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyy_yy = cbuffer.data(id_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyy_yz = cbuffer.data(id_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyy_zz = cbuffer.data(id_geom_1010_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxy_xx, g_x_0_x_0_xxxxy_xxy, g_x_0_x_0_xxxxy_xy, g_x_0_x_0_xxxxy_xyy, g_x_0_x_0_xxxxy_xyz, g_x_0_x_0_xxxxy_xz, g_x_0_x_0_xxxxy_yy, g_x_0_x_0_xxxxy_yyy, g_x_0_x_0_xxxxy_yyz, g_x_0_x_0_xxxxy_yz, g_x_0_x_0_xxxxy_yzz, g_x_0_x_0_xxxxy_zz, g_x_0_x_0_xxxxyy_xx, g_x_0_x_0_xxxxyy_xy, g_x_0_x_0_xxxxyy_xz, g_x_0_x_0_xxxxyy_yy, g_x_0_x_0_xxxxyy_yz, g_x_0_x_0_xxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxyy_xx[k] = -g_x_0_x_0_xxxxy_xx[k] * ab_y + g_x_0_x_0_xxxxy_xxy[k];

                g_x_0_x_0_xxxxyy_xy[k] = -g_x_0_x_0_xxxxy_xy[k] * ab_y + g_x_0_x_0_xxxxy_xyy[k];

                g_x_0_x_0_xxxxyy_xz[k] = -g_x_0_x_0_xxxxy_xz[k] * ab_y + g_x_0_x_0_xxxxy_xyz[k];

                g_x_0_x_0_xxxxyy_yy[k] = -g_x_0_x_0_xxxxy_yy[k] * ab_y + g_x_0_x_0_xxxxy_yyy[k];

                g_x_0_x_0_xxxxyy_yz[k] = -g_x_0_x_0_xxxxy_yz[k] * ab_y + g_x_0_x_0_xxxxy_yyz[k];

                g_x_0_x_0_xxxxyy_zz[k] = -g_x_0_x_0_xxxxy_zz[k] * ab_y + g_x_0_x_0_xxxxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxyz_xx = cbuffer.data(id_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyz_xy = cbuffer.data(id_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyz_xz = cbuffer.data(id_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyz_yy = cbuffer.data(id_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyz_yz = cbuffer.data(id_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyz_zz = cbuffer.data(id_geom_1010_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxyz_xx, g_x_0_x_0_xxxxyz_xy, g_x_0_x_0_xxxxyz_xz, g_x_0_x_0_xxxxyz_yy, g_x_0_x_0_xxxxyz_yz, g_x_0_x_0_xxxxyz_zz, g_x_0_x_0_xxxxz_xx, g_x_0_x_0_xxxxz_xxy, g_x_0_x_0_xxxxz_xy, g_x_0_x_0_xxxxz_xyy, g_x_0_x_0_xxxxz_xyz, g_x_0_x_0_xxxxz_xz, g_x_0_x_0_xxxxz_yy, g_x_0_x_0_xxxxz_yyy, g_x_0_x_0_xxxxz_yyz, g_x_0_x_0_xxxxz_yz, g_x_0_x_0_xxxxz_yzz, g_x_0_x_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxyz_xx[k] = -g_x_0_x_0_xxxxz_xx[k] * ab_y + g_x_0_x_0_xxxxz_xxy[k];

                g_x_0_x_0_xxxxyz_xy[k] = -g_x_0_x_0_xxxxz_xy[k] * ab_y + g_x_0_x_0_xxxxz_xyy[k];

                g_x_0_x_0_xxxxyz_xz[k] = -g_x_0_x_0_xxxxz_xz[k] * ab_y + g_x_0_x_0_xxxxz_xyz[k];

                g_x_0_x_0_xxxxyz_yy[k] = -g_x_0_x_0_xxxxz_yy[k] * ab_y + g_x_0_x_0_xxxxz_yyy[k];

                g_x_0_x_0_xxxxyz_yz[k] = -g_x_0_x_0_xxxxz_yz[k] * ab_y + g_x_0_x_0_xxxxz_yyz[k];

                g_x_0_x_0_xxxxyz_zz[k] = -g_x_0_x_0_xxxxz_zz[k] * ab_y + g_x_0_x_0_xxxxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxzz_xx = cbuffer.data(id_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxzz_xy = cbuffer.data(id_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxzz_xz = cbuffer.data(id_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxzz_yy = cbuffer.data(id_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxzz_yz = cbuffer.data(id_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxzz_zz = cbuffer.data(id_geom_1010_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxz_xx, g_x_0_x_0_xxxxz_xxz, g_x_0_x_0_xxxxz_xy, g_x_0_x_0_xxxxz_xyz, g_x_0_x_0_xxxxz_xz, g_x_0_x_0_xxxxz_xzz, g_x_0_x_0_xxxxz_yy, g_x_0_x_0_xxxxz_yyz, g_x_0_x_0_xxxxz_yz, g_x_0_x_0_xxxxz_yzz, g_x_0_x_0_xxxxz_zz, g_x_0_x_0_xxxxz_zzz, g_x_0_x_0_xxxxzz_xx, g_x_0_x_0_xxxxzz_xy, g_x_0_x_0_xxxxzz_xz, g_x_0_x_0_xxxxzz_yy, g_x_0_x_0_xxxxzz_yz, g_x_0_x_0_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxzz_xx[k] = -g_x_0_x_0_xxxxz_xx[k] * ab_z + g_x_0_x_0_xxxxz_xxz[k];

                g_x_0_x_0_xxxxzz_xy[k] = -g_x_0_x_0_xxxxz_xy[k] * ab_z + g_x_0_x_0_xxxxz_xyz[k];

                g_x_0_x_0_xxxxzz_xz[k] = -g_x_0_x_0_xxxxz_xz[k] * ab_z + g_x_0_x_0_xxxxz_xzz[k];

                g_x_0_x_0_xxxxzz_yy[k] = -g_x_0_x_0_xxxxz_yy[k] * ab_z + g_x_0_x_0_xxxxz_yyz[k];

                g_x_0_x_0_xxxxzz_yz[k] = -g_x_0_x_0_xxxxz_yz[k] * ab_z + g_x_0_x_0_xxxxz_yzz[k];

                g_x_0_x_0_xxxxzz_zz[k] = -g_x_0_x_0_xxxxz_zz[k] * ab_z + g_x_0_x_0_xxxxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxyyy_xx = cbuffer.data(id_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyy_xy = cbuffer.data(id_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyy_xz = cbuffer.data(id_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyy_yy = cbuffer.data(id_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyy_yz = cbuffer.data(id_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyy_zz = cbuffer.data(id_geom_1010_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxyy_xx, g_x_0_x_0_xxxyy_xxy, g_x_0_x_0_xxxyy_xy, g_x_0_x_0_xxxyy_xyy, g_x_0_x_0_xxxyy_xyz, g_x_0_x_0_xxxyy_xz, g_x_0_x_0_xxxyy_yy, g_x_0_x_0_xxxyy_yyy, g_x_0_x_0_xxxyy_yyz, g_x_0_x_0_xxxyy_yz, g_x_0_x_0_xxxyy_yzz, g_x_0_x_0_xxxyy_zz, g_x_0_x_0_xxxyyy_xx, g_x_0_x_0_xxxyyy_xy, g_x_0_x_0_xxxyyy_xz, g_x_0_x_0_xxxyyy_yy, g_x_0_x_0_xxxyyy_yz, g_x_0_x_0_xxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxyyy_xx[k] = -g_x_0_x_0_xxxyy_xx[k] * ab_y + g_x_0_x_0_xxxyy_xxy[k];

                g_x_0_x_0_xxxyyy_xy[k] = -g_x_0_x_0_xxxyy_xy[k] * ab_y + g_x_0_x_0_xxxyy_xyy[k];

                g_x_0_x_0_xxxyyy_xz[k] = -g_x_0_x_0_xxxyy_xz[k] * ab_y + g_x_0_x_0_xxxyy_xyz[k];

                g_x_0_x_0_xxxyyy_yy[k] = -g_x_0_x_0_xxxyy_yy[k] * ab_y + g_x_0_x_0_xxxyy_yyy[k];

                g_x_0_x_0_xxxyyy_yz[k] = -g_x_0_x_0_xxxyy_yz[k] * ab_y + g_x_0_x_0_xxxyy_yyz[k];

                g_x_0_x_0_xxxyyy_zz[k] = -g_x_0_x_0_xxxyy_zz[k] * ab_y + g_x_0_x_0_xxxyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxyyz_xx = cbuffer.data(id_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyz_xy = cbuffer.data(id_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyz_xz = cbuffer.data(id_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyz_yy = cbuffer.data(id_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyz_yz = cbuffer.data(id_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyz_zz = cbuffer.data(id_geom_1010_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxyyz_xx, g_x_0_x_0_xxxyyz_xy, g_x_0_x_0_xxxyyz_xz, g_x_0_x_0_xxxyyz_yy, g_x_0_x_0_xxxyyz_yz, g_x_0_x_0_xxxyyz_zz, g_x_0_x_0_xxxyz_xx, g_x_0_x_0_xxxyz_xxy, g_x_0_x_0_xxxyz_xy, g_x_0_x_0_xxxyz_xyy, g_x_0_x_0_xxxyz_xyz, g_x_0_x_0_xxxyz_xz, g_x_0_x_0_xxxyz_yy, g_x_0_x_0_xxxyz_yyy, g_x_0_x_0_xxxyz_yyz, g_x_0_x_0_xxxyz_yz, g_x_0_x_0_xxxyz_yzz, g_x_0_x_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxyyz_xx[k] = -g_x_0_x_0_xxxyz_xx[k] * ab_y + g_x_0_x_0_xxxyz_xxy[k];

                g_x_0_x_0_xxxyyz_xy[k] = -g_x_0_x_0_xxxyz_xy[k] * ab_y + g_x_0_x_0_xxxyz_xyy[k];

                g_x_0_x_0_xxxyyz_xz[k] = -g_x_0_x_0_xxxyz_xz[k] * ab_y + g_x_0_x_0_xxxyz_xyz[k];

                g_x_0_x_0_xxxyyz_yy[k] = -g_x_0_x_0_xxxyz_yy[k] * ab_y + g_x_0_x_0_xxxyz_yyy[k];

                g_x_0_x_0_xxxyyz_yz[k] = -g_x_0_x_0_xxxyz_yz[k] * ab_y + g_x_0_x_0_xxxyz_yyz[k];

                g_x_0_x_0_xxxyyz_zz[k] = -g_x_0_x_0_xxxyz_zz[k] * ab_y + g_x_0_x_0_xxxyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxyzz_xx = cbuffer.data(id_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyzz_xy = cbuffer.data(id_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyzz_xz = cbuffer.data(id_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyzz_yy = cbuffer.data(id_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyzz_yz = cbuffer.data(id_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyzz_zz = cbuffer.data(id_geom_1010_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxyzz_xx, g_x_0_x_0_xxxyzz_xy, g_x_0_x_0_xxxyzz_xz, g_x_0_x_0_xxxyzz_yy, g_x_0_x_0_xxxyzz_yz, g_x_0_x_0_xxxyzz_zz, g_x_0_x_0_xxxzz_xx, g_x_0_x_0_xxxzz_xxy, g_x_0_x_0_xxxzz_xy, g_x_0_x_0_xxxzz_xyy, g_x_0_x_0_xxxzz_xyz, g_x_0_x_0_xxxzz_xz, g_x_0_x_0_xxxzz_yy, g_x_0_x_0_xxxzz_yyy, g_x_0_x_0_xxxzz_yyz, g_x_0_x_0_xxxzz_yz, g_x_0_x_0_xxxzz_yzz, g_x_0_x_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxyzz_xx[k] = -g_x_0_x_0_xxxzz_xx[k] * ab_y + g_x_0_x_0_xxxzz_xxy[k];

                g_x_0_x_0_xxxyzz_xy[k] = -g_x_0_x_0_xxxzz_xy[k] * ab_y + g_x_0_x_0_xxxzz_xyy[k];

                g_x_0_x_0_xxxyzz_xz[k] = -g_x_0_x_0_xxxzz_xz[k] * ab_y + g_x_0_x_0_xxxzz_xyz[k];

                g_x_0_x_0_xxxyzz_yy[k] = -g_x_0_x_0_xxxzz_yy[k] * ab_y + g_x_0_x_0_xxxzz_yyy[k];

                g_x_0_x_0_xxxyzz_yz[k] = -g_x_0_x_0_xxxzz_yz[k] * ab_y + g_x_0_x_0_xxxzz_yyz[k];

                g_x_0_x_0_xxxyzz_zz[k] = -g_x_0_x_0_xxxzz_zz[k] * ab_y + g_x_0_x_0_xxxzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxzzz_xx = cbuffer.data(id_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzzz_xy = cbuffer.data(id_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzzz_xz = cbuffer.data(id_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzzz_yy = cbuffer.data(id_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzzz_yz = cbuffer.data(id_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzzz_zz = cbuffer.data(id_geom_1010_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxzz_xx, g_x_0_x_0_xxxzz_xxz, g_x_0_x_0_xxxzz_xy, g_x_0_x_0_xxxzz_xyz, g_x_0_x_0_xxxzz_xz, g_x_0_x_0_xxxzz_xzz, g_x_0_x_0_xxxzz_yy, g_x_0_x_0_xxxzz_yyz, g_x_0_x_0_xxxzz_yz, g_x_0_x_0_xxxzz_yzz, g_x_0_x_0_xxxzz_zz, g_x_0_x_0_xxxzz_zzz, g_x_0_x_0_xxxzzz_xx, g_x_0_x_0_xxxzzz_xy, g_x_0_x_0_xxxzzz_xz, g_x_0_x_0_xxxzzz_yy, g_x_0_x_0_xxxzzz_yz, g_x_0_x_0_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxzzz_xx[k] = -g_x_0_x_0_xxxzz_xx[k] * ab_z + g_x_0_x_0_xxxzz_xxz[k];

                g_x_0_x_0_xxxzzz_xy[k] = -g_x_0_x_0_xxxzz_xy[k] * ab_z + g_x_0_x_0_xxxzz_xyz[k];

                g_x_0_x_0_xxxzzz_xz[k] = -g_x_0_x_0_xxxzz_xz[k] * ab_z + g_x_0_x_0_xxxzz_xzz[k];

                g_x_0_x_0_xxxzzz_yy[k] = -g_x_0_x_0_xxxzz_yy[k] * ab_z + g_x_0_x_0_xxxzz_yyz[k];

                g_x_0_x_0_xxxzzz_yz[k] = -g_x_0_x_0_xxxzz_yz[k] * ab_z + g_x_0_x_0_xxxzz_yzz[k];

                g_x_0_x_0_xxxzzz_zz[k] = -g_x_0_x_0_xxxzz_zz[k] * ab_z + g_x_0_x_0_xxxzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyyyy_xx = cbuffer.data(id_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyy_xy = cbuffer.data(id_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyy_xz = cbuffer.data(id_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyy_yy = cbuffer.data(id_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyy_yz = cbuffer.data(id_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyy_zz = cbuffer.data(id_geom_1010_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyyy_xx, g_x_0_x_0_xxyyy_xxy, g_x_0_x_0_xxyyy_xy, g_x_0_x_0_xxyyy_xyy, g_x_0_x_0_xxyyy_xyz, g_x_0_x_0_xxyyy_xz, g_x_0_x_0_xxyyy_yy, g_x_0_x_0_xxyyy_yyy, g_x_0_x_0_xxyyy_yyz, g_x_0_x_0_xxyyy_yz, g_x_0_x_0_xxyyy_yzz, g_x_0_x_0_xxyyy_zz, g_x_0_x_0_xxyyyy_xx, g_x_0_x_0_xxyyyy_xy, g_x_0_x_0_xxyyyy_xz, g_x_0_x_0_xxyyyy_yy, g_x_0_x_0_xxyyyy_yz, g_x_0_x_0_xxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyyyy_xx[k] = -g_x_0_x_0_xxyyy_xx[k] * ab_y + g_x_0_x_0_xxyyy_xxy[k];

                g_x_0_x_0_xxyyyy_xy[k] = -g_x_0_x_0_xxyyy_xy[k] * ab_y + g_x_0_x_0_xxyyy_xyy[k];

                g_x_0_x_0_xxyyyy_xz[k] = -g_x_0_x_0_xxyyy_xz[k] * ab_y + g_x_0_x_0_xxyyy_xyz[k];

                g_x_0_x_0_xxyyyy_yy[k] = -g_x_0_x_0_xxyyy_yy[k] * ab_y + g_x_0_x_0_xxyyy_yyy[k];

                g_x_0_x_0_xxyyyy_yz[k] = -g_x_0_x_0_xxyyy_yz[k] * ab_y + g_x_0_x_0_xxyyy_yyz[k];

                g_x_0_x_0_xxyyyy_zz[k] = -g_x_0_x_0_xxyyy_zz[k] * ab_y + g_x_0_x_0_xxyyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyyyz_xx = cbuffer.data(id_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyz_xy = cbuffer.data(id_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyz_xz = cbuffer.data(id_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyz_yy = cbuffer.data(id_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyz_yz = cbuffer.data(id_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyz_zz = cbuffer.data(id_geom_1010_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyyyz_xx, g_x_0_x_0_xxyyyz_xy, g_x_0_x_0_xxyyyz_xz, g_x_0_x_0_xxyyyz_yy, g_x_0_x_0_xxyyyz_yz, g_x_0_x_0_xxyyyz_zz, g_x_0_x_0_xxyyz_xx, g_x_0_x_0_xxyyz_xxy, g_x_0_x_0_xxyyz_xy, g_x_0_x_0_xxyyz_xyy, g_x_0_x_0_xxyyz_xyz, g_x_0_x_0_xxyyz_xz, g_x_0_x_0_xxyyz_yy, g_x_0_x_0_xxyyz_yyy, g_x_0_x_0_xxyyz_yyz, g_x_0_x_0_xxyyz_yz, g_x_0_x_0_xxyyz_yzz, g_x_0_x_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyyyz_xx[k] = -g_x_0_x_0_xxyyz_xx[k] * ab_y + g_x_0_x_0_xxyyz_xxy[k];

                g_x_0_x_0_xxyyyz_xy[k] = -g_x_0_x_0_xxyyz_xy[k] * ab_y + g_x_0_x_0_xxyyz_xyy[k];

                g_x_0_x_0_xxyyyz_xz[k] = -g_x_0_x_0_xxyyz_xz[k] * ab_y + g_x_0_x_0_xxyyz_xyz[k];

                g_x_0_x_0_xxyyyz_yy[k] = -g_x_0_x_0_xxyyz_yy[k] * ab_y + g_x_0_x_0_xxyyz_yyy[k];

                g_x_0_x_0_xxyyyz_yz[k] = -g_x_0_x_0_xxyyz_yz[k] * ab_y + g_x_0_x_0_xxyyz_yyz[k];

                g_x_0_x_0_xxyyyz_zz[k] = -g_x_0_x_0_xxyyz_zz[k] * ab_y + g_x_0_x_0_xxyyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyyzz_xx = cbuffer.data(id_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyzz_xy = cbuffer.data(id_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyzz_xz = cbuffer.data(id_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyzz_yy = cbuffer.data(id_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyzz_yz = cbuffer.data(id_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyzz_zz = cbuffer.data(id_geom_1010_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyyzz_xx, g_x_0_x_0_xxyyzz_xy, g_x_0_x_0_xxyyzz_xz, g_x_0_x_0_xxyyzz_yy, g_x_0_x_0_xxyyzz_yz, g_x_0_x_0_xxyyzz_zz, g_x_0_x_0_xxyzz_xx, g_x_0_x_0_xxyzz_xxy, g_x_0_x_0_xxyzz_xy, g_x_0_x_0_xxyzz_xyy, g_x_0_x_0_xxyzz_xyz, g_x_0_x_0_xxyzz_xz, g_x_0_x_0_xxyzz_yy, g_x_0_x_0_xxyzz_yyy, g_x_0_x_0_xxyzz_yyz, g_x_0_x_0_xxyzz_yz, g_x_0_x_0_xxyzz_yzz, g_x_0_x_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyyzz_xx[k] = -g_x_0_x_0_xxyzz_xx[k] * ab_y + g_x_0_x_0_xxyzz_xxy[k];

                g_x_0_x_0_xxyyzz_xy[k] = -g_x_0_x_0_xxyzz_xy[k] * ab_y + g_x_0_x_0_xxyzz_xyy[k];

                g_x_0_x_0_xxyyzz_xz[k] = -g_x_0_x_0_xxyzz_xz[k] * ab_y + g_x_0_x_0_xxyzz_xyz[k];

                g_x_0_x_0_xxyyzz_yy[k] = -g_x_0_x_0_xxyzz_yy[k] * ab_y + g_x_0_x_0_xxyzz_yyy[k];

                g_x_0_x_0_xxyyzz_yz[k] = -g_x_0_x_0_xxyzz_yz[k] * ab_y + g_x_0_x_0_xxyzz_yyz[k];

                g_x_0_x_0_xxyyzz_zz[k] = -g_x_0_x_0_xxyzz_zz[k] * ab_y + g_x_0_x_0_xxyzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyzzz_xx = cbuffer.data(id_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzzz_xy = cbuffer.data(id_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzzz_xz = cbuffer.data(id_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzzz_yy = cbuffer.data(id_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzzz_yz = cbuffer.data(id_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzzz_zz = cbuffer.data(id_geom_1010_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyzzz_xx, g_x_0_x_0_xxyzzz_xy, g_x_0_x_0_xxyzzz_xz, g_x_0_x_0_xxyzzz_yy, g_x_0_x_0_xxyzzz_yz, g_x_0_x_0_xxyzzz_zz, g_x_0_x_0_xxzzz_xx, g_x_0_x_0_xxzzz_xxy, g_x_0_x_0_xxzzz_xy, g_x_0_x_0_xxzzz_xyy, g_x_0_x_0_xxzzz_xyz, g_x_0_x_0_xxzzz_xz, g_x_0_x_0_xxzzz_yy, g_x_0_x_0_xxzzz_yyy, g_x_0_x_0_xxzzz_yyz, g_x_0_x_0_xxzzz_yz, g_x_0_x_0_xxzzz_yzz, g_x_0_x_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyzzz_xx[k] = -g_x_0_x_0_xxzzz_xx[k] * ab_y + g_x_0_x_0_xxzzz_xxy[k];

                g_x_0_x_0_xxyzzz_xy[k] = -g_x_0_x_0_xxzzz_xy[k] * ab_y + g_x_0_x_0_xxzzz_xyy[k];

                g_x_0_x_0_xxyzzz_xz[k] = -g_x_0_x_0_xxzzz_xz[k] * ab_y + g_x_0_x_0_xxzzz_xyz[k];

                g_x_0_x_0_xxyzzz_yy[k] = -g_x_0_x_0_xxzzz_yy[k] * ab_y + g_x_0_x_0_xxzzz_yyy[k];

                g_x_0_x_0_xxyzzz_yz[k] = -g_x_0_x_0_xxzzz_yz[k] * ab_y + g_x_0_x_0_xxzzz_yyz[k];

                g_x_0_x_0_xxyzzz_zz[k] = -g_x_0_x_0_xxzzz_zz[k] * ab_y + g_x_0_x_0_xxzzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxzzzz_xx = cbuffer.data(id_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzzz_xy = cbuffer.data(id_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzzz_xz = cbuffer.data(id_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzzz_yy = cbuffer.data(id_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzzz_yz = cbuffer.data(id_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzzz_zz = cbuffer.data(id_geom_1010_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxzzz_xx, g_x_0_x_0_xxzzz_xxz, g_x_0_x_0_xxzzz_xy, g_x_0_x_0_xxzzz_xyz, g_x_0_x_0_xxzzz_xz, g_x_0_x_0_xxzzz_xzz, g_x_0_x_0_xxzzz_yy, g_x_0_x_0_xxzzz_yyz, g_x_0_x_0_xxzzz_yz, g_x_0_x_0_xxzzz_yzz, g_x_0_x_0_xxzzz_zz, g_x_0_x_0_xxzzz_zzz, g_x_0_x_0_xxzzzz_xx, g_x_0_x_0_xxzzzz_xy, g_x_0_x_0_xxzzzz_xz, g_x_0_x_0_xxzzzz_yy, g_x_0_x_0_xxzzzz_yz, g_x_0_x_0_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxzzzz_xx[k] = -g_x_0_x_0_xxzzz_xx[k] * ab_z + g_x_0_x_0_xxzzz_xxz[k];

                g_x_0_x_0_xxzzzz_xy[k] = -g_x_0_x_0_xxzzz_xy[k] * ab_z + g_x_0_x_0_xxzzz_xyz[k];

                g_x_0_x_0_xxzzzz_xz[k] = -g_x_0_x_0_xxzzz_xz[k] * ab_z + g_x_0_x_0_xxzzz_xzz[k];

                g_x_0_x_0_xxzzzz_yy[k] = -g_x_0_x_0_xxzzz_yy[k] * ab_z + g_x_0_x_0_xxzzz_yyz[k];

                g_x_0_x_0_xxzzzz_yz[k] = -g_x_0_x_0_xxzzz_yz[k] * ab_z + g_x_0_x_0_xxzzz_yzz[k];

                g_x_0_x_0_xxzzzz_zz[k] = -g_x_0_x_0_xxzzz_zz[k] * ab_z + g_x_0_x_0_xxzzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyyyy_xx = cbuffer.data(id_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyy_xy = cbuffer.data(id_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyy_xz = cbuffer.data(id_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyy_yy = cbuffer.data(id_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyy_yz = cbuffer.data(id_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyy_zz = cbuffer.data(id_geom_1010_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyyy_xx, g_x_0_x_0_xyyyy_xxy, g_x_0_x_0_xyyyy_xy, g_x_0_x_0_xyyyy_xyy, g_x_0_x_0_xyyyy_xyz, g_x_0_x_0_xyyyy_xz, g_x_0_x_0_xyyyy_yy, g_x_0_x_0_xyyyy_yyy, g_x_0_x_0_xyyyy_yyz, g_x_0_x_0_xyyyy_yz, g_x_0_x_0_xyyyy_yzz, g_x_0_x_0_xyyyy_zz, g_x_0_x_0_xyyyyy_xx, g_x_0_x_0_xyyyyy_xy, g_x_0_x_0_xyyyyy_xz, g_x_0_x_0_xyyyyy_yy, g_x_0_x_0_xyyyyy_yz, g_x_0_x_0_xyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyyyy_xx[k] = -g_x_0_x_0_xyyyy_xx[k] * ab_y + g_x_0_x_0_xyyyy_xxy[k];

                g_x_0_x_0_xyyyyy_xy[k] = -g_x_0_x_0_xyyyy_xy[k] * ab_y + g_x_0_x_0_xyyyy_xyy[k];

                g_x_0_x_0_xyyyyy_xz[k] = -g_x_0_x_0_xyyyy_xz[k] * ab_y + g_x_0_x_0_xyyyy_xyz[k];

                g_x_0_x_0_xyyyyy_yy[k] = -g_x_0_x_0_xyyyy_yy[k] * ab_y + g_x_0_x_0_xyyyy_yyy[k];

                g_x_0_x_0_xyyyyy_yz[k] = -g_x_0_x_0_xyyyy_yz[k] * ab_y + g_x_0_x_0_xyyyy_yyz[k];

                g_x_0_x_0_xyyyyy_zz[k] = -g_x_0_x_0_xyyyy_zz[k] * ab_y + g_x_0_x_0_xyyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyyyz_xx = cbuffer.data(id_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyz_xy = cbuffer.data(id_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyz_xz = cbuffer.data(id_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyz_yy = cbuffer.data(id_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyz_yz = cbuffer.data(id_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyz_zz = cbuffer.data(id_geom_1010_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyyyz_xx, g_x_0_x_0_xyyyyz_xy, g_x_0_x_0_xyyyyz_xz, g_x_0_x_0_xyyyyz_yy, g_x_0_x_0_xyyyyz_yz, g_x_0_x_0_xyyyyz_zz, g_x_0_x_0_xyyyz_xx, g_x_0_x_0_xyyyz_xxy, g_x_0_x_0_xyyyz_xy, g_x_0_x_0_xyyyz_xyy, g_x_0_x_0_xyyyz_xyz, g_x_0_x_0_xyyyz_xz, g_x_0_x_0_xyyyz_yy, g_x_0_x_0_xyyyz_yyy, g_x_0_x_0_xyyyz_yyz, g_x_0_x_0_xyyyz_yz, g_x_0_x_0_xyyyz_yzz, g_x_0_x_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyyyz_xx[k] = -g_x_0_x_0_xyyyz_xx[k] * ab_y + g_x_0_x_0_xyyyz_xxy[k];

                g_x_0_x_0_xyyyyz_xy[k] = -g_x_0_x_0_xyyyz_xy[k] * ab_y + g_x_0_x_0_xyyyz_xyy[k];

                g_x_0_x_0_xyyyyz_xz[k] = -g_x_0_x_0_xyyyz_xz[k] * ab_y + g_x_0_x_0_xyyyz_xyz[k];

                g_x_0_x_0_xyyyyz_yy[k] = -g_x_0_x_0_xyyyz_yy[k] * ab_y + g_x_0_x_0_xyyyz_yyy[k];

                g_x_0_x_0_xyyyyz_yz[k] = -g_x_0_x_0_xyyyz_yz[k] * ab_y + g_x_0_x_0_xyyyz_yyz[k];

                g_x_0_x_0_xyyyyz_zz[k] = -g_x_0_x_0_xyyyz_zz[k] * ab_y + g_x_0_x_0_xyyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyyzz_xx = cbuffer.data(id_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyzz_xy = cbuffer.data(id_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyzz_xz = cbuffer.data(id_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyzz_yy = cbuffer.data(id_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyzz_yz = cbuffer.data(id_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyzz_zz = cbuffer.data(id_geom_1010_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyyzz_xx, g_x_0_x_0_xyyyzz_xy, g_x_0_x_0_xyyyzz_xz, g_x_0_x_0_xyyyzz_yy, g_x_0_x_0_xyyyzz_yz, g_x_0_x_0_xyyyzz_zz, g_x_0_x_0_xyyzz_xx, g_x_0_x_0_xyyzz_xxy, g_x_0_x_0_xyyzz_xy, g_x_0_x_0_xyyzz_xyy, g_x_0_x_0_xyyzz_xyz, g_x_0_x_0_xyyzz_xz, g_x_0_x_0_xyyzz_yy, g_x_0_x_0_xyyzz_yyy, g_x_0_x_0_xyyzz_yyz, g_x_0_x_0_xyyzz_yz, g_x_0_x_0_xyyzz_yzz, g_x_0_x_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyyzz_xx[k] = -g_x_0_x_0_xyyzz_xx[k] * ab_y + g_x_0_x_0_xyyzz_xxy[k];

                g_x_0_x_0_xyyyzz_xy[k] = -g_x_0_x_0_xyyzz_xy[k] * ab_y + g_x_0_x_0_xyyzz_xyy[k];

                g_x_0_x_0_xyyyzz_xz[k] = -g_x_0_x_0_xyyzz_xz[k] * ab_y + g_x_0_x_0_xyyzz_xyz[k];

                g_x_0_x_0_xyyyzz_yy[k] = -g_x_0_x_0_xyyzz_yy[k] * ab_y + g_x_0_x_0_xyyzz_yyy[k];

                g_x_0_x_0_xyyyzz_yz[k] = -g_x_0_x_0_xyyzz_yz[k] * ab_y + g_x_0_x_0_xyyzz_yyz[k];

                g_x_0_x_0_xyyyzz_zz[k] = -g_x_0_x_0_xyyzz_zz[k] * ab_y + g_x_0_x_0_xyyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyzzz_xx = cbuffer.data(id_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzzz_xy = cbuffer.data(id_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzzz_xz = cbuffer.data(id_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzzz_yy = cbuffer.data(id_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzzz_yz = cbuffer.data(id_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzzz_zz = cbuffer.data(id_geom_1010_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyzzz_xx, g_x_0_x_0_xyyzzz_xy, g_x_0_x_0_xyyzzz_xz, g_x_0_x_0_xyyzzz_yy, g_x_0_x_0_xyyzzz_yz, g_x_0_x_0_xyyzzz_zz, g_x_0_x_0_xyzzz_xx, g_x_0_x_0_xyzzz_xxy, g_x_0_x_0_xyzzz_xy, g_x_0_x_0_xyzzz_xyy, g_x_0_x_0_xyzzz_xyz, g_x_0_x_0_xyzzz_xz, g_x_0_x_0_xyzzz_yy, g_x_0_x_0_xyzzz_yyy, g_x_0_x_0_xyzzz_yyz, g_x_0_x_0_xyzzz_yz, g_x_0_x_0_xyzzz_yzz, g_x_0_x_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyzzz_xx[k] = -g_x_0_x_0_xyzzz_xx[k] * ab_y + g_x_0_x_0_xyzzz_xxy[k];

                g_x_0_x_0_xyyzzz_xy[k] = -g_x_0_x_0_xyzzz_xy[k] * ab_y + g_x_0_x_0_xyzzz_xyy[k];

                g_x_0_x_0_xyyzzz_xz[k] = -g_x_0_x_0_xyzzz_xz[k] * ab_y + g_x_0_x_0_xyzzz_xyz[k];

                g_x_0_x_0_xyyzzz_yy[k] = -g_x_0_x_0_xyzzz_yy[k] * ab_y + g_x_0_x_0_xyzzz_yyy[k];

                g_x_0_x_0_xyyzzz_yz[k] = -g_x_0_x_0_xyzzz_yz[k] * ab_y + g_x_0_x_0_xyzzz_yyz[k];

                g_x_0_x_0_xyyzzz_zz[k] = -g_x_0_x_0_xyzzz_zz[k] * ab_y + g_x_0_x_0_xyzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyzzzz_xx = cbuffer.data(id_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzzz_xy = cbuffer.data(id_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzzz_xz = cbuffer.data(id_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzzz_yy = cbuffer.data(id_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzzz_yz = cbuffer.data(id_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzzz_zz = cbuffer.data(id_geom_1010_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyzzzz_xx, g_x_0_x_0_xyzzzz_xy, g_x_0_x_0_xyzzzz_xz, g_x_0_x_0_xyzzzz_yy, g_x_0_x_0_xyzzzz_yz, g_x_0_x_0_xyzzzz_zz, g_x_0_x_0_xzzzz_xx, g_x_0_x_0_xzzzz_xxy, g_x_0_x_0_xzzzz_xy, g_x_0_x_0_xzzzz_xyy, g_x_0_x_0_xzzzz_xyz, g_x_0_x_0_xzzzz_xz, g_x_0_x_0_xzzzz_yy, g_x_0_x_0_xzzzz_yyy, g_x_0_x_0_xzzzz_yyz, g_x_0_x_0_xzzzz_yz, g_x_0_x_0_xzzzz_yzz, g_x_0_x_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyzzzz_xx[k] = -g_x_0_x_0_xzzzz_xx[k] * ab_y + g_x_0_x_0_xzzzz_xxy[k];

                g_x_0_x_0_xyzzzz_xy[k] = -g_x_0_x_0_xzzzz_xy[k] * ab_y + g_x_0_x_0_xzzzz_xyy[k];

                g_x_0_x_0_xyzzzz_xz[k] = -g_x_0_x_0_xzzzz_xz[k] * ab_y + g_x_0_x_0_xzzzz_xyz[k];

                g_x_0_x_0_xyzzzz_yy[k] = -g_x_0_x_0_xzzzz_yy[k] * ab_y + g_x_0_x_0_xzzzz_yyy[k];

                g_x_0_x_0_xyzzzz_yz[k] = -g_x_0_x_0_xzzzz_yz[k] * ab_y + g_x_0_x_0_xzzzz_yyz[k];

                g_x_0_x_0_xyzzzz_zz[k] = -g_x_0_x_0_xzzzz_zz[k] * ab_y + g_x_0_x_0_xzzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xzzzzz_xx = cbuffer.data(id_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzzz_xy = cbuffer.data(id_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzzz_xz = cbuffer.data(id_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzzz_yy = cbuffer.data(id_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzzz_yz = cbuffer.data(id_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzzz_zz = cbuffer.data(id_geom_1010_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xzzzz_xx, g_x_0_x_0_xzzzz_xxz, g_x_0_x_0_xzzzz_xy, g_x_0_x_0_xzzzz_xyz, g_x_0_x_0_xzzzz_xz, g_x_0_x_0_xzzzz_xzz, g_x_0_x_0_xzzzz_yy, g_x_0_x_0_xzzzz_yyz, g_x_0_x_0_xzzzz_yz, g_x_0_x_0_xzzzz_yzz, g_x_0_x_0_xzzzz_zz, g_x_0_x_0_xzzzz_zzz, g_x_0_x_0_xzzzzz_xx, g_x_0_x_0_xzzzzz_xy, g_x_0_x_0_xzzzzz_xz, g_x_0_x_0_xzzzzz_yy, g_x_0_x_0_xzzzzz_yz, g_x_0_x_0_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xzzzzz_xx[k] = -g_x_0_x_0_xzzzz_xx[k] * ab_z + g_x_0_x_0_xzzzz_xxz[k];

                g_x_0_x_0_xzzzzz_xy[k] = -g_x_0_x_0_xzzzz_xy[k] * ab_z + g_x_0_x_0_xzzzz_xyz[k];

                g_x_0_x_0_xzzzzz_xz[k] = -g_x_0_x_0_xzzzz_xz[k] * ab_z + g_x_0_x_0_xzzzz_xzz[k];

                g_x_0_x_0_xzzzzz_yy[k] = -g_x_0_x_0_xzzzz_yy[k] * ab_z + g_x_0_x_0_xzzzz_yyz[k];

                g_x_0_x_0_xzzzzz_yz[k] = -g_x_0_x_0_xzzzz_yz[k] * ab_z + g_x_0_x_0_xzzzz_yzz[k];

                g_x_0_x_0_xzzzzz_zz[k] = -g_x_0_x_0_xzzzz_zz[k] * ab_z + g_x_0_x_0_xzzzz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyyyy_xx = cbuffer.data(id_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyy_xy = cbuffer.data(id_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyy_xz = cbuffer.data(id_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyy_yy = cbuffer.data(id_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyy_yz = cbuffer.data(id_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyy_zz = cbuffer.data(id_geom_1010_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyyy_xx, g_x_0_x_0_yyyyy_xxy, g_x_0_x_0_yyyyy_xy, g_x_0_x_0_yyyyy_xyy, g_x_0_x_0_yyyyy_xyz, g_x_0_x_0_yyyyy_xz, g_x_0_x_0_yyyyy_yy, g_x_0_x_0_yyyyy_yyy, g_x_0_x_0_yyyyy_yyz, g_x_0_x_0_yyyyy_yz, g_x_0_x_0_yyyyy_yzz, g_x_0_x_0_yyyyy_zz, g_x_0_x_0_yyyyyy_xx, g_x_0_x_0_yyyyyy_xy, g_x_0_x_0_yyyyyy_xz, g_x_0_x_0_yyyyyy_yy, g_x_0_x_0_yyyyyy_yz, g_x_0_x_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyyyy_xx[k] = -g_x_0_x_0_yyyyy_xx[k] * ab_y + g_x_0_x_0_yyyyy_xxy[k];

                g_x_0_x_0_yyyyyy_xy[k] = -g_x_0_x_0_yyyyy_xy[k] * ab_y + g_x_0_x_0_yyyyy_xyy[k];

                g_x_0_x_0_yyyyyy_xz[k] = -g_x_0_x_0_yyyyy_xz[k] * ab_y + g_x_0_x_0_yyyyy_xyz[k];

                g_x_0_x_0_yyyyyy_yy[k] = -g_x_0_x_0_yyyyy_yy[k] * ab_y + g_x_0_x_0_yyyyy_yyy[k];

                g_x_0_x_0_yyyyyy_yz[k] = -g_x_0_x_0_yyyyy_yz[k] * ab_y + g_x_0_x_0_yyyyy_yyz[k];

                g_x_0_x_0_yyyyyy_zz[k] = -g_x_0_x_0_yyyyy_zz[k] * ab_y + g_x_0_x_0_yyyyy_yzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyyyz_xx = cbuffer.data(id_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyz_xy = cbuffer.data(id_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyz_xz = cbuffer.data(id_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyz_yy = cbuffer.data(id_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyz_yz = cbuffer.data(id_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyz_zz = cbuffer.data(id_geom_1010_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyyyz_xx, g_x_0_x_0_yyyyyz_xy, g_x_0_x_0_yyyyyz_xz, g_x_0_x_0_yyyyyz_yy, g_x_0_x_0_yyyyyz_yz, g_x_0_x_0_yyyyyz_zz, g_x_0_x_0_yyyyz_xx, g_x_0_x_0_yyyyz_xxy, g_x_0_x_0_yyyyz_xy, g_x_0_x_0_yyyyz_xyy, g_x_0_x_0_yyyyz_xyz, g_x_0_x_0_yyyyz_xz, g_x_0_x_0_yyyyz_yy, g_x_0_x_0_yyyyz_yyy, g_x_0_x_0_yyyyz_yyz, g_x_0_x_0_yyyyz_yz, g_x_0_x_0_yyyyz_yzz, g_x_0_x_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyyyz_xx[k] = -g_x_0_x_0_yyyyz_xx[k] * ab_y + g_x_0_x_0_yyyyz_xxy[k];

                g_x_0_x_0_yyyyyz_xy[k] = -g_x_0_x_0_yyyyz_xy[k] * ab_y + g_x_0_x_0_yyyyz_xyy[k];

                g_x_0_x_0_yyyyyz_xz[k] = -g_x_0_x_0_yyyyz_xz[k] * ab_y + g_x_0_x_0_yyyyz_xyz[k];

                g_x_0_x_0_yyyyyz_yy[k] = -g_x_0_x_0_yyyyz_yy[k] * ab_y + g_x_0_x_0_yyyyz_yyy[k];

                g_x_0_x_0_yyyyyz_yz[k] = -g_x_0_x_0_yyyyz_yz[k] * ab_y + g_x_0_x_0_yyyyz_yyz[k];

                g_x_0_x_0_yyyyyz_zz[k] = -g_x_0_x_0_yyyyz_zz[k] * ab_y + g_x_0_x_0_yyyyz_yzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyyzz_xx = cbuffer.data(id_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyzz_xy = cbuffer.data(id_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyzz_xz = cbuffer.data(id_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyzz_yy = cbuffer.data(id_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyzz_yz = cbuffer.data(id_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyzz_zz = cbuffer.data(id_geom_1010_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyyzz_xx, g_x_0_x_0_yyyyzz_xy, g_x_0_x_0_yyyyzz_xz, g_x_0_x_0_yyyyzz_yy, g_x_0_x_0_yyyyzz_yz, g_x_0_x_0_yyyyzz_zz, g_x_0_x_0_yyyzz_xx, g_x_0_x_0_yyyzz_xxy, g_x_0_x_0_yyyzz_xy, g_x_0_x_0_yyyzz_xyy, g_x_0_x_0_yyyzz_xyz, g_x_0_x_0_yyyzz_xz, g_x_0_x_0_yyyzz_yy, g_x_0_x_0_yyyzz_yyy, g_x_0_x_0_yyyzz_yyz, g_x_0_x_0_yyyzz_yz, g_x_0_x_0_yyyzz_yzz, g_x_0_x_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyyzz_xx[k] = -g_x_0_x_0_yyyzz_xx[k] * ab_y + g_x_0_x_0_yyyzz_xxy[k];

                g_x_0_x_0_yyyyzz_xy[k] = -g_x_0_x_0_yyyzz_xy[k] * ab_y + g_x_0_x_0_yyyzz_xyy[k];

                g_x_0_x_0_yyyyzz_xz[k] = -g_x_0_x_0_yyyzz_xz[k] * ab_y + g_x_0_x_0_yyyzz_xyz[k];

                g_x_0_x_0_yyyyzz_yy[k] = -g_x_0_x_0_yyyzz_yy[k] * ab_y + g_x_0_x_0_yyyzz_yyy[k];

                g_x_0_x_0_yyyyzz_yz[k] = -g_x_0_x_0_yyyzz_yz[k] * ab_y + g_x_0_x_0_yyyzz_yyz[k];

                g_x_0_x_0_yyyyzz_zz[k] = -g_x_0_x_0_yyyzz_zz[k] * ab_y + g_x_0_x_0_yyyzz_yzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyzzz_xx = cbuffer.data(id_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzzz_xy = cbuffer.data(id_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzzz_xz = cbuffer.data(id_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzzz_yy = cbuffer.data(id_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzzz_yz = cbuffer.data(id_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzzz_zz = cbuffer.data(id_geom_1010_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyzzz_xx, g_x_0_x_0_yyyzzz_xy, g_x_0_x_0_yyyzzz_xz, g_x_0_x_0_yyyzzz_yy, g_x_0_x_0_yyyzzz_yz, g_x_0_x_0_yyyzzz_zz, g_x_0_x_0_yyzzz_xx, g_x_0_x_0_yyzzz_xxy, g_x_0_x_0_yyzzz_xy, g_x_0_x_0_yyzzz_xyy, g_x_0_x_0_yyzzz_xyz, g_x_0_x_0_yyzzz_xz, g_x_0_x_0_yyzzz_yy, g_x_0_x_0_yyzzz_yyy, g_x_0_x_0_yyzzz_yyz, g_x_0_x_0_yyzzz_yz, g_x_0_x_0_yyzzz_yzz, g_x_0_x_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyzzz_xx[k] = -g_x_0_x_0_yyzzz_xx[k] * ab_y + g_x_0_x_0_yyzzz_xxy[k];

                g_x_0_x_0_yyyzzz_xy[k] = -g_x_0_x_0_yyzzz_xy[k] * ab_y + g_x_0_x_0_yyzzz_xyy[k];

                g_x_0_x_0_yyyzzz_xz[k] = -g_x_0_x_0_yyzzz_xz[k] * ab_y + g_x_0_x_0_yyzzz_xyz[k];

                g_x_0_x_0_yyyzzz_yy[k] = -g_x_0_x_0_yyzzz_yy[k] * ab_y + g_x_0_x_0_yyzzz_yyy[k];

                g_x_0_x_0_yyyzzz_yz[k] = -g_x_0_x_0_yyzzz_yz[k] * ab_y + g_x_0_x_0_yyzzz_yyz[k];

                g_x_0_x_0_yyyzzz_zz[k] = -g_x_0_x_0_yyzzz_zz[k] * ab_y + g_x_0_x_0_yyzzz_yzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyzzzz_xx = cbuffer.data(id_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzzz_xy = cbuffer.data(id_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzzz_xz = cbuffer.data(id_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzzz_yy = cbuffer.data(id_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzzz_yz = cbuffer.data(id_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzzz_zz = cbuffer.data(id_geom_1010_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyzzzz_xx, g_x_0_x_0_yyzzzz_xy, g_x_0_x_0_yyzzzz_xz, g_x_0_x_0_yyzzzz_yy, g_x_0_x_0_yyzzzz_yz, g_x_0_x_0_yyzzzz_zz, g_x_0_x_0_yzzzz_xx, g_x_0_x_0_yzzzz_xxy, g_x_0_x_0_yzzzz_xy, g_x_0_x_0_yzzzz_xyy, g_x_0_x_0_yzzzz_xyz, g_x_0_x_0_yzzzz_xz, g_x_0_x_0_yzzzz_yy, g_x_0_x_0_yzzzz_yyy, g_x_0_x_0_yzzzz_yyz, g_x_0_x_0_yzzzz_yz, g_x_0_x_0_yzzzz_yzz, g_x_0_x_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyzzzz_xx[k] = -g_x_0_x_0_yzzzz_xx[k] * ab_y + g_x_0_x_0_yzzzz_xxy[k];

                g_x_0_x_0_yyzzzz_xy[k] = -g_x_0_x_0_yzzzz_xy[k] * ab_y + g_x_0_x_0_yzzzz_xyy[k];

                g_x_0_x_0_yyzzzz_xz[k] = -g_x_0_x_0_yzzzz_xz[k] * ab_y + g_x_0_x_0_yzzzz_xyz[k];

                g_x_0_x_0_yyzzzz_yy[k] = -g_x_0_x_0_yzzzz_yy[k] * ab_y + g_x_0_x_0_yzzzz_yyy[k];

                g_x_0_x_0_yyzzzz_yz[k] = -g_x_0_x_0_yzzzz_yz[k] * ab_y + g_x_0_x_0_yzzzz_yyz[k];

                g_x_0_x_0_yyzzzz_zz[k] = -g_x_0_x_0_yzzzz_zz[k] * ab_y + g_x_0_x_0_yzzzz_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yzzzzz_xx = cbuffer.data(id_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzzz_xy = cbuffer.data(id_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzzz_xz = cbuffer.data(id_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzzz_yy = cbuffer.data(id_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzzz_yz = cbuffer.data(id_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzzz_zz = cbuffer.data(id_geom_1010_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yzzzzz_xx, g_x_0_x_0_yzzzzz_xy, g_x_0_x_0_yzzzzz_xz, g_x_0_x_0_yzzzzz_yy, g_x_0_x_0_yzzzzz_yz, g_x_0_x_0_yzzzzz_zz, g_x_0_x_0_zzzzz_xx, g_x_0_x_0_zzzzz_xxy, g_x_0_x_0_zzzzz_xy, g_x_0_x_0_zzzzz_xyy, g_x_0_x_0_zzzzz_xyz, g_x_0_x_0_zzzzz_xz, g_x_0_x_0_zzzzz_yy, g_x_0_x_0_zzzzz_yyy, g_x_0_x_0_zzzzz_yyz, g_x_0_x_0_zzzzz_yz, g_x_0_x_0_zzzzz_yzz, g_x_0_x_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yzzzzz_xx[k] = -g_x_0_x_0_zzzzz_xx[k] * ab_y + g_x_0_x_0_zzzzz_xxy[k];

                g_x_0_x_0_yzzzzz_xy[k] = -g_x_0_x_0_zzzzz_xy[k] * ab_y + g_x_0_x_0_zzzzz_xyy[k];

                g_x_0_x_0_yzzzzz_xz[k] = -g_x_0_x_0_zzzzz_xz[k] * ab_y + g_x_0_x_0_zzzzz_xyz[k];

                g_x_0_x_0_yzzzzz_yy[k] = -g_x_0_x_0_zzzzz_yy[k] * ab_y + g_x_0_x_0_zzzzz_yyy[k];

                g_x_0_x_0_yzzzzz_yz[k] = -g_x_0_x_0_zzzzz_yz[k] * ab_y + g_x_0_x_0_zzzzz_yyz[k];

                g_x_0_x_0_yzzzzz_zz[k] = -g_x_0_x_0_zzzzz_zz[k] * ab_y + g_x_0_x_0_zzzzz_yzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_zzzzzz_xx = cbuffer.data(id_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzzz_xy = cbuffer.data(id_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzzz_xz = cbuffer.data(id_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzzz_yy = cbuffer.data(id_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzzz_yz = cbuffer.data(id_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzzz_zz = cbuffer.data(id_geom_1010_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_zzzzz_xx, g_x_0_x_0_zzzzz_xxz, g_x_0_x_0_zzzzz_xy, g_x_0_x_0_zzzzz_xyz, g_x_0_x_0_zzzzz_xz, g_x_0_x_0_zzzzz_xzz, g_x_0_x_0_zzzzz_yy, g_x_0_x_0_zzzzz_yyz, g_x_0_x_0_zzzzz_yz, g_x_0_x_0_zzzzz_yzz, g_x_0_x_0_zzzzz_zz, g_x_0_x_0_zzzzz_zzz, g_x_0_x_0_zzzzzz_xx, g_x_0_x_0_zzzzzz_xy, g_x_0_x_0_zzzzzz_xz, g_x_0_x_0_zzzzzz_yy, g_x_0_x_0_zzzzzz_yz, g_x_0_x_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_zzzzzz_xx[k] = -g_x_0_x_0_zzzzz_xx[k] * ab_z + g_x_0_x_0_zzzzz_xxz[k];

                g_x_0_x_0_zzzzzz_xy[k] = -g_x_0_x_0_zzzzz_xy[k] * ab_z + g_x_0_x_0_zzzzz_xyz[k];

                g_x_0_x_0_zzzzzz_xz[k] = -g_x_0_x_0_zzzzz_xz[k] * ab_z + g_x_0_x_0_zzzzz_xzz[k];

                g_x_0_x_0_zzzzzz_yy[k] = -g_x_0_x_0_zzzzz_yy[k] * ab_z + g_x_0_x_0_zzzzz_yyz[k];

                g_x_0_x_0_zzzzzz_yz[k] = -g_x_0_x_0_zzzzz_yz[k] * ab_z + g_x_0_x_0_zzzzz_yzz[k];

                g_x_0_x_0_zzzzzz_zz[k] = -g_x_0_x_0_zzzzz_zz[k] * ab_z + g_x_0_x_0_zzzzz_zzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxxx_xx = cbuffer.data(id_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxx_xy = cbuffer.data(id_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxx_xz = cbuffer.data(id_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxx_yy = cbuffer.data(id_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxx_yz = cbuffer.data(id_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxx_zz = cbuffer.data(id_geom_1010_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_xxxxx_xx, g_0_0_y_0_xxxxx_xy, g_0_0_y_0_xxxxx_xz, g_0_0_y_0_xxxxx_yy, g_0_0_y_0_xxxxx_yz, g_0_0_y_0_xxxxx_zz, g_x_0_y_0_xxxxx_xx, g_x_0_y_0_xxxxx_xxx, g_x_0_y_0_xxxxx_xxy, g_x_0_y_0_xxxxx_xxz, g_x_0_y_0_xxxxx_xy, g_x_0_y_0_xxxxx_xyy, g_x_0_y_0_xxxxx_xyz, g_x_0_y_0_xxxxx_xz, g_x_0_y_0_xxxxx_xzz, g_x_0_y_0_xxxxx_yy, g_x_0_y_0_xxxxx_yz, g_x_0_y_0_xxxxx_zz, g_x_0_y_0_xxxxxx_xx, g_x_0_y_0_xxxxxx_xy, g_x_0_y_0_xxxxxx_xz, g_x_0_y_0_xxxxxx_yy, g_x_0_y_0_xxxxxx_yz, g_x_0_y_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxxx_xx[k] = -g_0_0_y_0_xxxxx_xx[k] - g_x_0_y_0_xxxxx_xx[k] * ab_x + g_x_0_y_0_xxxxx_xxx[k];

                g_x_0_y_0_xxxxxx_xy[k] = -g_0_0_y_0_xxxxx_xy[k] - g_x_0_y_0_xxxxx_xy[k] * ab_x + g_x_0_y_0_xxxxx_xxy[k];

                g_x_0_y_0_xxxxxx_xz[k] = -g_0_0_y_0_xxxxx_xz[k] - g_x_0_y_0_xxxxx_xz[k] * ab_x + g_x_0_y_0_xxxxx_xxz[k];

                g_x_0_y_0_xxxxxx_yy[k] = -g_0_0_y_0_xxxxx_yy[k] - g_x_0_y_0_xxxxx_yy[k] * ab_x + g_x_0_y_0_xxxxx_xyy[k];

                g_x_0_y_0_xxxxxx_yz[k] = -g_0_0_y_0_xxxxx_yz[k] - g_x_0_y_0_xxxxx_yz[k] * ab_x + g_x_0_y_0_xxxxx_xyz[k];

                g_x_0_y_0_xxxxxx_zz[k] = -g_0_0_y_0_xxxxx_zz[k] - g_x_0_y_0_xxxxx_zz[k] * ab_x + g_x_0_y_0_xxxxx_xzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxxy_xx = cbuffer.data(id_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxy_xy = cbuffer.data(id_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxy_xz = cbuffer.data(id_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxy_yy = cbuffer.data(id_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxy_yz = cbuffer.data(id_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxy_zz = cbuffer.data(id_geom_1010_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxx_xx, g_x_0_y_0_xxxxx_xxy, g_x_0_y_0_xxxxx_xy, g_x_0_y_0_xxxxx_xyy, g_x_0_y_0_xxxxx_xyz, g_x_0_y_0_xxxxx_xz, g_x_0_y_0_xxxxx_yy, g_x_0_y_0_xxxxx_yyy, g_x_0_y_0_xxxxx_yyz, g_x_0_y_0_xxxxx_yz, g_x_0_y_0_xxxxx_yzz, g_x_0_y_0_xxxxx_zz, g_x_0_y_0_xxxxxy_xx, g_x_0_y_0_xxxxxy_xy, g_x_0_y_0_xxxxxy_xz, g_x_0_y_0_xxxxxy_yy, g_x_0_y_0_xxxxxy_yz, g_x_0_y_0_xxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxxy_xx[k] = -g_x_0_y_0_xxxxx_xx[k] * ab_y + g_x_0_y_0_xxxxx_xxy[k];

                g_x_0_y_0_xxxxxy_xy[k] = -g_x_0_y_0_xxxxx_xy[k] * ab_y + g_x_0_y_0_xxxxx_xyy[k];

                g_x_0_y_0_xxxxxy_xz[k] = -g_x_0_y_0_xxxxx_xz[k] * ab_y + g_x_0_y_0_xxxxx_xyz[k];

                g_x_0_y_0_xxxxxy_yy[k] = -g_x_0_y_0_xxxxx_yy[k] * ab_y + g_x_0_y_0_xxxxx_yyy[k];

                g_x_0_y_0_xxxxxy_yz[k] = -g_x_0_y_0_xxxxx_yz[k] * ab_y + g_x_0_y_0_xxxxx_yyz[k];

                g_x_0_y_0_xxxxxy_zz[k] = -g_x_0_y_0_xxxxx_zz[k] * ab_y + g_x_0_y_0_xxxxx_yzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxxz_xx = cbuffer.data(id_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxz_xy = cbuffer.data(id_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxz_xz = cbuffer.data(id_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxz_yy = cbuffer.data(id_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxz_yz = cbuffer.data(id_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxz_zz = cbuffer.data(id_geom_1010_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxx_xx, g_x_0_y_0_xxxxx_xxz, g_x_0_y_0_xxxxx_xy, g_x_0_y_0_xxxxx_xyz, g_x_0_y_0_xxxxx_xz, g_x_0_y_0_xxxxx_xzz, g_x_0_y_0_xxxxx_yy, g_x_0_y_0_xxxxx_yyz, g_x_0_y_0_xxxxx_yz, g_x_0_y_0_xxxxx_yzz, g_x_0_y_0_xxxxx_zz, g_x_0_y_0_xxxxx_zzz, g_x_0_y_0_xxxxxz_xx, g_x_0_y_0_xxxxxz_xy, g_x_0_y_0_xxxxxz_xz, g_x_0_y_0_xxxxxz_yy, g_x_0_y_0_xxxxxz_yz, g_x_0_y_0_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxxz_xx[k] = -g_x_0_y_0_xxxxx_xx[k] * ab_z + g_x_0_y_0_xxxxx_xxz[k];

                g_x_0_y_0_xxxxxz_xy[k] = -g_x_0_y_0_xxxxx_xy[k] * ab_z + g_x_0_y_0_xxxxx_xyz[k];

                g_x_0_y_0_xxxxxz_xz[k] = -g_x_0_y_0_xxxxx_xz[k] * ab_z + g_x_0_y_0_xxxxx_xzz[k];

                g_x_0_y_0_xxxxxz_yy[k] = -g_x_0_y_0_xxxxx_yy[k] * ab_z + g_x_0_y_0_xxxxx_yyz[k];

                g_x_0_y_0_xxxxxz_yz[k] = -g_x_0_y_0_xxxxx_yz[k] * ab_z + g_x_0_y_0_xxxxx_yzz[k];

                g_x_0_y_0_xxxxxz_zz[k] = -g_x_0_y_0_xxxxx_zz[k] * ab_z + g_x_0_y_0_xxxxx_zzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxyy_xx = cbuffer.data(id_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyy_xy = cbuffer.data(id_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyy_xz = cbuffer.data(id_geom_1010_off + 188 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyy_yy = cbuffer.data(id_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyy_yz = cbuffer.data(id_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyy_zz = cbuffer.data(id_geom_1010_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxy_xx, g_x_0_y_0_xxxxy_xxy, g_x_0_y_0_xxxxy_xy, g_x_0_y_0_xxxxy_xyy, g_x_0_y_0_xxxxy_xyz, g_x_0_y_0_xxxxy_xz, g_x_0_y_0_xxxxy_yy, g_x_0_y_0_xxxxy_yyy, g_x_0_y_0_xxxxy_yyz, g_x_0_y_0_xxxxy_yz, g_x_0_y_0_xxxxy_yzz, g_x_0_y_0_xxxxy_zz, g_x_0_y_0_xxxxyy_xx, g_x_0_y_0_xxxxyy_xy, g_x_0_y_0_xxxxyy_xz, g_x_0_y_0_xxxxyy_yy, g_x_0_y_0_xxxxyy_yz, g_x_0_y_0_xxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxyy_xx[k] = -g_x_0_y_0_xxxxy_xx[k] * ab_y + g_x_0_y_0_xxxxy_xxy[k];

                g_x_0_y_0_xxxxyy_xy[k] = -g_x_0_y_0_xxxxy_xy[k] * ab_y + g_x_0_y_0_xxxxy_xyy[k];

                g_x_0_y_0_xxxxyy_xz[k] = -g_x_0_y_0_xxxxy_xz[k] * ab_y + g_x_0_y_0_xxxxy_xyz[k];

                g_x_0_y_0_xxxxyy_yy[k] = -g_x_0_y_0_xxxxy_yy[k] * ab_y + g_x_0_y_0_xxxxy_yyy[k];

                g_x_0_y_0_xxxxyy_yz[k] = -g_x_0_y_0_xxxxy_yz[k] * ab_y + g_x_0_y_0_xxxxy_yyz[k];

                g_x_0_y_0_xxxxyy_zz[k] = -g_x_0_y_0_xxxxy_zz[k] * ab_y + g_x_0_y_0_xxxxy_yzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxyz_xx = cbuffer.data(id_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyz_xy = cbuffer.data(id_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyz_xz = cbuffer.data(id_geom_1010_off + 194 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyz_yy = cbuffer.data(id_geom_1010_off + 195 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyz_yz = cbuffer.data(id_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyz_zz = cbuffer.data(id_geom_1010_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxyz_xx, g_x_0_y_0_xxxxyz_xy, g_x_0_y_0_xxxxyz_xz, g_x_0_y_0_xxxxyz_yy, g_x_0_y_0_xxxxyz_yz, g_x_0_y_0_xxxxyz_zz, g_x_0_y_0_xxxxz_xx, g_x_0_y_0_xxxxz_xxy, g_x_0_y_0_xxxxz_xy, g_x_0_y_0_xxxxz_xyy, g_x_0_y_0_xxxxz_xyz, g_x_0_y_0_xxxxz_xz, g_x_0_y_0_xxxxz_yy, g_x_0_y_0_xxxxz_yyy, g_x_0_y_0_xxxxz_yyz, g_x_0_y_0_xxxxz_yz, g_x_0_y_0_xxxxz_yzz, g_x_0_y_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxyz_xx[k] = -g_x_0_y_0_xxxxz_xx[k] * ab_y + g_x_0_y_0_xxxxz_xxy[k];

                g_x_0_y_0_xxxxyz_xy[k] = -g_x_0_y_0_xxxxz_xy[k] * ab_y + g_x_0_y_0_xxxxz_xyy[k];

                g_x_0_y_0_xxxxyz_xz[k] = -g_x_0_y_0_xxxxz_xz[k] * ab_y + g_x_0_y_0_xxxxz_xyz[k];

                g_x_0_y_0_xxxxyz_yy[k] = -g_x_0_y_0_xxxxz_yy[k] * ab_y + g_x_0_y_0_xxxxz_yyy[k];

                g_x_0_y_0_xxxxyz_yz[k] = -g_x_0_y_0_xxxxz_yz[k] * ab_y + g_x_0_y_0_xxxxz_yyz[k];

                g_x_0_y_0_xxxxyz_zz[k] = -g_x_0_y_0_xxxxz_zz[k] * ab_y + g_x_0_y_0_xxxxz_yzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxzz_xx = cbuffer.data(id_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxzz_xy = cbuffer.data(id_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxzz_xz = cbuffer.data(id_geom_1010_off + 200 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxzz_yy = cbuffer.data(id_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxzz_yz = cbuffer.data(id_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxzz_zz = cbuffer.data(id_geom_1010_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxz_xx, g_x_0_y_0_xxxxz_xxz, g_x_0_y_0_xxxxz_xy, g_x_0_y_0_xxxxz_xyz, g_x_0_y_0_xxxxz_xz, g_x_0_y_0_xxxxz_xzz, g_x_0_y_0_xxxxz_yy, g_x_0_y_0_xxxxz_yyz, g_x_0_y_0_xxxxz_yz, g_x_0_y_0_xxxxz_yzz, g_x_0_y_0_xxxxz_zz, g_x_0_y_0_xxxxz_zzz, g_x_0_y_0_xxxxzz_xx, g_x_0_y_0_xxxxzz_xy, g_x_0_y_0_xxxxzz_xz, g_x_0_y_0_xxxxzz_yy, g_x_0_y_0_xxxxzz_yz, g_x_0_y_0_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxzz_xx[k] = -g_x_0_y_0_xxxxz_xx[k] * ab_z + g_x_0_y_0_xxxxz_xxz[k];

                g_x_0_y_0_xxxxzz_xy[k] = -g_x_0_y_0_xxxxz_xy[k] * ab_z + g_x_0_y_0_xxxxz_xyz[k];

                g_x_0_y_0_xxxxzz_xz[k] = -g_x_0_y_0_xxxxz_xz[k] * ab_z + g_x_0_y_0_xxxxz_xzz[k];

                g_x_0_y_0_xxxxzz_yy[k] = -g_x_0_y_0_xxxxz_yy[k] * ab_z + g_x_0_y_0_xxxxz_yyz[k];

                g_x_0_y_0_xxxxzz_yz[k] = -g_x_0_y_0_xxxxz_yz[k] * ab_z + g_x_0_y_0_xxxxz_yzz[k];

                g_x_0_y_0_xxxxzz_zz[k] = -g_x_0_y_0_xxxxz_zz[k] * ab_z + g_x_0_y_0_xxxxz_zzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxyyy_xx = cbuffer.data(id_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyy_xy = cbuffer.data(id_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyy_xz = cbuffer.data(id_geom_1010_off + 206 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyy_yy = cbuffer.data(id_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyy_yz = cbuffer.data(id_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyy_zz = cbuffer.data(id_geom_1010_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxyy_xx, g_x_0_y_0_xxxyy_xxy, g_x_0_y_0_xxxyy_xy, g_x_0_y_0_xxxyy_xyy, g_x_0_y_0_xxxyy_xyz, g_x_0_y_0_xxxyy_xz, g_x_0_y_0_xxxyy_yy, g_x_0_y_0_xxxyy_yyy, g_x_0_y_0_xxxyy_yyz, g_x_0_y_0_xxxyy_yz, g_x_0_y_0_xxxyy_yzz, g_x_0_y_0_xxxyy_zz, g_x_0_y_0_xxxyyy_xx, g_x_0_y_0_xxxyyy_xy, g_x_0_y_0_xxxyyy_xz, g_x_0_y_0_xxxyyy_yy, g_x_0_y_0_xxxyyy_yz, g_x_0_y_0_xxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxyyy_xx[k] = -g_x_0_y_0_xxxyy_xx[k] * ab_y + g_x_0_y_0_xxxyy_xxy[k];

                g_x_0_y_0_xxxyyy_xy[k] = -g_x_0_y_0_xxxyy_xy[k] * ab_y + g_x_0_y_0_xxxyy_xyy[k];

                g_x_0_y_0_xxxyyy_xz[k] = -g_x_0_y_0_xxxyy_xz[k] * ab_y + g_x_0_y_0_xxxyy_xyz[k];

                g_x_0_y_0_xxxyyy_yy[k] = -g_x_0_y_0_xxxyy_yy[k] * ab_y + g_x_0_y_0_xxxyy_yyy[k];

                g_x_0_y_0_xxxyyy_yz[k] = -g_x_0_y_0_xxxyy_yz[k] * ab_y + g_x_0_y_0_xxxyy_yyz[k];

                g_x_0_y_0_xxxyyy_zz[k] = -g_x_0_y_0_xxxyy_zz[k] * ab_y + g_x_0_y_0_xxxyy_yzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxyyz_xx = cbuffer.data(id_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyz_xy = cbuffer.data(id_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyz_xz = cbuffer.data(id_geom_1010_off + 212 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyz_yy = cbuffer.data(id_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyz_yz = cbuffer.data(id_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyz_zz = cbuffer.data(id_geom_1010_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxyyz_xx, g_x_0_y_0_xxxyyz_xy, g_x_0_y_0_xxxyyz_xz, g_x_0_y_0_xxxyyz_yy, g_x_0_y_0_xxxyyz_yz, g_x_0_y_0_xxxyyz_zz, g_x_0_y_0_xxxyz_xx, g_x_0_y_0_xxxyz_xxy, g_x_0_y_0_xxxyz_xy, g_x_0_y_0_xxxyz_xyy, g_x_0_y_0_xxxyz_xyz, g_x_0_y_0_xxxyz_xz, g_x_0_y_0_xxxyz_yy, g_x_0_y_0_xxxyz_yyy, g_x_0_y_0_xxxyz_yyz, g_x_0_y_0_xxxyz_yz, g_x_0_y_0_xxxyz_yzz, g_x_0_y_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxyyz_xx[k] = -g_x_0_y_0_xxxyz_xx[k] * ab_y + g_x_0_y_0_xxxyz_xxy[k];

                g_x_0_y_0_xxxyyz_xy[k] = -g_x_0_y_0_xxxyz_xy[k] * ab_y + g_x_0_y_0_xxxyz_xyy[k];

                g_x_0_y_0_xxxyyz_xz[k] = -g_x_0_y_0_xxxyz_xz[k] * ab_y + g_x_0_y_0_xxxyz_xyz[k];

                g_x_0_y_0_xxxyyz_yy[k] = -g_x_0_y_0_xxxyz_yy[k] * ab_y + g_x_0_y_0_xxxyz_yyy[k];

                g_x_0_y_0_xxxyyz_yz[k] = -g_x_0_y_0_xxxyz_yz[k] * ab_y + g_x_0_y_0_xxxyz_yyz[k];

                g_x_0_y_0_xxxyyz_zz[k] = -g_x_0_y_0_xxxyz_zz[k] * ab_y + g_x_0_y_0_xxxyz_yzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxyzz_xx = cbuffer.data(id_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyzz_xy = cbuffer.data(id_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyzz_xz = cbuffer.data(id_geom_1010_off + 218 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyzz_yy = cbuffer.data(id_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyzz_yz = cbuffer.data(id_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyzz_zz = cbuffer.data(id_geom_1010_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxyzz_xx, g_x_0_y_0_xxxyzz_xy, g_x_0_y_0_xxxyzz_xz, g_x_0_y_0_xxxyzz_yy, g_x_0_y_0_xxxyzz_yz, g_x_0_y_0_xxxyzz_zz, g_x_0_y_0_xxxzz_xx, g_x_0_y_0_xxxzz_xxy, g_x_0_y_0_xxxzz_xy, g_x_0_y_0_xxxzz_xyy, g_x_0_y_0_xxxzz_xyz, g_x_0_y_0_xxxzz_xz, g_x_0_y_0_xxxzz_yy, g_x_0_y_0_xxxzz_yyy, g_x_0_y_0_xxxzz_yyz, g_x_0_y_0_xxxzz_yz, g_x_0_y_0_xxxzz_yzz, g_x_0_y_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxyzz_xx[k] = -g_x_0_y_0_xxxzz_xx[k] * ab_y + g_x_0_y_0_xxxzz_xxy[k];

                g_x_0_y_0_xxxyzz_xy[k] = -g_x_0_y_0_xxxzz_xy[k] * ab_y + g_x_0_y_0_xxxzz_xyy[k];

                g_x_0_y_0_xxxyzz_xz[k] = -g_x_0_y_0_xxxzz_xz[k] * ab_y + g_x_0_y_0_xxxzz_xyz[k];

                g_x_0_y_0_xxxyzz_yy[k] = -g_x_0_y_0_xxxzz_yy[k] * ab_y + g_x_0_y_0_xxxzz_yyy[k];

                g_x_0_y_0_xxxyzz_yz[k] = -g_x_0_y_0_xxxzz_yz[k] * ab_y + g_x_0_y_0_xxxzz_yyz[k];

                g_x_0_y_0_xxxyzz_zz[k] = -g_x_0_y_0_xxxzz_zz[k] * ab_y + g_x_0_y_0_xxxzz_yzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxzzz_xx = cbuffer.data(id_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzzz_xy = cbuffer.data(id_geom_1010_off + 223 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzzz_xz = cbuffer.data(id_geom_1010_off + 224 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzzz_yy = cbuffer.data(id_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzzz_yz = cbuffer.data(id_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzzz_zz = cbuffer.data(id_geom_1010_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxzz_xx, g_x_0_y_0_xxxzz_xxz, g_x_0_y_0_xxxzz_xy, g_x_0_y_0_xxxzz_xyz, g_x_0_y_0_xxxzz_xz, g_x_0_y_0_xxxzz_xzz, g_x_0_y_0_xxxzz_yy, g_x_0_y_0_xxxzz_yyz, g_x_0_y_0_xxxzz_yz, g_x_0_y_0_xxxzz_yzz, g_x_0_y_0_xxxzz_zz, g_x_0_y_0_xxxzz_zzz, g_x_0_y_0_xxxzzz_xx, g_x_0_y_0_xxxzzz_xy, g_x_0_y_0_xxxzzz_xz, g_x_0_y_0_xxxzzz_yy, g_x_0_y_0_xxxzzz_yz, g_x_0_y_0_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxzzz_xx[k] = -g_x_0_y_0_xxxzz_xx[k] * ab_z + g_x_0_y_0_xxxzz_xxz[k];

                g_x_0_y_0_xxxzzz_xy[k] = -g_x_0_y_0_xxxzz_xy[k] * ab_z + g_x_0_y_0_xxxzz_xyz[k];

                g_x_0_y_0_xxxzzz_xz[k] = -g_x_0_y_0_xxxzz_xz[k] * ab_z + g_x_0_y_0_xxxzz_xzz[k];

                g_x_0_y_0_xxxzzz_yy[k] = -g_x_0_y_0_xxxzz_yy[k] * ab_z + g_x_0_y_0_xxxzz_yyz[k];

                g_x_0_y_0_xxxzzz_yz[k] = -g_x_0_y_0_xxxzz_yz[k] * ab_z + g_x_0_y_0_xxxzz_yzz[k];

                g_x_0_y_0_xxxzzz_zz[k] = -g_x_0_y_0_xxxzz_zz[k] * ab_z + g_x_0_y_0_xxxzz_zzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyyyy_xx = cbuffer.data(id_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyy_xy = cbuffer.data(id_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyy_xz = cbuffer.data(id_geom_1010_off + 230 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyy_yy = cbuffer.data(id_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyy_yz = cbuffer.data(id_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyy_zz = cbuffer.data(id_geom_1010_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyyy_xx, g_x_0_y_0_xxyyy_xxy, g_x_0_y_0_xxyyy_xy, g_x_0_y_0_xxyyy_xyy, g_x_0_y_0_xxyyy_xyz, g_x_0_y_0_xxyyy_xz, g_x_0_y_0_xxyyy_yy, g_x_0_y_0_xxyyy_yyy, g_x_0_y_0_xxyyy_yyz, g_x_0_y_0_xxyyy_yz, g_x_0_y_0_xxyyy_yzz, g_x_0_y_0_xxyyy_zz, g_x_0_y_0_xxyyyy_xx, g_x_0_y_0_xxyyyy_xy, g_x_0_y_0_xxyyyy_xz, g_x_0_y_0_xxyyyy_yy, g_x_0_y_0_xxyyyy_yz, g_x_0_y_0_xxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyyyy_xx[k] = -g_x_0_y_0_xxyyy_xx[k] * ab_y + g_x_0_y_0_xxyyy_xxy[k];

                g_x_0_y_0_xxyyyy_xy[k] = -g_x_0_y_0_xxyyy_xy[k] * ab_y + g_x_0_y_0_xxyyy_xyy[k];

                g_x_0_y_0_xxyyyy_xz[k] = -g_x_0_y_0_xxyyy_xz[k] * ab_y + g_x_0_y_0_xxyyy_xyz[k];

                g_x_0_y_0_xxyyyy_yy[k] = -g_x_0_y_0_xxyyy_yy[k] * ab_y + g_x_0_y_0_xxyyy_yyy[k];

                g_x_0_y_0_xxyyyy_yz[k] = -g_x_0_y_0_xxyyy_yz[k] * ab_y + g_x_0_y_0_xxyyy_yyz[k];

                g_x_0_y_0_xxyyyy_zz[k] = -g_x_0_y_0_xxyyy_zz[k] * ab_y + g_x_0_y_0_xxyyy_yzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyyyz_xx = cbuffer.data(id_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyz_xy = cbuffer.data(id_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyz_xz = cbuffer.data(id_geom_1010_off + 236 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyz_yy = cbuffer.data(id_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyz_yz = cbuffer.data(id_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyz_zz = cbuffer.data(id_geom_1010_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyyyz_xx, g_x_0_y_0_xxyyyz_xy, g_x_0_y_0_xxyyyz_xz, g_x_0_y_0_xxyyyz_yy, g_x_0_y_0_xxyyyz_yz, g_x_0_y_0_xxyyyz_zz, g_x_0_y_0_xxyyz_xx, g_x_0_y_0_xxyyz_xxy, g_x_0_y_0_xxyyz_xy, g_x_0_y_0_xxyyz_xyy, g_x_0_y_0_xxyyz_xyz, g_x_0_y_0_xxyyz_xz, g_x_0_y_0_xxyyz_yy, g_x_0_y_0_xxyyz_yyy, g_x_0_y_0_xxyyz_yyz, g_x_0_y_0_xxyyz_yz, g_x_0_y_0_xxyyz_yzz, g_x_0_y_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyyyz_xx[k] = -g_x_0_y_0_xxyyz_xx[k] * ab_y + g_x_0_y_0_xxyyz_xxy[k];

                g_x_0_y_0_xxyyyz_xy[k] = -g_x_0_y_0_xxyyz_xy[k] * ab_y + g_x_0_y_0_xxyyz_xyy[k];

                g_x_0_y_0_xxyyyz_xz[k] = -g_x_0_y_0_xxyyz_xz[k] * ab_y + g_x_0_y_0_xxyyz_xyz[k];

                g_x_0_y_0_xxyyyz_yy[k] = -g_x_0_y_0_xxyyz_yy[k] * ab_y + g_x_0_y_0_xxyyz_yyy[k];

                g_x_0_y_0_xxyyyz_yz[k] = -g_x_0_y_0_xxyyz_yz[k] * ab_y + g_x_0_y_0_xxyyz_yyz[k];

                g_x_0_y_0_xxyyyz_zz[k] = -g_x_0_y_0_xxyyz_zz[k] * ab_y + g_x_0_y_0_xxyyz_yzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyyzz_xx = cbuffer.data(id_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyzz_xy = cbuffer.data(id_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyzz_xz = cbuffer.data(id_geom_1010_off + 242 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyzz_yy = cbuffer.data(id_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyzz_yz = cbuffer.data(id_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyzz_zz = cbuffer.data(id_geom_1010_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyyzz_xx, g_x_0_y_0_xxyyzz_xy, g_x_0_y_0_xxyyzz_xz, g_x_0_y_0_xxyyzz_yy, g_x_0_y_0_xxyyzz_yz, g_x_0_y_0_xxyyzz_zz, g_x_0_y_0_xxyzz_xx, g_x_0_y_0_xxyzz_xxy, g_x_0_y_0_xxyzz_xy, g_x_0_y_0_xxyzz_xyy, g_x_0_y_0_xxyzz_xyz, g_x_0_y_0_xxyzz_xz, g_x_0_y_0_xxyzz_yy, g_x_0_y_0_xxyzz_yyy, g_x_0_y_0_xxyzz_yyz, g_x_0_y_0_xxyzz_yz, g_x_0_y_0_xxyzz_yzz, g_x_0_y_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyyzz_xx[k] = -g_x_0_y_0_xxyzz_xx[k] * ab_y + g_x_0_y_0_xxyzz_xxy[k];

                g_x_0_y_0_xxyyzz_xy[k] = -g_x_0_y_0_xxyzz_xy[k] * ab_y + g_x_0_y_0_xxyzz_xyy[k];

                g_x_0_y_0_xxyyzz_xz[k] = -g_x_0_y_0_xxyzz_xz[k] * ab_y + g_x_0_y_0_xxyzz_xyz[k];

                g_x_0_y_0_xxyyzz_yy[k] = -g_x_0_y_0_xxyzz_yy[k] * ab_y + g_x_0_y_0_xxyzz_yyy[k];

                g_x_0_y_0_xxyyzz_yz[k] = -g_x_0_y_0_xxyzz_yz[k] * ab_y + g_x_0_y_0_xxyzz_yyz[k];

                g_x_0_y_0_xxyyzz_zz[k] = -g_x_0_y_0_xxyzz_zz[k] * ab_y + g_x_0_y_0_xxyzz_yzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyzzz_xx = cbuffer.data(id_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzzz_xy = cbuffer.data(id_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzzz_xz = cbuffer.data(id_geom_1010_off + 248 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzzz_yy = cbuffer.data(id_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzzz_yz = cbuffer.data(id_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzzz_zz = cbuffer.data(id_geom_1010_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyzzz_xx, g_x_0_y_0_xxyzzz_xy, g_x_0_y_0_xxyzzz_xz, g_x_0_y_0_xxyzzz_yy, g_x_0_y_0_xxyzzz_yz, g_x_0_y_0_xxyzzz_zz, g_x_0_y_0_xxzzz_xx, g_x_0_y_0_xxzzz_xxy, g_x_0_y_0_xxzzz_xy, g_x_0_y_0_xxzzz_xyy, g_x_0_y_0_xxzzz_xyz, g_x_0_y_0_xxzzz_xz, g_x_0_y_0_xxzzz_yy, g_x_0_y_0_xxzzz_yyy, g_x_0_y_0_xxzzz_yyz, g_x_0_y_0_xxzzz_yz, g_x_0_y_0_xxzzz_yzz, g_x_0_y_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyzzz_xx[k] = -g_x_0_y_0_xxzzz_xx[k] * ab_y + g_x_0_y_0_xxzzz_xxy[k];

                g_x_0_y_0_xxyzzz_xy[k] = -g_x_0_y_0_xxzzz_xy[k] * ab_y + g_x_0_y_0_xxzzz_xyy[k];

                g_x_0_y_0_xxyzzz_xz[k] = -g_x_0_y_0_xxzzz_xz[k] * ab_y + g_x_0_y_0_xxzzz_xyz[k];

                g_x_0_y_0_xxyzzz_yy[k] = -g_x_0_y_0_xxzzz_yy[k] * ab_y + g_x_0_y_0_xxzzz_yyy[k];

                g_x_0_y_0_xxyzzz_yz[k] = -g_x_0_y_0_xxzzz_yz[k] * ab_y + g_x_0_y_0_xxzzz_yyz[k];

                g_x_0_y_0_xxyzzz_zz[k] = -g_x_0_y_0_xxzzz_zz[k] * ab_y + g_x_0_y_0_xxzzz_yzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxzzzz_xx = cbuffer.data(id_geom_1010_off + 252 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzzz_xy = cbuffer.data(id_geom_1010_off + 253 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzzz_xz = cbuffer.data(id_geom_1010_off + 254 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzzz_yy = cbuffer.data(id_geom_1010_off + 255 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzzz_yz = cbuffer.data(id_geom_1010_off + 256 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzzz_zz = cbuffer.data(id_geom_1010_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxzzz_xx, g_x_0_y_0_xxzzz_xxz, g_x_0_y_0_xxzzz_xy, g_x_0_y_0_xxzzz_xyz, g_x_0_y_0_xxzzz_xz, g_x_0_y_0_xxzzz_xzz, g_x_0_y_0_xxzzz_yy, g_x_0_y_0_xxzzz_yyz, g_x_0_y_0_xxzzz_yz, g_x_0_y_0_xxzzz_yzz, g_x_0_y_0_xxzzz_zz, g_x_0_y_0_xxzzz_zzz, g_x_0_y_0_xxzzzz_xx, g_x_0_y_0_xxzzzz_xy, g_x_0_y_0_xxzzzz_xz, g_x_0_y_0_xxzzzz_yy, g_x_0_y_0_xxzzzz_yz, g_x_0_y_0_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxzzzz_xx[k] = -g_x_0_y_0_xxzzz_xx[k] * ab_z + g_x_0_y_0_xxzzz_xxz[k];

                g_x_0_y_0_xxzzzz_xy[k] = -g_x_0_y_0_xxzzz_xy[k] * ab_z + g_x_0_y_0_xxzzz_xyz[k];

                g_x_0_y_0_xxzzzz_xz[k] = -g_x_0_y_0_xxzzz_xz[k] * ab_z + g_x_0_y_0_xxzzz_xzz[k];

                g_x_0_y_0_xxzzzz_yy[k] = -g_x_0_y_0_xxzzz_yy[k] * ab_z + g_x_0_y_0_xxzzz_yyz[k];

                g_x_0_y_0_xxzzzz_yz[k] = -g_x_0_y_0_xxzzz_yz[k] * ab_z + g_x_0_y_0_xxzzz_yzz[k];

                g_x_0_y_0_xxzzzz_zz[k] = -g_x_0_y_0_xxzzz_zz[k] * ab_z + g_x_0_y_0_xxzzz_zzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyyyy_xx = cbuffer.data(id_geom_1010_off + 258 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyy_xy = cbuffer.data(id_geom_1010_off + 259 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyy_xz = cbuffer.data(id_geom_1010_off + 260 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyy_yy = cbuffer.data(id_geom_1010_off + 261 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyy_yz = cbuffer.data(id_geom_1010_off + 262 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyy_zz = cbuffer.data(id_geom_1010_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyyy_xx, g_x_0_y_0_xyyyy_xxy, g_x_0_y_0_xyyyy_xy, g_x_0_y_0_xyyyy_xyy, g_x_0_y_0_xyyyy_xyz, g_x_0_y_0_xyyyy_xz, g_x_0_y_0_xyyyy_yy, g_x_0_y_0_xyyyy_yyy, g_x_0_y_0_xyyyy_yyz, g_x_0_y_0_xyyyy_yz, g_x_0_y_0_xyyyy_yzz, g_x_0_y_0_xyyyy_zz, g_x_0_y_0_xyyyyy_xx, g_x_0_y_0_xyyyyy_xy, g_x_0_y_0_xyyyyy_xz, g_x_0_y_0_xyyyyy_yy, g_x_0_y_0_xyyyyy_yz, g_x_0_y_0_xyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyyyy_xx[k] = -g_x_0_y_0_xyyyy_xx[k] * ab_y + g_x_0_y_0_xyyyy_xxy[k];

                g_x_0_y_0_xyyyyy_xy[k] = -g_x_0_y_0_xyyyy_xy[k] * ab_y + g_x_0_y_0_xyyyy_xyy[k];

                g_x_0_y_0_xyyyyy_xz[k] = -g_x_0_y_0_xyyyy_xz[k] * ab_y + g_x_0_y_0_xyyyy_xyz[k];

                g_x_0_y_0_xyyyyy_yy[k] = -g_x_0_y_0_xyyyy_yy[k] * ab_y + g_x_0_y_0_xyyyy_yyy[k];

                g_x_0_y_0_xyyyyy_yz[k] = -g_x_0_y_0_xyyyy_yz[k] * ab_y + g_x_0_y_0_xyyyy_yyz[k];

                g_x_0_y_0_xyyyyy_zz[k] = -g_x_0_y_0_xyyyy_zz[k] * ab_y + g_x_0_y_0_xyyyy_yzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyyyz_xx = cbuffer.data(id_geom_1010_off + 264 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyz_xy = cbuffer.data(id_geom_1010_off + 265 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyz_xz = cbuffer.data(id_geom_1010_off + 266 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyz_yy = cbuffer.data(id_geom_1010_off + 267 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyz_yz = cbuffer.data(id_geom_1010_off + 268 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyz_zz = cbuffer.data(id_geom_1010_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyyyz_xx, g_x_0_y_0_xyyyyz_xy, g_x_0_y_0_xyyyyz_xz, g_x_0_y_0_xyyyyz_yy, g_x_0_y_0_xyyyyz_yz, g_x_0_y_0_xyyyyz_zz, g_x_0_y_0_xyyyz_xx, g_x_0_y_0_xyyyz_xxy, g_x_0_y_0_xyyyz_xy, g_x_0_y_0_xyyyz_xyy, g_x_0_y_0_xyyyz_xyz, g_x_0_y_0_xyyyz_xz, g_x_0_y_0_xyyyz_yy, g_x_0_y_0_xyyyz_yyy, g_x_0_y_0_xyyyz_yyz, g_x_0_y_0_xyyyz_yz, g_x_0_y_0_xyyyz_yzz, g_x_0_y_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyyyz_xx[k] = -g_x_0_y_0_xyyyz_xx[k] * ab_y + g_x_0_y_0_xyyyz_xxy[k];

                g_x_0_y_0_xyyyyz_xy[k] = -g_x_0_y_0_xyyyz_xy[k] * ab_y + g_x_0_y_0_xyyyz_xyy[k];

                g_x_0_y_0_xyyyyz_xz[k] = -g_x_0_y_0_xyyyz_xz[k] * ab_y + g_x_0_y_0_xyyyz_xyz[k];

                g_x_0_y_0_xyyyyz_yy[k] = -g_x_0_y_0_xyyyz_yy[k] * ab_y + g_x_0_y_0_xyyyz_yyy[k];

                g_x_0_y_0_xyyyyz_yz[k] = -g_x_0_y_0_xyyyz_yz[k] * ab_y + g_x_0_y_0_xyyyz_yyz[k];

                g_x_0_y_0_xyyyyz_zz[k] = -g_x_0_y_0_xyyyz_zz[k] * ab_y + g_x_0_y_0_xyyyz_yzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyyzz_xx = cbuffer.data(id_geom_1010_off + 270 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyzz_xy = cbuffer.data(id_geom_1010_off + 271 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyzz_xz = cbuffer.data(id_geom_1010_off + 272 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyzz_yy = cbuffer.data(id_geom_1010_off + 273 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyzz_yz = cbuffer.data(id_geom_1010_off + 274 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyzz_zz = cbuffer.data(id_geom_1010_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyyzz_xx, g_x_0_y_0_xyyyzz_xy, g_x_0_y_0_xyyyzz_xz, g_x_0_y_0_xyyyzz_yy, g_x_0_y_0_xyyyzz_yz, g_x_0_y_0_xyyyzz_zz, g_x_0_y_0_xyyzz_xx, g_x_0_y_0_xyyzz_xxy, g_x_0_y_0_xyyzz_xy, g_x_0_y_0_xyyzz_xyy, g_x_0_y_0_xyyzz_xyz, g_x_0_y_0_xyyzz_xz, g_x_0_y_0_xyyzz_yy, g_x_0_y_0_xyyzz_yyy, g_x_0_y_0_xyyzz_yyz, g_x_0_y_0_xyyzz_yz, g_x_0_y_0_xyyzz_yzz, g_x_0_y_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyyzz_xx[k] = -g_x_0_y_0_xyyzz_xx[k] * ab_y + g_x_0_y_0_xyyzz_xxy[k];

                g_x_0_y_0_xyyyzz_xy[k] = -g_x_0_y_0_xyyzz_xy[k] * ab_y + g_x_0_y_0_xyyzz_xyy[k];

                g_x_0_y_0_xyyyzz_xz[k] = -g_x_0_y_0_xyyzz_xz[k] * ab_y + g_x_0_y_0_xyyzz_xyz[k];

                g_x_0_y_0_xyyyzz_yy[k] = -g_x_0_y_0_xyyzz_yy[k] * ab_y + g_x_0_y_0_xyyzz_yyy[k];

                g_x_0_y_0_xyyyzz_yz[k] = -g_x_0_y_0_xyyzz_yz[k] * ab_y + g_x_0_y_0_xyyzz_yyz[k];

                g_x_0_y_0_xyyyzz_zz[k] = -g_x_0_y_0_xyyzz_zz[k] * ab_y + g_x_0_y_0_xyyzz_yzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyzzz_xx = cbuffer.data(id_geom_1010_off + 276 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzzz_xy = cbuffer.data(id_geom_1010_off + 277 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzzz_xz = cbuffer.data(id_geom_1010_off + 278 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzzz_yy = cbuffer.data(id_geom_1010_off + 279 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzzz_yz = cbuffer.data(id_geom_1010_off + 280 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzzz_zz = cbuffer.data(id_geom_1010_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyzzz_xx, g_x_0_y_0_xyyzzz_xy, g_x_0_y_0_xyyzzz_xz, g_x_0_y_0_xyyzzz_yy, g_x_0_y_0_xyyzzz_yz, g_x_0_y_0_xyyzzz_zz, g_x_0_y_0_xyzzz_xx, g_x_0_y_0_xyzzz_xxy, g_x_0_y_0_xyzzz_xy, g_x_0_y_0_xyzzz_xyy, g_x_0_y_0_xyzzz_xyz, g_x_0_y_0_xyzzz_xz, g_x_0_y_0_xyzzz_yy, g_x_0_y_0_xyzzz_yyy, g_x_0_y_0_xyzzz_yyz, g_x_0_y_0_xyzzz_yz, g_x_0_y_0_xyzzz_yzz, g_x_0_y_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyzzz_xx[k] = -g_x_0_y_0_xyzzz_xx[k] * ab_y + g_x_0_y_0_xyzzz_xxy[k];

                g_x_0_y_0_xyyzzz_xy[k] = -g_x_0_y_0_xyzzz_xy[k] * ab_y + g_x_0_y_0_xyzzz_xyy[k];

                g_x_0_y_0_xyyzzz_xz[k] = -g_x_0_y_0_xyzzz_xz[k] * ab_y + g_x_0_y_0_xyzzz_xyz[k];

                g_x_0_y_0_xyyzzz_yy[k] = -g_x_0_y_0_xyzzz_yy[k] * ab_y + g_x_0_y_0_xyzzz_yyy[k];

                g_x_0_y_0_xyyzzz_yz[k] = -g_x_0_y_0_xyzzz_yz[k] * ab_y + g_x_0_y_0_xyzzz_yyz[k];

                g_x_0_y_0_xyyzzz_zz[k] = -g_x_0_y_0_xyzzz_zz[k] * ab_y + g_x_0_y_0_xyzzz_yzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyzzzz_xx = cbuffer.data(id_geom_1010_off + 282 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzzz_xy = cbuffer.data(id_geom_1010_off + 283 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzzz_xz = cbuffer.data(id_geom_1010_off + 284 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzzz_yy = cbuffer.data(id_geom_1010_off + 285 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzzz_yz = cbuffer.data(id_geom_1010_off + 286 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzzz_zz = cbuffer.data(id_geom_1010_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyzzzz_xx, g_x_0_y_0_xyzzzz_xy, g_x_0_y_0_xyzzzz_xz, g_x_0_y_0_xyzzzz_yy, g_x_0_y_0_xyzzzz_yz, g_x_0_y_0_xyzzzz_zz, g_x_0_y_0_xzzzz_xx, g_x_0_y_0_xzzzz_xxy, g_x_0_y_0_xzzzz_xy, g_x_0_y_0_xzzzz_xyy, g_x_0_y_0_xzzzz_xyz, g_x_0_y_0_xzzzz_xz, g_x_0_y_0_xzzzz_yy, g_x_0_y_0_xzzzz_yyy, g_x_0_y_0_xzzzz_yyz, g_x_0_y_0_xzzzz_yz, g_x_0_y_0_xzzzz_yzz, g_x_0_y_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyzzzz_xx[k] = -g_x_0_y_0_xzzzz_xx[k] * ab_y + g_x_0_y_0_xzzzz_xxy[k];

                g_x_0_y_0_xyzzzz_xy[k] = -g_x_0_y_0_xzzzz_xy[k] * ab_y + g_x_0_y_0_xzzzz_xyy[k];

                g_x_0_y_0_xyzzzz_xz[k] = -g_x_0_y_0_xzzzz_xz[k] * ab_y + g_x_0_y_0_xzzzz_xyz[k];

                g_x_0_y_0_xyzzzz_yy[k] = -g_x_0_y_0_xzzzz_yy[k] * ab_y + g_x_0_y_0_xzzzz_yyy[k];

                g_x_0_y_0_xyzzzz_yz[k] = -g_x_0_y_0_xzzzz_yz[k] * ab_y + g_x_0_y_0_xzzzz_yyz[k];

                g_x_0_y_0_xyzzzz_zz[k] = -g_x_0_y_0_xzzzz_zz[k] * ab_y + g_x_0_y_0_xzzzz_yzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xzzzzz_xx = cbuffer.data(id_geom_1010_off + 288 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzzz_xy = cbuffer.data(id_geom_1010_off + 289 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzzz_xz = cbuffer.data(id_geom_1010_off + 290 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzzz_yy = cbuffer.data(id_geom_1010_off + 291 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzzz_yz = cbuffer.data(id_geom_1010_off + 292 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzzz_zz = cbuffer.data(id_geom_1010_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xzzzz_xx, g_x_0_y_0_xzzzz_xxz, g_x_0_y_0_xzzzz_xy, g_x_0_y_0_xzzzz_xyz, g_x_0_y_0_xzzzz_xz, g_x_0_y_0_xzzzz_xzz, g_x_0_y_0_xzzzz_yy, g_x_0_y_0_xzzzz_yyz, g_x_0_y_0_xzzzz_yz, g_x_0_y_0_xzzzz_yzz, g_x_0_y_0_xzzzz_zz, g_x_0_y_0_xzzzz_zzz, g_x_0_y_0_xzzzzz_xx, g_x_0_y_0_xzzzzz_xy, g_x_0_y_0_xzzzzz_xz, g_x_0_y_0_xzzzzz_yy, g_x_0_y_0_xzzzzz_yz, g_x_0_y_0_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xzzzzz_xx[k] = -g_x_0_y_0_xzzzz_xx[k] * ab_z + g_x_0_y_0_xzzzz_xxz[k];

                g_x_0_y_0_xzzzzz_xy[k] = -g_x_0_y_0_xzzzz_xy[k] * ab_z + g_x_0_y_0_xzzzz_xyz[k];

                g_x_0_y_0_xzzzzz_xz[k] = -g_x_0_y_0_xzzzz_xz[k] * ab_z + g_x_0_y_0_xzzzz_xzz[k];

                g_x_0_y_0_xzzzzz_yy[k] = -g_x_0_y_0_xzzzz_yy[k] * ab_z + g_x_0_y_0_xzzzz_yyz[k];

                g_x_0_y_0_xzzzzz_yz[k] = -g_x_0_y_0_xzzzz_yz[k] * ab_z + g_x_0_y_0_xzzzz_yzz[k];

                g_x_0_y_0_xzzzzz_zz[k] = -g_x_0_y_0_xzzzz_zz[k] * ab_z + g_x_0_y_0_xzzzz_zzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyyyy_xx = cbuffer.data(id_geom_1010_off + 294 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyy_xy = cbuffer.data(id_geom_1010_off + 295 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyy_xz = cbuffer.data(id_geom_1010_off + 296 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyy_yy = cbuffer.data(id_geom_1010_off + 297 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyy_yz = cbuffer.data(id_geom_1010_off + 298 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyy_zz = cbuffer.data(id_geom_1010_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyyy_xx, g_x_0_y_0_yyyyy_xxy, g_x_0_y_0_yyyyy_xy, g_x_0_y_0_yyyyy_xyy, g_x_0_y_0_yyyyy_xyz, g_x_0_y_0_yyyyy_xz, g_x_0_y_0_yyyyy_yy, g_x_0_y_0_yyyyy_yyy, g_x_0_y_0_yyyyy_yyz, g_x_0_y_0_yyyyy_yz, g_x_0_y_0_yyyyy_yzz, g_x_0_y_0_yyyyy_zz, g_x_0_y_0_yyyyyy_xx, g_x_0_y_0_yyyyyy_xy, g_x_0_y_0_yyyyyy_xz, g_x_0_y_0_yyyyyy_yy, g_x_0_y_0_yyyyyy_yz, g_x_0_y_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyyyy_xx[k] = -g_x_0_y_0_yyyyy_xx[k] * ab_y + g_x_0_y_0_yyyyy_xxy[k];

                g_x_0_y_0_yyyyyy_xy[k] = -g_x_0_y_0_yyyyy_xy[k] * ab_y + g_x_0_y_0_yyyyy_xyy[k];

                g_x_0_y_0_yyyyyy_xz[k] = -g_x_0_y_0_yyyyy_xz[k] * ab_y + g_x_0_y_0_yyyyy_xyz[k];

                g_x_0_y_0_yyyyyy_yy[k] = -g_x_0_y_0_yyyyy_yy[k] * ab_y + g_x_0_y_0_yyyyy_yyy[k];

                g_x_0_y_0_yyyyyy_yz[k] = -g_x_0_y_0_yyyyy_yz[k] * ab_y + g_x_0_y_0_yyyyy_yyz[k];

                g_x_0_y_0_yyyyyy_zz[k] = -g_x_0_y_0_yyyyy_zz[k] * ab_y + g_x_0_y_0_yyyyy_yzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyyyz_xx = cbuffer.data(id_geom_1010_off + 300 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyz_xy = cbuffer.data(id_geom_1010_off + 301 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyz_xz = cbuffer.data(id_geom_1010_off + 302 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyz_yy = cbuffer.data(id_geom_1010_off + 303 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyz_yz = cbuffer.data(id_geom_1010_off + 304 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyz_zz = cbuffer.data(id_geom_1010_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyyyz_xx, g_x_0_y_0_yyyyyz_xy, g_x_0_y_0_yyyyyz_xz, g_x_0_y_0_yyyyyz_yy, g_x_0_y_0_yyyyyz_yz, g_x_0_y_0_yyyyyz_zz, g_x_0_y_0_yyyyz_xx, g_x_0_y_0_yyyyz_xxy, g_x_0_y_0_yyyyz_xy, g_x_0_y_0_yyyyz_xyy, g_x_0_y_0_yyyyz_xyz, g_x_0_y_0_yyyyz_xz, g_x_0_y_0_yyyyz_yy, g_x_0_y_0_yyyyz_yyy, g_x_0_y_0_yyyyz_yyz, g_x_0_y_0_yyyyz_yz, g_x_0_y_0_yyyyz_yzz, g_x_0_y_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyyyz_xx[k] = -g_x_0_y_0_yyyyz_xx[k] * ab_y + g_x_0_y_0_yyyyz_xxy[k];

                g_x_0_y_0_yyyyyz_xy[k] = -g_x_0_y_0_yyyyz_xy[k] * ab_y + g_x_0_y_0_yyyyz_xyy[k];

                g_x_0_y_0_yyyyyz_xz[k] = -g_x_0_y_0_yyyyz_xz[k] * ab_y + g_x_0_y_0_yyyyz_xyz[k];

                g_x_0_y_0_yyyyyz_yy[k] = -g_x_0_y_0_yyyyz_yy[k] * ab_y + g_x_0_y_0_yyyyz_yyy[k];

                g_x_0_y_0_yyyyyz_yz[k] = -g_x_0_y_0_yyyyz_yz[k] * ab_y + g_x_0_y_0_yyyyz_yyz[k];

                g_x_0_y_0_yyyyyz_zz[k] = -g_x_0_y_0_yyyyz_zz[k] * ab_y + g_x_0_y_0_yyyyz_yzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyyzz_xx = cbuffer.data(id_geom_1010_off + 306 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyzz_xy = cbuffer.data(id_geom_1010_off + 307 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyzz_xz = cbuffer.data(id_geom_1010_off + 308 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyzz_yy = cbuffer.data(id_geom_1010_off + 309 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyzz_yz = cbuffer.data(id_geom_1010_off + 310 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyzz_zz = cbuffer.data(id_geom_1010_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyyzz_xx, g_x_0_y_0_yyyyzz_xy, g_x_0_y_0_yyyyzz_xz, g_x_0_y_0_yyyyzz_yy, g_x_0_y_0_yyyyzz_yz, g_x_0_y_0_yyyyzz_zz, g_x_0_y_0_yyyzz_xx, g_x_0_y_0_yyyzz_xxy, g_x_0_y_0_yyyzz_xy, g_x_0_y_0_yyyzz_xyy, g_x_0_y_0_yyyzz_xyz, g_x_0_y_0_yyyzz_xz, g_x_0_y_0_yyyzz_yy, g_x_0_y_0_yyyzz_yyy, g_x_0_y_0_yyyzz_yyz, g_x_0_y_0_yyyzz_yz, g_x_0_y_0_yyyzz_yzz, g_x_0_y_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyyzz_xx[k] = -g_x_0_y_0_yyyzz_xx[k] * ab_y + g_x_0_y_0_yyyzz_xxy[k];

                g_x_0_y_0_yyyyzz_xy[k] = -g_x_0_y_0_yyyzz_xy[k] * ab_y + g_x_0_y_0_yyyzz_xyy[k];

                g_x_0_y_0_yyyyzz_xz[k] = -g_x_0_y_0_yyyzz_xz[k] * ab_y + g_x_0_y_0_yyyzz_xyz[k];

                g_x_0_y_0_yyyyzz_yy[k] = -g_x_0_y_0_yyyzz_yy[k] * ab_y + g_x_0_y_0_yyyzz_yyy[k];

                g_x_0_y_0_yyyyzz_yz[k] = -g_x_0_y_0_yyyzz_yz[k] * ab_y + g_x_0_y_0_yyyzz_yyz[k];

                g_x_0_y_0_yyyyzz_zz[k] = -g_x_0_y_0_yyyzz_zz[k] * ab_y + g_x_0_y_0_yyyzz_yzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyzzz_xx = cbuffer.data(id_geom_1010_off + 312 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzzz_xy = cbuffer.data(id_geom_1010_off + 313 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzzz_xz = cbuffer.data(id_geom_1010_off + 314 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzzz_yy = cbuffer.data(id_geom_1010_off + 315 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzzz_yz = cbuffer.data(id_geom_1010_off + 316 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzzz_zz = cbuffer.data(id_geom_1010_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyzzz_xx, g_x_0_y_0_yyyzzz_xy, g_x_0_y_0_yyyzzz_xz, g_x_0_y_0_yyyzzz_yy, g_x_0_y_0_yyyzzz_yz, g_x_0_y_0_yyyzzz_zz, g_x_0_y_0_yyzzz_xx, g_x_0_y_0_yyzzz_xxy, g_x_0_y_0_yyzzz_xy, g_x_0_y_0_yyzzz_xyy, g_x_0_y_0_yyzzz_xyz, g_x_0_y_0_yyzzz_xz, g_x_0_y_0_yyzzz_yy, g_x_0_y_0_yyzzz_yyy, g_x_0_y_0_yyzzz_yyz, g_x_0_y_0_yyzzz_yz, g_x_0_y_0_yyzzz_yzz, g_x_0_y_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyzzz_xx[k] = -g_x_0_y_0_yyzzz_xx[k] * ab_y + g_x_0_y_0_yyzzz_xxy[k];

                g_x_0_y_0_yyyzzz_xy[k] = -g_x_0_y_0_yyzzz_xy[k] * ab_y + g_x_0_y_0_yyzzz_xyy[k];

                g_x_0_y_0_yyyzzz_xz[k] = -g_x_0_y_0_yyzzz_xz[k] * ab_y + g_x_0_y_0_yyzzz_xyz[k];

                g_x_0_y_0_yyyzzz_yy[k] = -g_x_0_y_0_yyzzz_yy[k] * ab_y + g_x_0_y_0_yyzzz_yyy[k];

                g_x_0_y_0_yyyzzz_yz[k] = -g_x_0_y_0_yyzzz_yz[k] * ab_y + g_x_0_y_0_yyzzz_yyz[k];

                g_x_0_y_0_yyyzzz_zz[k] = -g_x_0_y_0_yyzzz_zz[k] * ab_y + g_x_0_y_0_yyzzz_yzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyzzzz_xx = cbuffer.data(id_geom_1010_off + 318 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzzz_xy = cbuffer.data(id_geom_1010_off + 319 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzzz_xz = cbuffer.data(id_geom_1010_off + 320 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzzz_yy = cbuffer.data(id_geom_1010_off + 321 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzzz_yz = cbuffer.data(id_geom_1010_off + 322 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzzz_zz = cbuffer.data(id_geom_1010_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyzzzz_xx, g_x_0_y_0_yyzzzz_xy, g_x_0_y_0_yyzzzz_xz, g_x_0_y_0_yyzzzz_yy, g_x_0_y_0_yyzzzz_yz, g_x_0_y_0_yyzzzz_zz, g_x_0_y_0_yzzzz_xx, g_x_0_y_0_yzzzz_xxy, g_x_0_y_0_yzzzz_xy, g_x_0_y_0_yzzzz_xyy, g_x_0_y_0_yzzzz_xyz, g_x_0_y_0_yzzzz_xz, g_x_0_y_0_yzzzz_yy, g_x_0_y_0_yzzzz_yyy, g_x_0_y_0_yzzzz_yyz, g_x_0_y_0_yzzzz_yz, g_x_0_y_0_yzzzz_yzz, g_x_0_y_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyzzzz_xx[k] = -g_x_0_y_0_yzzzz_xx[k] * ab_y + g_x_0_y_0_yzzzz_xxy[k];

                g_x_0_y_0_yyzzzz_xy[k] = -g_x_0_y_0_yzzzz_xy[k] * ab_y + g_x_0_y_0_yzzzz_xyy[k];

                g_x_0_y_0_yyzzzz_xz[k] = -g_x_0_y_0_yzzzz_xz[k] * ab_y + g_x_0_y_0_yzzzz_xyz[k];

                g_x_0_y_0_yyzzzz_yy[k] = -g_x_0_y_0_yzzzz_yy[k] * ab_y + g_x_0_y_0_yzzzz_yyy[k];

                g_x_0_y_0_yyzzzz_yz[k] = -g_x_0_y_0_yzzzz_yz[k] * ab_y + g_x_0_y_0_yzzzz_yyz[k];

                g_x_0_y_0_yyzzzz_zz[k] = -g_x_0_y_0_yzzzz_zz[k] * ab_y + g_x_0_y_0_yzzzz_yzz[k];
            }

            /// Set up 324-330 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yzzzzz_xx = cbuffer.data(id_geom_1010_off + 324 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzzz_xy = cbuffer.data(id_geom_1010_off + 325 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzzz_xz = cbuffer.data(id_geom_1010_off + 326 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzzz_yy = cbuffer.data(id_geom_1010_off + 327 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzzz_yz = cbuffer.data(id_geom_1010_off + 328 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzzz_zz = cbuffer.data(id_geom_1010_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yzzzzz_xx, g_x_0_y_0_yzzzzz_xy, g_x_0_y_0_yzzzzz_xz, g_x_0_y_0_yzzzzz_yy, g_x_0_y_0_yzzzzz_yz, g_x_0_y_0_yzzzzz_zz, g_x_0_y_0_zzzzz_xx, g_x_0_y_0_zzzzz_xxy, g_x_0_y_0_zzzzz_xy, g_x_0_y_0_zzzzz_xyy, g_x_0_y_0_zzzzz_xyz, g_x_0_y_0_zzzzz_xz, g_x_0_y_0_zzzzz_yy, g_x_0_y_0_zzzzz_yyy, g_x_0_y_0_zzzzz_yyz, g_x_0_y_0_zzzzz_yz, g_x_0_y_0_zzzzz_yzz, g_x_0_y_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yzzzzz_xx[k] = -g_x_0_y_0_zzzzz_xx[k] * ab_y + g_x_0_y_0_zzzzz_xxy[k];

                g_x_0_y_0_yzzzzz_xy[k] = -g_x_0_y_0_zzzzz_xy[k] * ab_y + g_x_0_y_0_zzzzz_xyy[k];

                g_x_0_y_0_yzzzzz_xz[k] = -g_x_0_y_0_zzzzz_xz[k] * ab_y + g_x_0_y_0_zzzzz_xyz[k];

                g_x_0_y_0_yzzzzz_yy[k] = -g_x_0_y_0_zzzzz_yy[k] * ab_y + g_x_0_y_0_zzzzz_yyy[k];

                g_x_0_y_0_yzzzzz_yz[k] = -g_x_0_y_0_zzzzz_yz[k] * ab_y + g_x_0_y_0_zzzzz_yyz[k];

                g_x_0_y_0_yzzzzz_zz[k] = -g_x_0_y_0_zzzzz_zz[k] * ab_y + g_x_0_y_0_zzzzz_yzz[k];
            }

            /// Set up 330-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_zzzzzz_xx = cbuffer.data(id_geom_1010_off + 330 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzzz_xy = cbuffer.data(id_geom_1010_off + 331 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzzz_xz = cbuffer.data(id_geom_1010_off + 332 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzzz_yy = cbuffer.data(id_geom_1010_off + 333 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzzz_yz = cbuffer.data(id_geom_1010_off + 334 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzzz_zz = cbuffer.data(id_geom_1010_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_zzzzz_xx, g_x_0_y_0_zzzzz_xxz, g_x_0_y_0_zzzzz_xy, g_x_0_y_0_zzzzz_xyz, g_x_0_y_0_zzzzz_xz, g_x_0_y_0_zzzzz_xzz, g_x_0_y_0_zzzzz_yy, g_x_0_y_0_zzzzz_yyz, g_x_0_y_0_zzzzz_yz, g_x_0_y_0_zzzzz_yzz, g_x_0_y_0_zzzzz_zz, g_x_0_y_0_zzzzz_zzz, g_x_0_y_0_zzzzzz_xx, g_x_0_y_0_zzzzzz_xy, g_x_0_y_0_zzzzzz_xz, g_x_0_y_0_zzzzzz_yy, g_x_0_y_0_zzzzzz_yz, g_x_0_y_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_zzzzzz_xx[k] = -g_x_0_y_0_zzzzz_xx[k] * ab_z + g_x_0_y_0_zzzzz_xxz[k];

                g_x_0_y_0_zzzzzz_xy[k] = -g_x_0_y_0_zzzzz_xy[k] * ab_z + g_x_0_y_0_zzzzz_xyz[k];

                g_x_0_y_0_zzzzzz_xz[k] = -g_x_0_y_0_zzzzz_xz[k] * ab_z + g_x_0_y_0_zzzzz_xzz[k];

                g_x_0_y_0_zzzzzz_yy[k] = -g_x_0_y_0_zzzzz_yy[k] * ab_z + g_x_0_y_0_zzzzz_yyz[k];

                g_x_0_y_0_zzzzzz_yz[k] = -g_x_0_y_0_zzzzz_yz[k] * ab_z + g_x_0_y_0_zzzzz_yzz[k];

                g_x_0_y_0_zzzzzz_zz[k] = -g_x_0_y_0_zzzzz_zz[k] * ab_z + g_x_0_y_0_zzzzz_zzz[k];
            }

            /// Set up 336-342 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxxx_xx = cbuffer.data(id_geom_1010_off + 336 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxx_xy = cbuffer.data(id_geom_1010_off + 337 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxx_xz = cbuffer.data(id_geom_1010_off + 338 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxx_yy = cbuffer.data(id_geom_1010_off + 339 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxx_yz = cbuffer.data(id_geom_1010_off + 340 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxx_zz = cbuffer.data(id_geom_1010_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_xxxxx_xx, g_0_0_z_0_xxxxx_xy, g_0_0_z_0_xxxxx_xz, g_0_0_z_0_xxxxx_yy, g_0_0_z_0_xxxxx_yz, g_0_0_z_0_xxxxx_zz, g_x_0_z_0_xxxxx_xx, g_x_0_z_0_xxxxx_xxx, g_x_0_z_0_xxxxx_xxy, g_x_0_z_0_xxxxx_xxz, g_x_0_z_0_xxxxx_xy, g_x_0_z_0_xxxxx_xyy, g_x_0_z_0_xxxxx_xyz, g_x_0_z_0_xxxxx_xz, g_x_0_z_0_xxxxx_xzz, g_x_0_z_0_xxxxx_yy, g_x_0_z_0_xxxxx_yz, g_x_0_z_0_xxxxx_zz, g_x_0_z_0_xxxxxx_xx, g_x_0_z_0_xxxxxx_xy, g_x_0_z_0_xxxxxx_xz, g_x_0_z_0_xxxxxx_yy, g_x_0_z_0_xxxxxx_yz, g_x_0_z_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxxx_xx[k] = -g_0_0_z_0_xxxxx_xx[k] - g_x_0_z_0_xxxxx_xx[k] * ab_x + g_x_0_z_0_xxxxx_xxx[k];

                g_x_0_z_0_xxxxxx_xy[k] = -g_0_0_z_0_xxxxx_xy[k] - g_x_0_z_0_xxxxx_xy[k] * ab_x + g_x_0_z_0_xxxxx_xxy[k];

                g_x_0_z_0_xxxxxx_xz[k] = -g_0_0_z_0_xxxxx_xz[k] - g_x_0_z_0_xxxxx_xz[k] * ab_x + g_x_0_z_0_xxxxx_xxz[k];

                g_x_0_z_0_xxxxxx_yy[k] = -g_0_0_z_0_xxxxx_yy[k] - g_x_0_z_0_xxxxx_yy[k] * ab_x + g_x_0_z_0_xxxxx_xyy[k];

                g_x_0_z_0_xxxxxx_yz[k] = -g_0_0_z_0_xxxxx_yz[k] - g_x_0_z_0_xxxxx_yz[k] * ab_x + g_x_0_z_0_xxxxx_xyz[k];

                g_x_0_z_0_xxxxxx_zz[k] = -g_0_0_z_0_xxxxx_zz[k] - g_x_0_z_0_xxxxx_zz[k] * ab_x + g_x_0_z_0_xxxxx_xzz[k];
            }

            /// Set up 342-348 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxxy_xx = cbuffer.data(id_geom_1010_off + 342 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxy_xy = cbuffer.data(id_geom_1010_off + 343 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxy_xz = cbuffer.data(id_geom_1010_off + 344 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxy_yy = cbuffer.data(id_geom_1010_off + 345 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxy_yz = cbuffer.data(id_geom_1010_off + 346 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxy_zz = cbuffer.data(id_geom_1010_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxx_xx, g_x_0_z_0_xxxxx_xxy, g_x_0_z_0_xxxxx_xy, g_x_0_z_0_xxxxx_xyy, g_x_0_z_0_xxxxx_xyz, g_x_0_z_0_xxxxx_xz, g_x_0_z_0_xxxxx_yy, g_x_0_z_0_xxxxx_yyy, g_x_0_z_0_xxxxx_yyz, g_x_0_z_0_xxxxx_yz, g_x_0_z_0_xxxxx_yzz, g_x_0_z_0_xxxxx_zz, g_x_0_z_0_xxxxxy_xx, g_x_0_z_0_xxxxxy_xy, g_x_0_z_0_xxxxxy_xz, g_x_0_z_0_xxxxxy_yy, g_x_0_z_0_xxxxxy_yz, g_x_0_z_0_xxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxxy_xx[k] = -g_x_0_z_0_xxxxx_xx[k] * ab_y + g_x_0_z_0_xxxxx_xxy[k];

                g_x_0_z_0_xxxxxy_xy[k] = -g_x_0_z_0_xxxxx_xy[k] * ab_y + g_x_0_z_0_xxxxx_xyy[k];

                g_x_0_z_0_xxxxxy_xz[k] = -g_x_0_z_0_xxxxx_xz[k] * ab_y + g_x_0_z_0_xxxxx_xyz[k];

                g_x_0_z_0_xxxxxy_yy[k] = -g_x_0_z_0_xxxxx_yy[k] * ab_y + g_x_0_z_0_xxxxx_yyy[k];

                g_x_0_z_0_xxxxxy_yz[k] = -g_x_0_z_0_xxxxx_yz[k] * ab_y + g_x_0_z_0_xxxxx_yyz[k];

                g_x_0_z_0_xxxxxy_zz[k] = -g_x_0_z_0_xxxxx_zz[k] * ab_y + g_x_0_z_0_xxxxx_yzz[k];
            }

            /// Set up 348-354 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxxz_xx = cbuffer.data(id_geom_1010_off + 348 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxz_xy = cbuffer.data(id_geom_1010_off + 349 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxz_xz = cbuffer.data(id_geom_1010_off + 350 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxz_yy = cbuffer.data(id_geom_1010_off + 351 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxz_yz = cbuffer.data(id_geom_1010_off + 352 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxz_zz = cbuffer.data(id_geom_1010_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxx_xx, g_x_0_z_0_xxxxx_xxz, g_x_0_z_0_xxxxx_xy, g_x_0_z_0_xxxxx_xyz, g_x_0_z_0_xxxxx_xz, g_x_0_z_0_xxxxx_xzz, g_x_0_z_0_xxxxx_yy, g_x_0_z_0_xxxxx_yyz, g_x_0_z_0_xxxxx_yz, g_x_0_z_0_xxxxx_yzz, g_x_0_z_0_xxxxx_zz, g_x_0_z_0_xxxxx_zzz, g_x_0_z_0_xxxxxz_xx, g_x_0_z_0_xxxxxz_xy, g_x_0_z_0_xxxxxz_xz, g_x_0_z_0_xxxxxz_yy, g_x_0_z_0_xxxxxz_yz, g_x_0_z_0_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxxz_xx[k] = -g_x_0_z_0_xxxxx_xx[k] * ab_z + g_x_0_z_0_xxxxx_xxz[k];

                g_x_0_z_0_xxxxxz_xy[k] = -g_x_0_z_0_xxxxx_xy[k] * ab_z + g_x_0_z_0_xxxxx_xyz[k];

                g_x_0_z_0_xxxxxz_xz[k] = -g_x_0_z_0_xxxxx_xz[k] * ab_z + g_x_0_z_0_xxxxx_xzz[k];

                g_x_0_z_0_xxxxxz_yy[k] = -g_x_0_z_0_xxxxx_yy[k] * ab_z + g_x_0_z_0_xxxxx_yyz[k];

                g_x_0_z_0_xxxxxz_yz[k] = -g_x_0_z_0_xxxxx_yz[k] * ab_z + g_x_0_z_0_xxxxx_yzz[k];

                g_x_0_z_0_xxxxxz_zz[k] = -g_x_0_z_0_xxxxx_zz[k] * ab_z + g_x_0_z_0_xxxxx_zzz[k];
            }

            /// Set up 354-360 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxyy_xx = cbuffer.data(id_geom_1010_off + 354 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyy_xy = cbuffer.data(id_geom_1010_off + 355 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyy_xz = cbuffer.data(id_geom_1010_off + 356 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyy_yy = cbuffer.data(id_geom_1010_off + 357 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyy_yz = cbuffer.data(id_geom_1010_off + 358 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyy_zz = cbuffer.data(id_geom_1010_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxy_xx, g_x_0_z_0_xxxxy_xxy, g_x_0_z_0_xxxxy_xy, g_x_0_z_0_xxxxy_xyy, g_x_0_z_0_xxxxy_xyz, g_x_0_z_0_xxxxy_xz, g_x_0_z_0_xxxxy_yy, g_x_0_z_0_xxxxy_yyy, g_x_0_z_0_xxxxy_yyz, g_x_0_z_0_xxxxy_yz, g_x_0_z_0_xxxxy_yzz, g_x_0_z_0_xxxxy_zz, g_x_0_z_0_xxxxyy_xx, g_x_0_z_0_xxxxyy_xy, g_x_0_z_0_xxxxyy_xz, g_x_0_z_0_xxxxyy_yy, g_x_0_z_0_xxxxyy_yz, g_x_0_z_0_xxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxyy_xx[k] = -g_x_0_z_0_xxxxy_xx[k] * ab_y + g_x_0_z_0_xxxxy_xxy[k];

                g_x_0_z_0_xxxxyy_xy[k] = -g_x_0_z_0_xxxxy_xy[k] * ab_y + g_x_0_z_0_xxxxy_xyy[k];

                g_x_0_z_0_xxxxyy_xz[k] = -g_x_0_z_0_xxxxy_xz[k] * ab_y + g_x_0_z_0_xxxxy_xyz[k];

                g_x_0_z_0_xxxxyy_yy[k] = -g_x_0_z_0_xxxxy_yy[k] * ab_y + g_x_0_z_0_xxxxy_yyy[k];

                g_x_0_z_0_xxxxyy_yz[k] = -g_x_0_z_0_xxxxy_yz[k] * ab_y + g_x_0_z_0_xxxxy_yyz[k];

                g_x_0_z_0_xxxxyy_zz[k] = -g_x_0_z_0_xxxxy_zz[k] * ab_y + g_x_0_z_0_xxxxy_yzz[k];
            }

            /// Set up 360-366 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxyz_xx = cbuffer.data(id_geom_1010_off + 360 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyz_xy = cbuffer.data(id_geom_1010_off + 361 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyz_xz = cbuffer.data(id_geom_1010_off + 362 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyz_yy = cbuffer.data(id_geom_1010_off + 363 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyz_yz = cbuffer.data(id_geom_1010_off + 364 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyz_zz = cbuffer.data(id_geom_1010_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxyz_xx, g_x_0_z_0_xxxxyz_xy, g_x_0_z_0_xxxxyz_xz, g_x_0_z_0_xxxxyz_yy, g_x_0_z_0_xxxxyz_yz, g_x_0_z_0_xxxxyz_zz, g_x_0_z_0_xxxxz_xx, g_x_0_z_0_xxxxz_xxy, g_x_0_z_0_xxxxz_xy, g_x_0_z_0_xxxxz_xyy, g_x_0_z_0_xxxxz_xyz, g_x_0_z_0_xxxxz_xz, g_x_0_z_0_xxxxz_yy, g_x_0_z_0_xxxxz_yyy, g_x_0_z_0_xxxxz_yyz, g_x_0_z_0_xxxxz_yz, g_x_0_z_0_xxxxz_yzz, g_x_0_z_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxyz_xx[k] = -g_x_0_z_0_xxxxz_xx[k] * ab_y + g_x_0_z_0_xxxxz_xxy[k];

                g_x_0_z_0_xxxxyz_xy[k] = -g_x_0_z_0_xxxxz_xy[k] * ab_y + g_x_0_z_0_xxxxz_xyy[k];

                g_x_0_z_0_xxxxyz_xz[k] = -g_x_0_z_0_xxxxz_xz[k] * ab_y + g_x_0_z_0_xxxxz_xyz[k];

                g_x_0_z_0_xxxxyz_yy[k] = -g_x_0_z_0_xxxxz_yy[k] * ab_y + g_x_0_z_0_xxxxz_yyy[k];

                g_x_0_z_0_xxxxyz_yz[k] = -g_x_0_z_0_xxxxz_yz[k] * ab_y + g_x_0_z_0_xxxxz_yyz[k];

                g_x_0_z_0_xxxxyz_zz[k] = -g_x_0_z_0_xxxxz_zz[k] * ab_y + g_x_0_z_0_xxxxz_yzz[k];
            }

            /// Set up 366-372 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxzz_xx = cbuffer.data(id_geom_1010_off + 366 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxzz_xy = cbuffer.data(id_geom_1010_off + 367 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxzz_xz = cbuffer.data(id_geom_1010_off + 368 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxzz_yy = cbuffer.data(id_geom_1010_off + 369 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxzz_yz = cbuffer.data(id_geom_1010_off + 370 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxzz_zz = cbuffer.data(id_geom_1010_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxz_xx, g_x_0_z_0_xxxxz_xxz, g_x_0_z_0_xxxxz_xy, g_x_0_z_0_xxxxz_xyz, g_x_0_z_0_xxxxz_xz, g_x_0_z_0_xxxxz_xzz, g_x_0_z_0_xxxxz_yy, g_x_0_z_0_xxxxz_yyz, g_x_0_z_0_xxxxz_yz, g_x_0_z_0_xxxxz_yzz, g_x_0_z_0_xxxxz_zz, g_x_0_z_0_xxxxz_zzz, g_x_0_z_0_xxxxzz_xx, g_x_0_z_0_xxxxzz_xy, g_x_0_z_0_xxxxzz_xz, g_x_0_z_0_xxxxzz_yy, g_x_0_z_0_xxxxzz_yz, g_x_0_z_0_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxzz_xx[k] = -g_x_0_z_0_xxxxz_xx[k] * ab_z + g_x_0_z_0_xxxxz_xxz[k];

                g_x_0_z_0_xxxxzz_xy[k] = -g_x_0_z_0_xxxxz_xy[k] * ab_z + g_x_0_z_0_xxxxz_xyz[k];

                g_x_0_z_0_xxxxzz_xz[k] = -g_x_0_z_0_xxxxz_xz[k] * ab_z + g_x_0_z_0_xxxxz_xzz[k];

                g_x_0_z_0_xxxxzz_yy[k] = -g_x_0_z_0_xxxxz_yy[k] * ab_z + g_x_0_z_0_xxxxz_yyz[k];

                g_x_0_z_0_xxxxzz_yz[k] = -g_x_0_z_0_xxxxz_yz[k] * ab_z + g_x_0_z_0_xxxxz_yzz[k];

                g_x_0_z_0_xxxxzz_zz[k] = -g_x_0_z_0_xxxxz_zz[k] * ab_z + g_x_0_z_0_xxxxz_zzz[k];
            }

            /// Set up 372-378 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxyyy_xx = cbuffer.data(id_geom_1010_off + 372 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyy_xy = cbuffer.data(id_geom_1010_off + 373 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyy_xz = cbuffer.data(id_geom_1010_off + 374 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyy_yy = cbuffer.data(id_geom_1010_off + 375 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyy_yz = cbuffer.data(id_geom_1010_off + 376 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyy_zz = cbuffer.data(id_geom_1010_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxyy_xx, g_x_0_z_0_xxxyy_xxy, g_x_0_z_0_xxxyy_xy, g_x_0_z_0_xxxyy_xyy, g_x_0_z_0_xxxyy_xyz, g_x_0_z_0_xxxyy_xz, g_x_0_z_0_xxxyy_yy, g_x_0_z_0_xxxyy_yyy, g_x_0_z_0_xxxyy_yyz, g_x_0_z_0_xxxyy_yz, g_x_0_z_0_xxxyy_yzz, g_x_0_z_0_xxxyy_zz, g_x_0_z_0_xxxyyy_xx, g_x_0_z_0_xxxyyy_xy, g_x_0_z_0_xxxyyy_xz, g_x_0_z_0_xxxyyy_yy, g_x_0_z_0_xxxyyy_yz, g_x_0_z_0_xxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxyyy_xx[k] = -g_x_0_z_0_xxxyy_xx[k] * ab_y + g_x_0_z_0_xxxyy_xxy[k];

                g_x_0_z_0_xxxyyy_xy[k] = -g_x_0_z_0_xxxyy_xy[k] * ab_y + g_x_0_z_0_xxxyy_xyy[k];

                g_x_0_z_0_xxxyyy_xz[k] = -g_x_0_z_0_xxxyy_xz[k] * ab_y + g_x_0_z_0_xxxyy_xyz[k];

                g_x_0_z_0_xxxyyy_yy[k] = -g_x_0_z_0_xxxyy_yy[k] * ab_y + g_x_0_z_0_xxxyy_yyy[k];

                g_x_0_z_0_xxxyyy_yz[k] = -g_x_0_z_0_xxxyy_yz[k] * ab_y + g_x_0_z_0_xxxyy_yyz[k];

                g_x_0_z_0_xxxyyy_zz[k] = -g_x_0_z_0_xxxyy_zz[k] * ab_y + g_x_0_z_0_xxxyy_yzz[k];
            }

            /// Set up 378-384 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxyyz_xx = cbuffer.data(id_geom_1010_off + 378 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyz_xy = cbuffer.data(id_geom_1010_off + 379 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyz_xz = cbuffer.data(id_geom_1010_off + 380 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyz_yy = cbuffer.data(id_geom_1010_off + 381 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyz_yz = cbuffer.data(id_geom_1010_off + 382 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyz_zz = cbuffer.data(id_geom_1010_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxyyz_xx, g_x_0_z_0_xxxyyz_xy, g_x_0_z_0_xxxyyz_xz, g_x_0_z_0_xxxyyz_yy, g_x_0_z_0_xxxyyz_yz, g_x_0_z_0_xxxyyz_zz, g_x_0_z_0_xxxyz_xx, g_x_0_z_0_xxxyz_xxy, g_x_0_z_0_xxxyz_xy, g_x_0_z_0_xxxyz_xyy, g_x_0_z_0_xxxyz_xyz, g_x_0_z_0_xxxyz_xz, g_x_0_z_0_xxxyz_yy, g_x_0_z_0_xxxyz_yyy, g_x_0_z_0_xxxyz_yyz, g_x_0_z_0_xxxyz_yz, g_x_0_z_0_xxxyz_yzz, g_x_0_z_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxyyz_xx[k] = -g_x_0_z_0_xxxyz_xx[k] * ab_y + g_x_0_z_0_xxxyz_xxy[k];

                g_x_0_z_0_xxxyyz_xy[k] = -g_x_0_z_0_xxxyz_xy[k] * ab_y + g_x_0_z_0_xxxyz_xyy[k];

                g_x_0_z_0_xxxyyz_xz[k] = -g_x_0_z_0_xxxyz_xz[k] * ab_y + g_x_0_z_0_xxxyz_xyz[k];

                g_x_0_z_0_xxxyyz_yy[k] = -g_x_0_z_0_xxxyz_yy[k] * ab_y + g_x_0_z_0_xxxyz_yyy[k];

                g_x_0_z_0_xxxyyz_yz[k] = -g_x_0_z_0_xxxyz_yz[k] * ab_y + g_x_0_z_0_xxxyz_yyz[k];

                g_x_0_z_0_xxxyyz_zz[k] = -g_x_0_z_0_xxxyz_zz[k] * ab_y + g_x_0_z_0_xxxyz_yzz[k];
            }

            /// Set up 384-390 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxyzz_xx = cbuffer.data(id_geom_1010_off + 384 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyzz_xy = cbuffer.data(id_geom_1010_off + 385 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyzz_xz = cbuffer.data(id_geom_1010_off + 386 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyzz_yy = cbuffer.data(id_geom_1010_off + 387 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyzz_yz = cbuffer.data(id_geom_1010_off + 388 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyzz_zz = cbuffer.data(id_geom_1010_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxyzz_xx, g_x_0_z_0_xxxyzz_xy, g_x_0_z_0_xxxyzz_xz, g_x_0_z_0_xxxyzz_yy, g_x_0_z_0_xxxyzz_yz, g_x_0_z_0_xxxyzz_zz, g_x_0_z_0_xxxzz_xx, g_x_0_z_0_xxxzz_xxy, g_x_0_z_0_xxxzz_xy, g_x_0_z_0_xxxzz_xyy, g_x_0_z_0_xxxzz_xyz, g_x_0_z_0_xxxzz_xz, g_x_0_z_0_xxxzz_yy, g_x_0_z_0_xxxzz_yyy, g_x_0_z_0_xxxzz_yyz, g_x_0_z_0_xxxzz_yz, g_x_0_z_0_xxxzz_yzz, g_x_0_z_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxyzz_xx[k] = -g_x_0_z_0_xxxzz_xx[k] * ab_y + g_x_0_z_0_xxxzz_xxy[k];

                g_x_0_z_0_xxxyzz_xy[k] = -g_x_0_z_0_xxxzz_xy[k] * ab_y + g_x_0_z_0_xxxzz_xyy[k];

                g_x_0_z_0_xxxyzz_xz[k] = -g_x_0_z_0_xxxzz_xz[k] * ab_y + g_x_0_z_0_xxxzz_xyz[k];

                g_x_0_z_0_xxxyzz_yy[k] = -g_x_0_z_0_xxxzz_yy[k] * ab_y + g_x_0_z_0_xxxzz_yyy[k];

                g_x_0_z_0_xxxyzz_yz[k] = -g_x_0_z_0_xxxzz_yz[k] * ab_y + g_x_0_z_0_xxxzz_yyz[k];

                g_x_0_z_0_xxxyzz_zz[k] = -g_x_0_z_0_xxxzz_zz[k] * ab_y + g_x_0_z_0_xxxzz_yzz[k];
            }

            /// Set up 390-396 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxzzz_xx = cbuffer.data(id_geom_1010_off + 390 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzzz_xy = cbuffer.data(id_geom_1010_off + 391 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzzz_xz = cbuffer.data(id_geom_1010_off + 392 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzzz_yy = cbuffer.data(id_geom_1010_off + 393 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzzz_yz = cbuffer.data(id_geom_1010_off + 394 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzzz_zz = cbuffer.data(id_geom_1010_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxzz_xx, g_x_0_z_0_xxxzz_xxz, g_x_0_z_0_xxxzz_xy, g_x_0_z_0_xxxzz_xyz, g_x_0_z_0_xxxzz_xz, g_x_0_z_0_xxxzz_xzz, g_x_0_z_0_xxxzz_yy, g_x_0_z_0_xxxzz_yyz, g_x_0_z_0_xxxzz_yz, g_x_0_z_0_xxxzz_yzz, g_x_0_z_0_xxxzz_zz, g_x_0_z_0_xxxzz_zzz, g_x_0_z_0_xxxzzz_xx, g_x_0_z_0_xxxzzz_xy, g_x_0_z_0_xxxzzz_xz, g_x_0_z_0_xxxzzz_yy, g_x_0_z_0_xxxzzz_yz, g_x_0_z_0_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxzzz_xx[k] = -g_x_0_z_0_xxxzz_xx[k] * ab_z + g_x_0_z_0_xxxzz_xxz[k];

                g_x_0_z_0_xxxzzz_xy[k] = -g_x_0_z_0_xxxzz_xy[k] * ab_z + g_x_0_z_0_xxxzz_xyz[k];

                g_x_0_z_0_xxxzzz_xz[k] = -g_x_0_z_0_xxxzz_xz[k] * ab_z + g_x_0_z_0_xxxzz_xzz[k];

                g_x_0_z_0_xxxzzz_yy[k] = -g_x_0_z_0_xxxzz_yy[k] * ab_z + g_x_0_z_0_xxxzz_yyz[k];

                g_x_0_z_0_xxxzzz_yz[k] = -g_x_0_z_0_xxxzz_yz[k] * ab_z + g_x_0_z_0_xxxzz_yzz[k];

                g_x_0_z_0_xxxzzz_zz[k] = -g_x_0_z_0_xxxzz_zz[k] * ab_z + g_x_0_z_0_xxxzz_zzz[k];
            }

            /// Set up 396-402 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyyyy_xx = cbuffer.data(id_geom_1010_off + 396 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyy_xy = cbuffer.data(id_geom_1010_off + 397 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyy_xz = cbuffer.data(id_geom_1010_off + 398 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyy_yy = cbuffer.data(id_geom_1010_off + 399 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyy_yz = cbuffer.data(id_geom_1010_off + 400 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyy_zz = cbuffer.data(id_geom_1010_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyyy_xx, g_x_0_z_0_xxyyy_xxy, g_x_0_z_0_xxyyy_xy, g_x_0_z_0_xxyyy_xyy, g_x_0_z_0_xxyyy_xyz, g_x_0_z_0_xxyyy_xz, g_x_0_z_0_xxyyy_yy, g_x_0_z_0_xxyyy_yyy, g_x_0_z_0_xxyyy_yyz, g_x_0_z_0_xxyyy_yz, g_x_0_z_0_xxyyy_yzz, g_x_0_z_0_xxyyy_zz, g_x_0_z_0_xxyyyy_xx, g_x_0_z_0_xxyyyy_xy, g_x_0_z_0_xxyyyy_xz, g_x_0_z_0_xxyyyy_yy, g_x_0_z_0_xxyyyy_yz, g_x_0_z_0_xxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyyyy_xx[k] = -g_x_0_z_0_xxyyy_xx[k] * ab_y + g_x_0_z_0_xxyyy_xxy[k];

                g_x_0_z_0_xxyyyy_xy[k] = -g_x_0_z_0_xxyyy_xy[k] * ab_y + g_x_0_z_0_xxyyy_xyy[k];

                g_x_0_z_0_xxyyyy_xz[k] = -g_x_0_z_0_xxyyy_xz[k] * ab_y + g_x_0_z_0_xxyyy_xyz[k];

                g_x_0_z_0_xxyyyy_yy[k] = -g_x_0_z_0_xxyyy_yy[k] * ab_y + g_x_0_z_0_xxyyy_yyy[k];

                g_x_0_z_0_xxyyyy_yz[k] = -g_x_0_z_0_xxyyy_yz[k] * ab_y + g_x_0_z_0_xxyyy_yyz[k];

                g_x_0_z_0_xxyyyy_zz[k] = -g_x_0_z_0_xxyyy_zz[k] * ab_y + g_x_0_z_0_xxyyy_yzz[k];
            }

            /// Set up 402-408 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyyyz_xx = cbuffer.data(id_geom_1010_off + 402 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyz_xy = cbuffer.data(id_geom_1010_off + 403 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyz_xz = cbuffer.data(id_geom_1010_off + 404 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyz_yy = cbuffer.data(id_geom_1010_off + 405 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyz_yz = cbuffer.data(id_geom_1010_off + 406 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyz_zz = cbuffer.data(id_geom_1010_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyyyz_xx, g_x_0_z_0_xxyyyz_xy, g_x_0_z_0_xxyyyz_xz, g_x_0_z_0_xxyyyz_yy, g_x_0_z_0_xxyyyz_yz, g_x_0_z_0_xxyyyz_zz, g_x_0_z_0_xxyyz_xx, g_x_0_z_0_xxyyz_xxy, g_x_0_z_0_xxyyz_xy, g_x_0_z_0_xxyyz_xyy, g_x_0_z_0_xxyyz_xyz, g_x_0_z_0_xxyyz_xz, g_x_0_z_0_xxyyz_yy, g_x_0_z_0_xxyyz_yyy, g_x_0_z_0_xxyyz_yyz, g_x_0_z_0_xxyyz_yz, g_x_0_z_0_xxyyz_yzz, g_x_0_z_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyyyz_xx[k] = -g_x_0_z_0_xxyyz_xx[k] * ab_y + g_x_0_z_0_xxyyz_xxy[k];

                g_x_0_z_0_xxyyyz_xy[k] = -g_x_0_z_0_xxyyz_xy[k] * ab_y + g_x_0_z_0_xxyyz_xyy[k];

                g_x_0_z_0_xxyyyz_xz[k] = -g_x_0_z_0_xxyyz_xz[k] * ab_y + g_x_0_z_0_xxyyz_xyz[k];

                g_x_0_z_0_xxyyyz_yy[k] = -g_x_0_z_0_xxyyz_yy[k] * ab_y + g_x_0_z_0_xxyyz_yyy[k];

                g_x_0_z_0_xxyyyz_yz[k] = -g_x_0_z_0_xxyyz_yz[k] * ab_y + g_x_0_z_0_xxyyz_yyz[k];

                g_x_0_z_0_xxyyyz_zz[k] = -g_x_0_z_0_xxyyz_zz[k] * ab_y + g_x_0_z_0_xxyyz_yzz[k];
            }

            /// Set up 408-414 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyyzz_xx = cbuffer.data(id_geom_1010_off + 408 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyzz_xy = cbuffer.data(id_geom_1010_off + 409 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyzz_xz = cbuffer.data(id_geom_1010_off + 410 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyzz_yy = cbuffer.data(id_geom_1010_off + 411 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyzz_yz = cbuffer.data(id_geom_1010_off + 412 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyzz_zz = cbuffer.data(id_geom_1010_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyyzz_xx, g_x_0_z_0_xxyyzz_xy, g_x_0_z_0_xxyyzz_xz, g_x_0_z_0_xxyyzz_yy, g_x_0_z_0_xxyyzz_yz, g_x_0_z_0_xxyyzz_zz, g_x_0_z_0_xxyzz_xx, g_x_0_z_0_xxyzz_xxy, g_x_0_z_0_xxyzz_xy, g_x_0_z_0_xxyzz_xyy, g_x_0_z_0_xxyzz_xyz, g_x_0_z_0_xxyzz_xz, g_x_0_z_0_xxyzz_yy, g_x_0_z_0_xxyzz_yyy, g_x_0_z_0_xxyzz_yyz, g_x_0_z_0_xxyzz_yz, g_x_0_z_0_xxyzz_yzz, g_x_0_z_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyyzz_xx[k] = -g_x_0_z_0_xxyzz_xx[k] * ab_y + g_x_0_z_0_xxyzz_xxy[k];

                g_x_0_z_0_xxyyzz_xy[k] = -g_x_0_z_0_xxyzz_xy[k] * ab_y + g_x_0_z_0_xxyzz_xyy[k];

                g_x_0_z_0_xxyyzz_xz[k] = -g_x_0_z_0_xxyzz_xz[k] * ab_y + g_x_0_z_0_xxyzz_xyz[k];

                g_x_0_z_0_xxyyzz_yy[k] = -g_x_0_z_0_xxyzz_yy[k] * ab_y + g_x_0_z_0_xxyzz_yyy[k];

                g_x_0_z_0_xxyyzz_yz[k] = -g_x_0_z_0_xxyzz_yz[k] * ab_y + g_x_0_z_0_xxyzz_yyz[k];

                g_x_0_z_0_xxyyzz_zz[k] = -g_x_0_z_0_xxyzz_zz[k] * ab_y + g_x_0_z_0_xxyzz_yzz[k];
            }

            /// Set up 414-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyzzz_xx = cbuffer.data(id_geom_1010_off + 414 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzzz_xy = cbuffer.data(id_geom_1010_off + 415 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzzz_xz = cbuffer.data(id_geom_1010_off + 416 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzzz_yy = cbuffer.data(id_geom_1010_off + 417 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzzz_yz = cbuffer.data(id_geom_1010_off + 418 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzzz_zz = cbuffer.data(id_geom_1010_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyzzz_xx, g_x_0_z_0_xxyzzz_xy, g_x_0_z_0_xxyzzz_xz, g_x_0_z_0_xxyzzz_yy, g_x_0_z_0_xxyzzz_yz, g_x_0_z_0_xxyzzz_zz, g_x_0_z_0_xxzzz_xx, g_x_0_z_0_xxzzz_xxy, g_x_0_z_0_xxzzz_xy, g_x_0_z_0_xxzzz_xyy, g_x_0_z_0_xxzzz_xyz, g_x_0_z_0_xxzzz_xz, g_x_0_z_0_xxzzz_yy, g_x_0_z_0_xxzzz_yyy, g_x_0_z_0_xxzzz_yyz, g_x_0_z_0_xxzzz_yz, g_x_0_z_0_xxzzz_yzz, g_x_0_z_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyzzz_xx[k] = -g_x_0_z_0_xxzzz_xx[k] * ab_y + g_x_0_z_0_xxzzz_xxy[k];

                g_x_0_z_0_xxyzzz_xy[k] = -g_x_0_z_0_xxzzz_xy[k] * ab_y + g_x_0_z_0_xxzzz_xyy[k];

                g_x_0_z_0_xxyzzz_xz[k] = -g_x_0_z_0_xxzzz_xz[k] * ab_y + g_x_0_z_0_xxzzz_xyz[k];

                g_x_0_z_0_xxyzzz_yy[k] = -g_x_0_z_0_xxzzz_yy[k] * ab_y + g_x_0_z_0_xxzzz_yyy[k];

                g_x_0_z_0_xxyzzz_yz[k] = -g_x_0_z_0_xxzzz_yz[k] * ab_y + g_x_0_z_0_xxzzz_yyz[k];

                g_x_0_z_0_xxyzzz_zz[k] = -g_x_0_z_0_xxzzz_zz[k] * ab_y + g_x_0_z_0_xxzzz_yzz[k];
            }

            /// Set up 420-426 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxzzzz_xx = cbuffer.data(id_geom_1010_off + 420 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzzz_xy = cbuffer.data(id_geom_1010_off + 421 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzzz_xz = cbuffer.data(id_geom_1010_off + 422 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzzz_yy = cbuffer.data(id_geom_1010_off + 423 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzzz_yz = cbuffer.data(id_geom_1010_off + 424 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzzz_zz = cbuffer.data(id_geom_1010_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxzzz_xx, g_x_0_z_0_xxzzz_xxz, g_x_0_z_0_xxzzz_xy, g_x_0_z_0_xxzzz_xyz, g_x_0_z_0_xxzzz_xz, g_x_0_z_0_xxzzz_xzz, g_x_0_z_0_xxzzz_yy, g_x_0_z_0_xxzzz_yyz, g_x_0_z_0_xxzzz_yz, g_x_0_z_0_xxzzz_yzz, g_x_0_z_0_xxzzz_zz, g_x_0_z_0_xxzzz_zzz, g_x_0_z_0_xxzzzz_xx, g_x_0_z_0_xxzzzz_xy, g_x_0_z_0_xxzzzz_xz, g_x_0_z_0_xxzzzz_yy, g_x_0_z_0_xxzzzz_yz, g_x_0_z_0_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxzzzz_xx[k] = -g_x_0_z_0_xxzzz_xx[k] * ab_z + g_x_0_z_0_xxzzz_xxz[k];

                g_x_0_z_0_xxzzzz_xy[k] = -g_x_0_z_0_xxzzz_xy[k] * ab_z + g_x_0_z_0_xxzzz_xyz[k];

                g_x_0_z_0_xxzzzz_xz[k] = -g_x_0_z_0_xxzzz_xz[k] * ab_z + g_x_0_z_0_xxzzz_xzz[k];

                g_x_0_z_0_xxzzzz_yy[k] = -g_x_0_z_0_xxzzz_yy[k] * ab_z + g_x_0_z_0_xxzzz_yyz[k];

                g_x_0_z_0_xxzzzz_yz[k] = -g_x_0_z_0_xxzzz_yz[k] * ab_z + g_x_0_z_0_xxzzz_yzz[k];

                g_x_0_z_0_xxzzzz_zz[k] = -g_x_0_z_0_xxzzz_zz[k] * ab_z + g_x_0_z_0_xxzzz_zzz[k];
            }

            /// Set up 426-432 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyyyy_xx = cbuffer.data(id_geom_1010_off + 426 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyy_xy = cbuffer.data(id_geom_1010_off + 427 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyy_xz = cbuffer.data(id_geom_1010_off + 428 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyy_yy = cbuffer.data(id_geom_1010_off + 429 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyy_yz = cbuffer.data(id_geom_1010_off + 430 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyy_zz = cbuffer.data(id_geom_1010_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyyy_xx, g_x_0_z_0_xyyyy_xxy, g_x_0_z_0_xyyyy_xy, g_x_0_z_0_xyyyy_xyy, g_x_0_z_0_xyyyy_xyz, g_x_0_z_0_xyyyy_xz, g_x_0_z_0_xyyyy_yy, g_x_0_z_0_xyyyy_yyy, g_x_0_z_0_xyyyy_yyz, g_x_0_z_0_xyyyy_yz, g_x_0_z_0_xyyyy_yzz, g_x_0_z_0_xyyyy_zz, g_x_0_z_0_xyyyyy_xx, g_x_0_z_0_xyyyyy_xy, g_x_0_z_0_xyyyyy_xz, g_x_0_z_0_xyyyyy_yy, g_x_0_z_0_xyyyyy_yz, g_x_0_z_0_xyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyyyy_xx[k] = -g_x_0_z_0_xyyyy_xx[k] * ab_y + g_x_0_z_0_xyyyy_xxy[k];

                g_x_0_z_0_xyyyyy_xy[k] = -g_x_0_z_0_xyyyy_xy[k] * ab_y + g_x_0_z_0_xyyyy_xyy[k];

                g_x_0_z_0_xyyyyy_xz[k] = -g_x_0_z_0_xyyyy_xz[k] * ab_y + g_x_0_z_0_xyyyy_xyz[k];

                g_x_0_z_0_xyyyyy_yy[k] = -g_x_0_z_0_xyyyy_yy[k] * ab_y + g_x_0_z_0_xyyyy_yyy[k];

                g_x_0_z_0_xyyyyy_yz[k] = -g_x_0_z_0_xyyyy_yz[k] * ab_y + g_x_0_z_0_xyyyy_yyz[k];

                g_x_0_z_0_xyyyyy_zz[k] = -g_x_0_z_0_xyyyy_zz[k] * ab_y + g_x_0_z_0_xyyyy_yzz[k];
            }

            /// Set up 432-438 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyyyz_xx = cbuffer.data(id_geom_1010_off + 432 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyz_xy = cbuffer.data(id_geom_1010_off + 433 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyz_xz = cbuffer.data(id_geom_1010_off + 434 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyz_yy = cbuffer.data(id_geom_1010_off + 435 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyz_yz = cbuffer.data(id_geom_1010_off + 436 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyz_zz = cbuffer.data(id_geom_1010_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyyyz_xx, g_x_0_z_0_xyyyyz_xy, g_x_0_z_0_xyyyyz_xz, g_x_0_z_0_xyyyyz_yy, g_x_0_z_0_xyyyyz_yz, g_x_0_z_0_xyyyyz_zz, g_x_0_z_0_xyyyz_xx, g_x_0_z_0_xyyyz_xxy, g_x_0_z_0_xyyyz_xy, g_x_0_z_0_xyyyz_xyy, g_x_0_z_0_xyyyz_xyz, g_x_0_z_0_xyyyz_xz, g_x_0_z_0_xyyyz_yy, g_x_0_z_0_xyyyz_yyy, g_x_0_z_0_xyyyz_yyz, g_x_0_z_0_xyyyz_yz, g_x_0_z_0_xyyyz_yzz, g_x_0_z_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyyyz_xx[k] = -g_x_0_z_0_xyyyz_xx[k] * ab_y + g_x_0_z_0_xyyyz_xxy[k];

                g_x_0_z_0_xyyyyz_xy[k] = -g_x_0_z_0_xyyyz_xy[k] * ab_y + g_x_0_z_0_xyyyz_xyy[k];

                g_x_0_z_0_xyyyyz_xz[k] = -g_x_0_z_0_xyyyz_xz[k] * ab_y + g_x_0_z_0_xyyyz_xyz[k];

                g_x_0_z_0_xyyyyz_yy[k] = -g_x_0_z_0_xyyyz_yy[k] * ab_y + g_x_0_z_0_xyyyz_yyy[k];

                g_x_0_z_0_xyyyyz_yz[k] = -g_x_0_z_0_xyyyz_yz[k] * ab_y + g_x_0_z_0_xyyyz_yyz[k];

                g_x_0_z_0_xyyyyz_zz[k] = -g_x_0_z_0_xyyyz_zz[k] * ab_y + g_x_0_z_0_xyyyz_yzz[k];
            }

            /// Set up 438-444 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyyzz_xx = cbuffer.data(id_geom_1010_off + 438 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyzz_xy = cbuffer.data(id_geom_1010_off + 439 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyzz_xz = cbuffer.data(id_geom_1010_off + 440 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyzz_yy = cbuffer.data(id_geom_1010_off + 441 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyzz_yz = cbuffer.data(id_geom_1010_off + 442 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyzz_zz = cbuffer.data(id_geom_1010_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyyzz_xx, g_x_0_z_0_xyyyzz_xy, g_x_0_z_0_xyyyzz_xz, g_x_0_z_0_xyyyzz_yy, g_x_0_z_0_xyyyzz_yz, g_x_0_z_0_xyyyzz_zz, g_x_0_z_0_xyyzz_xx, g_x_0_z_0_xyyzz_xxy, g_x_0_z_0_xyyzz_xy, g_x_0_z_0_xyyzz_xyy, g_x_0_z_0_xyyzz_xyz, g_x_0_z_0_xyyzz_xz, g_x_0_z_0_xyyzz_yy, g_x_0_z_0_xyyzz_yyy, g_x_0_z_0_xyyzz_yyz, g_x_0_z_0_xyyzz_yz, g_x_0_z_0_xyyzz_yzz, g_x_0_z_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyyzz_xx[k] = -g_x_0_z_0_xyyzz_xx[k] * ab_y + g_x_0_z_0_xyyzz_xxy[k];

                g_x_0_z_0_xyyyzz_xy[k] = -g_x_0_z_0_xyyzz_xy[k] * ab_y + g_x_0_z_0_xyyzz_xyy[k];

                g_x_0_z_0_xyyyzz_xz[k] = -g_x_0_z_0_xyyzz_xz[k] * ab_y + g_x_0_z_0_xyyzz_xyz[k];

                g_x_0_z_0_xyyyzz_yy[k] = -g_x_0_z_0_xyyzz_yy[k] * ab_y + g_x_0_z_0_xyyzz_yyy[k];

                g_x_0_z_0_xyyyzz_yz[k] = -g_x_0_z_0_xyyzz_yz[k] * ab_y + g_x_0_z_0_xyyzz_yyz[k];

                g_x_0_z_0_xyyyzz_zz[k] = -g_x_0_z_0_xyyzz_zz[k] * ab_y + g_x_0_z_0_xyyzz_yzz[k];
            }

            /// Set up 444-450 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyzzz_xx = cbuffer.data(id_geom_1010_off + 444 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzzz_xy = cbuffer.data(id_geom_1010_off + 445 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzzz_xz = cbuffer.data(id_geom_1010_off + 446 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzzz_yy = cbuffer.data(id_geom_1010_off + 447 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzzz_yz = cbuffer.data(id_geom_1010_off + 448 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzzz_zz = cbuffer.data(id_geom_1010_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyzzz_xx, g_x_0_z_0_xyyzzz_xy, g_x_0_z_0_xyyzzz_xz, g_x_0_z_0_xyyzzz_yy, g_x_0_z_0_xyyzzz_yz, g_x_0_z_0_xyyzzz_zz, g_x_0_z_0_xyzzz_xx, g_x_0_z_0_xyzzz_xxy, g_x_0_z_0_xyzzz_xy, g_x_0_z_0_xyzzz_xyy, g_x_0_z_0_xyzzz_xyz, g_x_0_z_0_xyzzz_xz, g_x_0_z_0_xyzzz_yy, g_x_0_z_0_xyzzz_yyy, g_x_0_z_0_xyzzz_yyz, g_x_0_z_0_xyzzz_yz, g_x_0_z_0_xyzzz_yzz, g_x_0_z_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyzzz_xx[k] = -g_x_0_z_0_xyzzz_xx[k] * ab_y + g_x_0_z_0_xyzzz_xxy[k];

                g_x_0_z_0_xyyzzz_xy[k] = -g_x_0_z_0_xyzzz_xy[k] * ab_y + g_x_0_z_0_xyzzz_xyy[k];

                g_x_0_z_0_xyyzzz_xz[k] = -g_x_0_z_0_xyzzz_xz[k] * ab_y + g_x_0_z_0_xyzzz_xyz[k];

                g_x_0_z_0_xyyzzz_yy[k] = -g_x_0_z_0_xyzzz_yy[k] * ab_y + g_x_0_z_0_xyzzz_yyy[k];

                g_x_0_z_0_xyyzzz_yz[k] = -g_x_0_z_0_xyzzz_yz[k] * ab_y + g_x_0_z_0_xyzzz_yyz[k];

                g_x_0_z_0_xyyzzz_zz[k] = -g_x_0_z_0_xyzzz_zz[k] * ab_y + g_x_0_z_0_xyzzz_yzz[k];
            }

            /// Set up 450-456 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyzzzz_xx = cbuffer.data(id_geom_1010_off + 450 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzzz_xy = cbuffer.data(id_geom_1010_off + 451 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzzz_xz = cbuffer.data(id_geom_1010_off + 452 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzzz_yy = cbuffer.data(id_geom_1010_off + 453 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzzz_yz = cbuffer.data(id_geom_1010_off + 454 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzzz_zz = cbuffer.data(id_geom_1010_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyzzzz_xx, g_x_0_z_0_xyzzzz_xy, g_x_0_z_0_xyzzzz_xz, g_x_0_z_0_xyzzzz_yy, g_x_0_z_0_xyzzzz_yz, g_x_0_z_0_xyzzzz_zz, g_x_0_z_0_xzzzz_xx, g_x_0_z_0_xzzzz_xxy, g_x_0_z_0_xzzzz_xy, g_x_0_z_0_xzzzz_xyy, g_x_0_z_0_xzzzz_xyz, g_x_0_z_0_xzzzz_xz, g_x_0_z_0_xzzzz_yy, g_x_0_z_0_xzzzz_yyy, g_x_0_z_0_xzzzz_yyz, g_x_0_z_0_xzzzz_yz, g_x_0_z_0_xzzzz_yzz, g_x_0_z_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyzzzz_xx[k] = -g_x_0_z_0_xzzzz_xx[k] * ab_y + g_x_0_z_0_xzzzz_xxy[k];

                g_x_0_z_0_xyzzzz_xy[k] = -g_x_0_z_0_xzzzz_xy[k] * ab_y + g_x_0_z_0_xzzzz_xyy[k];

                g_x_0_z_0_xyzzzz_xz[k] = -g_x_0_z_0_xzzzz_xz[k] * ab_y + g_x_0_z_0_xzzzz_xyz[k];

                g_x_0_z_0_xyzzzz_yy[k] = -g_x_0_z_0_xzzzz_yy[k] * ab_y + g_x_0_z_0_xzzzz_yyy[k];

                g_x_0_z_0_xyzzzz_yz[k] = -g_x_0_z_0_xzzzz_yz[k] * ab_y + g_x_0_z_0_xzzzz_yyz[k];

                g_x_0_z_0_xyzzzz_zz[k] = -g_x_0_z_0_xzzzz_zz[k] * ab_y + g_x_0_z_0_xzzzz_yzz[k];
            }

            /// Set up 456-462 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xzzzzz_xx = cbuffer.data(id_geom_1010_off + 456 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzzz_xy = cbuffer.data(id_geom_1010_off + 457 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzzz_xz = cbuffer.data(id_geom_1010_off + 458 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzzz_yy = cbuffer.data(id_geom_1010_off + 459 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzzz_yz = cbuffer.data(id_geom_1010_off + 460 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzzz_zz = cbuffer.data(id_geom_1010_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xzzzz_xx, g_x_0_z_0_xzzzz_xxz, g_x_0_z_0_xzzzz_xy, g_x_0_z_0_xzzzz_xyz, g_x_0_z_0_xzzzz_xz, g_x_0_z_0_xzzzz_xzz, g_x_0_z_0_xzzzz_yy, g_x_0_z_0_xzzzz_yyz, g_x_0_z_0_xzzzz_yz, g_x_0_z_0_xzzzz_yzz, g_x_0_z_0_xzzzz_zz, g_x_0_z_0_xzzzz_zzz, g_x_0_z_0_xzzzzz_xx, g_x_0_z_0_xzzzzz_xy, g_x_0_z_0_xzzzzz_xz, g_x_0_z_0_xzzzzz_yy, g_x_0_z_0_xzzzzz_yz, g_x_0_z_0_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xzzzzz_xx[k] = -g_x_0_z_0_xzzzz_xx[k] * ab_z + g_x_0_z_0_xzzzz_xxz[k];

                g_x_0_z_0_xzzzzz_xy[k] = -g_x_0_z_0_xzzzz_xy[k] * ab_z + g_x_0_z_0_xzzzz_xyz[k];

                g_x_0_z_0_xzzzzz_xz[k] = -g_x_0_z_0_xzzzz_xz[k] * ab_z + g_x_0_z_0_xzzzz_xzz[k];

                g_x_0_z_0_xzzzzz_yy[k] = -g_x_0_z_0_xzzzz_yy[k] * ab_z + g_x_0_z_0_xzzzz_yyz[k];

                g_x_0_z_0_xzzzzz_yz[k] = -g_x_0_z_0_xzzzz_yz[k] * ab_z + g_x_0_z_0_xzzzz_yzz[k];

                g_x_0_z_0_xzzzzz_zz[k] = -g_x_0_z_0_xzzzz_zz[k] * ab_z + g_x_0_z_0_xzzzz_zzz[k];
            }

            /// Set up 462-468 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyyyy_xx = cbuffer.data(id_geom_1010_off + 462 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyy_xy = cbuffer.data(id_geom_1010_off + 463 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyy_xz = cbuffer.data(id_geom_1010_off + 464 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyy_yy = cbuffer.data(id_geom_1010_off + 465 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyy_yz = cbuffer.data(id_geom_1010_off + 466 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyy_zz = cbuffer.data(id_geom_1010_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyyy_xx, g_x_0_z_0_yyyyy_xxy, g_x_0_z_0_yyyyy_xy, g_x_0_z_0_yyyyy_xyy, g_x_0_z_0_yyyyy_xyz, g_x_0_z_0_yyyyy_xz, g_x_0_z_0_yyyyy_yy, g_x_0_z_0_yyyyy_yyy, g_x_0_z_0_yyyyy_yyz, g_x_0_z_0_yyyyy_yz, g_x_0_z_0_yyyyy_yzz, g_x_0_z_0_yyyyy_zz, g_x_0_z_0_yyyyyy_xx, g_x_0_z_0_yyyyyy_xy, g_x_0_z_0_yyyyyy_xz, g_x_0_z_0_yyyyyy_yy, g_x_0_z_0_yyyyyy_yz, g_x_0_z_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyyyy_xx[k] = -g_x_0_z_0_yyyyy_xx[k] * ab_y + g_x_0_z_0_yyyyy_xxy[k];

                g_x_0_z_0_yyyyyy_xy[k] = -g_x_0_z_0_yyyyy_xy[k] * ab_y + g_x_0_z_0_yyyyy_xyy[k];

                g_x_0_z_0_yyyyyy_xz[k] = -g_x_0_z_0_yyyyy_xz[k] * ab_y + g_x_0_z_0_yyyyy_xyz[k];

                g_x_0_z_0_yyyyyy_yy[k] = -g_x_0_z_0_yyyyy_yy[k] * ab_y + g_x_0_z_0_yyyyy_yyy[k];

                g_x_0_z_0_yyyyyy_yz[k] = -g_x_0_z_0_yyyyy_yz[k] * ab_y + g_x_0_z_0_yyyyy_yyz[k];

                g_x_0_z_0_yyyyyy_zz[k] = -g_x_0_z_0_yyyyy_zz[k] * ab_y + g_x_0_z_0_yyyyy_yzz[k];
            }

            /// Set up 468-474 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyyyz_xx = cbuffer.data(id_geom_1010_off + 468 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyz_xy = cbuffer.data(id_geom_1010_off + 469 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyz_xz = cbuffer.data(id_geom_1010_off + 470 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyz_yy = cbuffer.data(id_geom_1010_off + 471 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyz_yz = cbuffer.data(id_geom_1010_off + 472 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyz_zz = cbuffer.data(id_geom_1010_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyyyz_xx, g_x_0_z_0_yyyyyz_xy, g_x_0_z_0_yyyyyz_xz, g_x_0_z_0_yyyyyz_yy, g_x_0_z_0_yyyyyz_yz, g_x_0_z_0_yyyyyz_zz, g_x_0_z_0_yyyyz_xx, g_x_0_z_0_yyyyz_xxy, g_x_0_z_0_yyyyz_xy, g_x_0_z_0_yyyyz_xyy, g_x_0_z_0_yyyyz_xyz, g_x_0_z_0_yyyyz_xz, g_x_0_z_0_yyyyz_yy, g_x_0_z_0_yyyyz_yyy, g_x_0_z_0_yyyyz_yyz, g_x_0_z_0_yyyyz_yz, g_x_0_z_0_yyyyz_yzz, g_x_0_z_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyyyz_xx[k] = -g_x_0_z_0_yyyyz_xx[k] * ab_y + g_x_0_z_0_yyyyz_xxy[k];

                g_x_0_z_0_yyyyyz_xy[k] = -g_x_0_z_0_yyyyz_xy[k] * ab_y + g_x_0_z_0_yyyyz_xyy[k];

                g_x_0_z_0_yyyyyz_xz[k] = -g_x_0_z_0_yyyyz_xz[k] * ab_y + g_x_0_z_0_yyyyz_xyz[k];

                g_x_0_z_0_yyyyyz_yy[k] = -g_x_0_z_0_yyyyz_yy[k] * ab_y + g_x_0_z_0_yyyyz_yyy[k];

                g_x_0_z_0_yyyyyz_yz[k] = -g_x_0_z_0_yyyyz_yz[k] * ab_y + g_x_0_z_0_yyyyz_yyz[k];

                g_x_0_z_0_yyyyyz_zz[k] = -g_x_0_z_0_yyyyz_zz[k] * ab_y + g_x_0_z_0_yyyyz_yzz[k];
            }

            /// Set up 474-480 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyyzz_xx = cbuffer.data(id_geom_1010_off + 474 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyzz_xy = cbuffer.data(id_geom_1010_off + 475 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyzz_xz = cbuffer.data(id_geom_1010_off + 476 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyzz_yy = cbuffer.data(id_geom_1010_off + 477 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyzz_yz = cbuffer.data(id_geom_1010_off + 478 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyzz_zz = cbuffer.data(id_geom_1010_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyyzz_xx, g_x_0_z_0_yyyyzz_xy, g_x_0_z_0_yyyyzz_xz, g_x_0_z_0_yyyyzz_yy, g_x_0_z_0_yyyyzz_yz, g_x_0_z_0_yyyyzz_zz, g_x_0_z_0_yyyzz_xx, g_x_0_z_0_yyyzz_xxy, g_x_0_z_0_yyyzz_xy, g_x_0_z_0_yyyzz_xyy, g_x_0_z_0_yyyzz_xyz, g_x_0_z_0_yyyzz_xz, g_x_0_z_0_yyyzz_yy, g_x_0_z_0_yyyzz_yyy, g_x_0_z_0_yyyzz_yyz, g_x_0_z_0_yyyzz_yz, g_x_0_z_0_yyyzz_yzz, g_x_0_z_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyyzz_xx[k] = -g_x_0_z_0_yyyzz_xx[k] * ab_y + g_x_0_z_0_yyyzz_xxy[k];

                g_x_0_z_0_yyyyzz_xy[k] = -g_x_0_z_0_yyyzz_xy[k] * ab_y + g_x_0_z_0_yyyzz_xyy[k];

                g_x_0_z_0_yyyyzz_xz[k] = -g_x_0_z_0_yyyzz_xz[k] * ab_y + g_x_0_z_0_yyyzz_xyz[k];

                g_x_0_z_0_yyyyzz_yy[k] = -g_x_0_z_0_yyyzz_yy[k] * ab_y + g_x_0_z_0_yyyzz_yyy[k];

                g_x_0_z_0_yyyyzz_yz[k] = -g_x_0_z_0_yyyzz_yz[k] * ab_y + g_x_0_z_0_yyyzz_yyz[k];

                g_x_0_z_0_yyyyzz_zz[k] = -g_x_0_z_0_yyyzz_zz[k] * ab_y + g_x_0_z_0_yyyzz_yzz[k];
            }

            /// Set up 480-486 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyzzz_xx = cbuffer.data(id_geom_1010_off + 480 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzzz_xy = cbuffer.data(id_geom_1010_off + 481 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzzz_xz = cbuffer.data(id_geom_1010_off + 482 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzzz_yy = cbuffer.data(id_geom_1010_off + 483 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzzz_yz = cbuffer.data(id_geom_1010_off + 484 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzzz_zz = cbuffer.data(id_geom_1010_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyzzz_xx, g_x_0_z_0_yyyzzz_xy, g_x_0_z_0_yyyzzz_xz, g_x_0_z_0_yyyzzz_yy, g_x_0_z_0_yyyzzz_yz, g_x_0_z_0_yyyzzz_zz, g_x_0_z_0_yyzzz_xx, g_x_0_z_0_yyzzz_xxy, g_x_0_z_0_yyzzz_xy, g_x_0_z_0_yyzzz_xyy, g_x_0_z_0_yyzzz_xyz, g_x_0_z_0_yyzzz_xz, g_x_0_z_0_yyzzz_yy, g_x_0_z_0_yyzzz_yyy, g_x_0_z_0_yyzzz_yyz, g_x_0_z_0_yyzzz_yz, g_x_0_z_0_yyzzz_yzz, g_x_0_z_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyzzz_xx[k] = -g_x_0_z_0_yyzzz_xx[k] * ab_y + g_x_0_z_0_yyzzz_xxy[k];

                g_x_0_z_0_yyyzzz_xy[k] = -g_x_0_z_0_yyzzz_xy[k] * ab_y + g_x_0_z_0_yyzzz_xyy[k];

                g_x_0_z_0_yyyzzz_xz[k] = -g_x_0_z_0_yyzzz_xz[k] * ab_y + g_x_0_z_0_yyzzz_xyz[k];

                g_x_0_z_0_yyyzzz_yy[k] = -g_x_0_z_0_yyzzz_yy[k] * ab_y + g_x_0_z_0_yyzzz_yyy[k];

                g_x_0_z_0_yyyzzz_yz[k] = -g_x_0_z_0_yyzzz_yz[k] * ab_y + g_x_0_z_0_yyzzz_yyz[k];

                g_x_0_z_0_yyyzzz_zz[k] = -g_x_0_z_0_yyzzz_zz[k] * ab_y + g_x_0_z_0_yyzzz_yzz[k];
            }

            /// Set up 486-492 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyzzzz_xx = cbuffer.data(id_geom_1010_off + 486 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzzz_xy = cbuffer.data(id_geom_1010_off + 487 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzzz_xz = cbuffer.data(id_geom_1010_off + 488 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzzz_yy = cbuffer.data(id_geom_1010_off + 489 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzzz_yz = cbuffer.data(id_geom_1010_off + 490 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzzz_zz = cbuffer.data(id_geom_1010_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyzzzz_xx, g_x_0_z_0_yyzzzz_xy, g_x_0_z_0_yyzzzz_xz, g_x_0_z_0_yyzzzz_yy, g_x_0_z_0_yyzzzz_yz, g_x_0_z_0_yyzzzz_zz, g_x_0_z_0_yzzzz_xx, g_x_0_z_0_yzzzz_xxy, g_x_0_z_0_yzzzz_xy, g_x_0_z_0_yzzzz_xyy, g_x_0_z_0_yzzzz_xyz, g_x_0_z_0_yzzzz_xz, g_x_0_z_0_yzzzz_yy, g_x_0_z_0_yzzzz_yyy, g_x_0_z_0_yzzzz_yyz, g_x_0_z_0_yzzzz_yz, g_x_0_z_0_yzzzz_yzz, g_x_0_z_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyzzzz_xx[k] = -g_x_0_z_0_yzzzz_xx[k] * ab_y + g_x_0_z_0_yzzzz_xxy[k];

                g_x_0_z_0_yyzzzz_xy[k] = -g_x_0_z_0_yzzzz_xy[k] * ab_y + g_x_0_z_0_yzzzz_xyy[k];

                g_x_0_z_0_yyzzzz_xz[k] = -g_x_0_z_0_yzzzz_xz[k] * ab_y + g_x_0_z_0_yzzzz_xyz[k];

                g_x_0_z_0_yyzzzz_yy[k] = -g_x_0_z_0_yzzzz_yy[k] * ab_y + g_x_0_z_0_yzzzz_yyy[k];

                g_x_0_z_0_yyzzzz_yz[k] = -g_x_0_z_0_yzzzz_yz[k] * ab_y + g_x_0_z_0_yzzzz_yyz[k];

                g_x_0_z_0_yyzzzz_zz[k] = -g_x_0_z_0_yzzzz_zz[k] * ab_y + g_x_0_z_0_yzzzz_yzz[k];
            }

            /// Set up 492-498 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yzzzzz_xx = cbuffer.data(id_geom_1010_off + 492 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzzz_xy = cbuffer.data(id_geom_1010_off + 493 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzzz_xz = cbuffer.data(id_geom_1010_off + 494 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzzz_yy = cbuffer.data(id_geom_1010_off + 495 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzzz_yz = cbuffer.data(id_geom_1010_off + 496 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzzz_zz = cbuffer.data(id_geom_1010_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yzzzzz_xx, g_x_0_z_0_yzzzzz_xy, g_x_0_z_0_yzzzzz_xz, g_x_0_z_0_yzzzzz_yy, g_x_0_z_0_yzzzzz_yz, g_x_0_z_0_yzzzzz_zz, g_x_0_z_0_zzzzz_xx, g_x_0_z_0_zzzzz_xxy, g_x_0_z_0_zzzzz_xy, g_x_0_z_0_zzzzz_xyy, g_x_0_z_0_zzzzz_xyz, g_x_0_z_0_zzzzz_xz, g_x_0_z_0_zzzzz_yy, g_x_0_z_0_zzzzz_yyy, g_x_0_z_0_zzzzz_yyz, g_x_0_z_0_zzzzz_yz, g_x_0_z_0_zzzzz_yzz, g_x_0_z_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yzzzzz_xx[k] = -g_x_0_z_0_zzzzz_xx[k] * ab_y + g_x_0_z_0_zzzzz_xxy[k];

                g_x_0_z_0_yzzzzz_xy[k] = -g_x_0_z_0_zzzzz_xy[k] * ab_y + g_x_0_z_0_zzzzz_xyy[k];

                g_x_0_z_0_yzzzzz_xz[k] = -g_x_0_z_0_zzzzz_xz[k] * ab_y + g_x_0_z_0_zzzzz_xyz[k];

                g_x_0_z_0_yzzzzz_yy[k] = -g_x_0_z_0_zzzzz_yy[k] * ab_y + g_x_0_z_0_zzzzz_yyy[k];

                g_x_0_z_0_yzzzzz_yz[k] = -g_x_0_z_0_zzzzz_yz[k] * ab_y + g_x_0_z_0_zzzzz_yyz[k];

                g_x_0_z_0_yzzzzz_zz[k] = -g_x_0_z_0_zzzzz_zz[k] * ab_y + g_x_0_z_0_zzzzz_yzz[k];
            }

            /// Set up 498-504 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_zzzzzz_xx = cbuffer.data(id_geom_1010_off + 498 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzzz_xy = cbuffer.data(id_geom_1010_off + 499 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzzz_xz = cbuffer.data(id_geom_1010_off + 500 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzzz_yy = cbuffer.data(id_geom_1010_off + 501 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzzz_yz = cbuffer.data(id_geom_1010_off + 502 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzzz_zz = cbuffer.data(id_geom_1010_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_zzzzz_xx, g_x_0_z_0_zzzzz_xxz, g_x_0_z_0_zzzzz_xy, g_x_0_z_0_zzzzz_xyz, g_x_0_z_0_zzzzz_xz, g_x_0_z_0_zzzzz_xzz, g_x_0_z_0_zzzzz_yy, g_x_0_z_0_zzzzz_yyz, g_x_0_z_0_zzzzz_yz, g_x_0_z_0_zzzzz_yzz, g_x_0_z_0_zzzzz_zz, g_x_0_z_0_zzzzz_zzz, g_x_0_z_0_zzzzzz_xx, g_x_0_z_0_zzzzzz_xy, g_x_0_z_0_zzzzzz_xz, g_x_0_z_0_zzzzzz_yy, g_x_0_z_0_zzzzzz_yz, g_x_0_z_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_zzzzzz_xx[k] = -g_x_0_z_0_zzzzz_xx[k] * ab_z + g_x_0_z_0_zzzzz_xxz[k];

                g_x_0_z_0_zzzzzz_xy[k] = -g_x_0_z_0_zzzzz_xy[k] * ab_z + g_x_0_z_0_zzzzz_xyz[k];

                g_x_0_z_0_zzzzzz_xz[k] = -g_x_0_z_0_zzzzz_xz[k] * ab_z + g_x_0_z_0_zzzzz_xzz[k];

                g_x_0_z_0_zzzzzz_yy[k] = -g_x_0_z_0_zzzzz_yy[k] * ab_z + g_x_0_z_0_zzzzz_yyz[k];

                g_x_0_z_0_zzzzzz_yz[k] = -g_x_0_z_0_zzzzz_yz[k] * ab_z + g_x_0_z_0_zzzzz_yzz[k];

                g_x_0_z_0_zzzzzz_zz[k] = -g_x_0_z_0_zzzzz_zz[k] * ab_z + g_x_0_z_0_zzzzz_zzz[k];
            }

            /// Set up 504-510 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxxx_xx = cbuffer.data(id_geom_1010_off + 504 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxx_xy = cbuffer.data(id_geom_1010_off + 505 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxx_xz = cbuffer.data(id_geom_1010_off + 506 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxx_yy = cbuffer.data(id_geom_1010_off + 507 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxx_yz = cbuffer.data(id_geom_1010_off + 508 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxx_zz = cbuffer.data(id_geom_1010_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxx_xx, g_y_0_x_0_xxxxx_xxx, g_y_0_x_0_xxxxx_xxy, g_y_0_x_0_xxxxx_xxz, g_y_0_x_0_xxxxx_xy, g_y_0_x_0_xxxxx_xyy, g_y_0_x_0_xxxxx_xyz, g_y_0_x_0_xxxxx_xz, g_y_0_x_0_xxxxx_xzz, g_y_0_x_0_xxxxx_yy, g_y_0_x_0_xxxxx_yz, g_y_0_x_0_xxxxx_zz, g_y_0_x_0_xxxxxx_xx, g_y_0_x_0_xxxxxx_xy, g_y_0_x_0_xxxxxx_xz, g_y_0_x_0_xxxxxx_yy, g_y_0_x_0_xxxxxx_yz, g_y_0_x_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxxx_xx[k] = -g_y_0_x_0_xxxxx_xx[k] * ab_x + g_y_0_x_0_xxxxx_xxx[k];

                g_y_0_x_0_xxxxxx_xy[k] = -g_y_0_x_0_xxxxx_xy[k] * ab_x + g_y_0_x_0_xxxxx_xxy[k];

                g_y_0_x_0_xxxxxx_xz[k] = -g_y_0_x_0_xxxxx_xz[k] * ab_x + g_y_0_x_0_xxxxx_xxz[k];

                g_y_0_x_0_xxxxxx_yy[k] = -g_y_0_x_0_xxxxx_yy[k] * ab_x + g_y_0_x_0_xxxxx_xyy[k];

                g_y_0_x_0_xxxxxx_yz[k] = -g_y_0_x_0_xxxxx_yz[k] * ab_x + g_y_0_x_0_xxxxx_xyz[k];

                g_y_0_x_0_xxxxxx_zz[k] = -g_y_0_x_0_xxxxx_zz[k] * ab_x + g_y_0_x_0_xxxxx_xzz[k];
            }

            /// Set up 510-516 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxxy_xx = cbuffer.data(id_geom_1010_off + 510 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxy_xy = cbuffer.data(id_geom_1010_off + 511 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxy_xz = cbuffer.data(id_geom_1010_off + 512 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxy_yy = cbuffer.data(id_geom_1010_off + 513 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxy_yz = cbuffer.data(id_geom_1010_off + 514 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxy_zz = cbuffer.data(id_geom_1010_off + 515 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxxy_xx, g_y_0_x_0_xxxxxy_xy, g_y_0_x_0_xxxxxy_xz, g_y_0_x_0_xxxxxy_yy, g_y_0_x_0_xxxxxy_yz, g_y_0_x_0_xxxxxy_zz, g_y_0_x_0_xxxxy_xx, g_y_0_x_0_xxxxy_xxx, g_y_0_x_0_xxxxy_xxy, g_y_0_x_0_xxxxy_xxz, g_y_0_x_0_xxxxy_xy, g_y_0_x_0_xxxxy_xyy, g_y_0_x_0_xxxxy_xyz, g_y_0_x_0_xxxxy_xz, g_y_0_x_0_xxxxy_xzz, g_y_0_x_0_xxxxy_yy, g_y_0_x_0_xxxxy_yz, g_y_0_x_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxxy_xx[k] = -g_y_0_x_0_xxxxy_xx[k] * ab_x + g_y_0_x_0_xxxxy_xxx[k];

                g_y_0_x_0_xxxxxy_xy[k] = -g_y_0_x_0_xxxxy_xy[k] * ab_x + g_y_0_x_0_xxxxy_xxy[k];

                g_y_0_x_0_xxxxxy_xz[k] = -g_y_0_x_0_xxxxy_xz[k] * ab_x + g_y_0_x_0_xxxxy_xxz[k];

                g_y_0_x_0_xxxxxy_yy[k] = -g_y_0_x_0_xxxxy_yy[k] * ab_x + g_y_0_x_0_xxxxy_xyy[k];

                g_y_0_x_0_xxxxxy_yz[k] = -g_y_0_x_0_xxxxy_yz[k] * ab_x + g_y_0_x_0_xxxxy_xyz[k];

                g_y_0_x_0_xxxxxy_zz[k] = -g_y_0_x_0_xxxxy_zz[k] * ab_x + g_y_0_x_0_xxxxy_xzz[k];
            }

            /// Set up 516-522 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxxz_xx = cbuffer.data(id_geom_1010_off + 516 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxz_xy = cbuffer.data(id_geom_1010_off + 517 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxz_xz = cbuffer.data(id_geom_1010_off + 518 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxz_yy = cbuffer.data(id_geom_1010_off + 519 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxz_yz = cbuffer.data(id_geom_1010_off + 520 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxz_zz = cbuffer.data(id_geom_1010_off + 521 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxxz_xx, g_y_0_x_0_xxxxxz_xy, g_y_0_x_0_xxxxxz_xz, g_y_0_x_0_xxxxxz_yy, g_y_0_x_0_xxxxxz_yz, g_y_0_x_0_xxxxxz_zz, g_y_0_x_0_xxxxz_xx, g_y_0_x_0_xxxxz_xxx, g_y_0_x_0_xxxxz_xxy, g_y_0_x_0_xxxxz_xxz, g_y_0_x_0_xxxxz_xy, g_y_0_x_0_xxxxz_xyy, g_y_0_x_0_xxxxz_xyz, g_y_0_x_0_xxxxz_xz, g_y_0_x_0_xxxxz_xzz, g_y_0_x_0_xxxxz_yy, g_y_0_x_0_xxxxz_yz, g_y_0_x_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxxz_xx[k] = -g_y_0_x_0_xxxxz_xx[k] * ab_x + g_y_0_x_0_xxxxz_xxx[k];

                g_y_0_x_0_xxxxxz_xy[k] = -g_y_0_x_0_xxxxz_xy[k] * ab_x + g_y_0_x_0_xxxxz_xxy[k];

                g_y_0_x_0_xxxxxz_xz[k] = -g_y_0_x_0_xxxxz_xz[k] * ab_x + g_y_0_x_0_xxxxz_xxz[k];

                g_y_0_x_0_xxxxxz_yy[k] = -g_y_0_x_0_xxxxz_yy[k] * ab_x + g_y_0_x_0_xxxxz_xyy[k];

                g_y_0_x_0_xxxxxz_yz[k] = -g_y_0_x_0_xxxxz_yz[k] * ab_x + g_y_0_x_0_xxxxz_xyz[k];

                g_y_0_x_0_xxxxxz_zz[k] = -g_y_0_x_0_xxxxz_zz[k] * ab_x + g_y_0_x_0_xxxxz_xzz[k];
            }

            /// Set up 522-528 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxyy_xx = cbuffer.data(id_geom_1010_off + 522 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyy_xy = cbuffer.data(id_geom_1010_off + 523 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyy_xz = cbuffer.data(id_geom_1010_off + 524 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyy_yy = cbuffer.data(id_geom_1010_off + 525 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyy_yz = cbuffer.data(id_geom_1010_off + 526 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyy_zz = cbuffer.data(id_geom_1010_off + 527 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxyy_xx, g_y_0_x_0_xxxxyy_xy, g_y_0_x_0_xxxxyy_xz, g_y_0_x_0_xxxxyy_yy, g_y_0_x_0_xxxxyy_yz, g_y_0_x_0_xxxxyy_zz, g_y_0_x_0_xxxyy_xx, g_y_0_x_0_xxxyy_xxx, g_y_0_x_0_xxxyy_xxy, g_y_0_x_0_xxxyy_xxz, g_y_0_x_0_xxxyy_xy, g_y_0_x_0_xxxyy_xyy, g_y_0_x_0_xxxyy_xyz, g_y_0_x_0_xxxyy_xz, g_y_0_x_0_xxxyy_xzz, g_y_0_x_0_xxxyy_yy, g_y_0_x_0_xxxyy_yz, g_y_0_x_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxyy_xx[k] = -g_y_0_x_0_xxxyy_xx[k] * ab_x + g_y_0_x_0_xxxyy_xxx[k];

                g_y_0_x_0_xxxxyy_xy[k] = -g_y_0_x_0_xxxyy_xy[k] * ab_x + g_y_0_x_0_xxxyy_xxy[k];

                g_y_0_x_0_xxxxyy_xz[k] = -g_y_0_x_0_xxxyy_xz[k] * ab_x + g_y_0_x_0_xxxyy_xxz[k];

                g_y_0_x_0_xxxxyy_yy[k] = -g_y_0_x_0_xxxyy_yy[k] * ab_x + g_y_0_x_0_xxxyy_xyy[k];

                g_y_0_x_0_xxxxyy_yz[k] = -g_y_0_x_0_xxxyy_yz[k] * ab_x + g_y_0_x_0_xxxyy_xyz[k];

                g_y_0_x_0_xxxxyy_zz[k] = -g_y_0_x_0_xxxyy_zz[k] * ab_x + g_y_0_x_0_xxxyy_xzz[k];
            }

            /// Set up 528-534 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxyz_xx = cbuffer.data(id_geom_1010_off + 528 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyz_xy = cbuffer.data(id_geom_1010_off + 529 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyz_xz = cbuffer.data(id_geom_1010_off + 530 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyz_yy = cbuffer.data(id_geom_1010_off + 531 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyz_yz = cbuffer.data(id_geom_1010_off + 532 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyz_zz = cbuffer.data(id_geom_1010_off + 533 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxyz_xx, g_y_0_x_0_xxxxyz_xy, g_y_0_x_0_xxxxyz_xz, g_y_0_x_0_xxxxyz_yy, g_y_0_x_0_xxxxyz_yz, g_y_0_x_0_xxxxyz_zz, g_y_0_x_0_xxxyz_xx, g_y_0_x_0_xxxyz_xxx, g_y_0_x_0_xxxyz_xxy, g_y_0_x_0_xxxyz_xxz, g_y_0_x_0_xxxyz_xy, g_y_0_x_0_xxxyz_xyy, g_y_0_x_0_xxxyz_xyz, g_y_0_x_0_xxxyz_xz, g_y_0_x_0_xxxyz_xzz, g_y_0_x_0_xxxyz_yy, g_y_0_x_0_xxxyz_yz, g_y_0_x_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxyz_xx[k] = -g_y_0_x_0_xxxyz_xx[k] * ab_x + g_y_0_x_0_xxxyz_xxx[k];

                g_y_0_x_0_xxxxyz_xy[k] = -g_y_0_x_0_xxxyz_xy[k] * ab_x + g_y_0_x_0_xxxyz_xxy[k];

                g_y_0_x_0_xxxxyz_xz[k] = -g_y_0_x_0_xxxyz_xz[k] * ab_x + g_y_0_x_0_xxxyz_xxz[k];

                g_y_0_x_0_xxxxyz_yy[k] = -g_y_0_x_0_xxxyz_yy[k] * ab_x + g_y_0_x_0_xxxyz_xyy[k];

                g_y_0_x_0_xxxxyz_yz[k] = -g_y_0_x_0_xxxyz_yz[k] * ab_x + g_y_0_x_0_xxxyz_xyz[k];

                g_y_0_x_0_xxxxyz_zz[k] = -g_y_0_x_0_xxxyz_zz[k] * ab_x + g_y_0_x_0_xxxyz_xzz[k];
            }

            /// Set up 534-540 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxzz_xx = cbuffer.data(id_geom_1010_off + 534 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxzz_xy = cbuffer.data(id_geom_1010_off + 535 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxzz_xz = cbuffer.data(id_geom_1010_off + 536 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxzz_yy = cbuffer.data(id_geom_1010_off + 537 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxzz_yz = cbuffer.data(id_geom_1010_off + 538 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxzz_zz = cbuffer.data(id_geom_1010_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxzz_xx, g_y_0_x_0_xxxxzz_xy, g_y_0_x_0_xxxxzz_xz, g_y_0_x_0_xxxxzz_yy, g_y_0_x_0_xxxxzz_yz, g_y_0_x_0_xxxxzz_zz, g_y_0_x_0_xxxzz_xx, g_y_0_x_0_xxxzz_xxx, g_y_0_x_0_xxxzz_xxy, g_y_0_x_0_xxxzz_xxz, g_y_0_x_0_xxxzz_xy, g_y_0_x_0_xxxzz_xyy, g_y_0_x_0_xxxzz_xyz, g_y_0_x_0_xxxzz_xz, g_y_0_x_0_xxxzz_xzz, g_y_0_x_0_xxxzz_yy, g_y_0_x_0_xxxzz_yz, g_y_0_x_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxzz_xx[k] = -g_y_0_x_0_xxxzz_xx[k] * ab_x + g_y_0_x_0_xxxzz_xxx[k];

                g_y_0_x_0_xxxxzz_xy[k] = -g_y_0_x_0_xxxzz_xy[k] * ab_x + g_y_0_x_0_xxxzz_xxy[k];

                g_y_0_x_0_xxxxzz_xz[k] = -g_y_0_x_0_xxxzz_xz[k] * ab_x + g_y_0_x_0_xxxzz_xxz[k];

                g_y_0_x_0_xxxxzz_yy[k] = -g_y_0_x_0_xxxzz_yy[k] * ab_x + g_y_0_x_0_xxxzz_xyy[k];

                g_y_0_x_0_xxxxzz_yz[k] = -g_y_0_x_0_xxxzz_yz[k] * ab_x + g_y_0_x_0_xxxzz_xyz[k];

                g_y_0_x_0_xxxxzz_zz[k] = -g_y_0_x_0_xxxzz_zz[k] * ab_x + g_y_0_x_0_xxxzz_xzz[k];
            }

            /// Set up 540-546 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxyyy_xx = cbuffer.data(id_geom_1010_off + 540 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyy_xy = cbuffer.data(id_geom_1010_off + 541 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyy_xz = cbuffer.data(id_geom_1010_off + 542 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyy_yy = cbuffer.data(id_geom_1010_off + 543 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyy_yz = cbuffer.data(id_geom_1010_off + 544 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyy_zz = cbuffer.data(id_geom_1010_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxyyy_xx, g_y_0_x_0_xxxyyy_xy, g_y_0_x_0_xxxyyy_xz, g_y_0_x_0_xxxyyy_yy, g_y_0_x_0_xxxyyy_yz, g_y_0_x_0_xxxyyy_zz, g_y_0_x_0_xxyyy_xx, g_y_0_x_0_xxyyy_xxx, g_y_0_x_0_xxyyy_xxy, g_y_0_x_0_xxyyy_xxz, g_y_0_x_0_xxyyy_xy, g_y_0_x_0_xxyyy_xyy, g_y_0_x_0_xxyyy_xyz, g_y_0_x_0_xxyyy_xz, g_y_0_x_0_xxyyy_xzz, g_y_0_x_0_xxyyy_yy, g_y_0_x_0_xxyyy_yz, g_y_0_x_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxyyy_xx[k] = -g_y_0_x_0_xxyyy_xx[k] * ab_x + g_y_0_x_0_xxyyy_xxx[k];

                g_y_0_x_0_xxxyyy_xy[k] = -g_y_0_x_0_xxyyy_xy[k] * ab_x + g_y_0_x_0_xxyyy_xxy[k];

                g_y_0_x_0_xxxyyy_xz[k] = -g_y_0_x_0_xxyyy_xz[k] * ab_x + g_y_0_x_0_xxyyy_xxz[k];

                g_y_0_x_0_xxxyyy_yy[k] = -g_y_0_x_0_xxyyy_yy[k] * ab_x + g_y_0_x_0_xxyyy_xyy[k];

                g_y_0_x_0_xxxyyy_yz[k] = -g_y_0_x_0_xxyyy_yz[k] * ab_x + g_y_0_x_0_xxyyy_xyz[k];

                g_y_0_x_0_xxxyyy_zz[k] = -g_y_0_x_0_xxyyy_zz[k] * ab_x + g_y_0_x_0_xxyyy_xzz[k];
            }

            /// Set up 546-552 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxyyz_xx = cbuffer.data(id_geom_1010_off + 546 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyz_xy = cbuffer.data(id_geom_1010_off + 547 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyz_xz = cbuffer.data(id_geom_1010_off + 548 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyz_yy = cbuffer.data(id_geom_1010_off + 549 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyz_yz = cbuffer.data(id_geom_1010_off + 550 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyz_zz = cbuffer.data(id_geom_1010_off + 551 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxyyz_xx, g_y_0_x_0_xxxyyz_xy, g_y_0_x_0_xxxyyz_xz, g_y_0_x_0_xxxyyz_yy, g_y_0_x_0_xxxyyz_yz, g_y_0_x_0_xxxyyz_zz, g_y_0_x_0_xxyyz_xx, g_y_0_x_0_xxyyz_xxx, g_y_0_x_0_xxyyz_xxy, g_y_0_x_0_xxyyz_xxz, g_y_0_x_0_xxyyz_xy, g_y_0_x_0_xxyyz_xyy, g_y_0_x_0_xxyyz_xyz, g_y_0_x_0_xxyyz_xz, g_y_0_x_0_xxyyz_xzz, g_y_0_x_0_xxyyz_yy, g_y_0_x_0_xxyyz_yz, g_y_0_x_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxyyz_xx[k] = -g_y_0_x_0_xxyyz_xx[k] * ab_x + g_y_0_x_0_xxyyz_xxx[k];

                g_y_0_x_0_xxxyyz_xy[k] = -g_y_0_x_0_xxyyz_xy[k] * ab_x + g_y_0_x_0_xxyyz_xxy[k];

                g_y_0_x_0_xxxyyz_xz[k] = -g_y_0_x_0_xxyyz_xz[k] * ab_x + g_y_0_x_0_xxyyz_xxz[k];

                g_y_0_x_0_xxxyyz_yy[k] = -g_y_0_x_0_xxyyz_yy[k] * ab_x + g_y_0_x_0_xxyyz_xyy[k];

                g_y_0_x_0_xxxyyz_yz[k] = -g_y_0_x_0_xxyyz_yz[k] * ab_x + g_y_0_x_0_xxyyz_xyz[k];

                g_y_0_x_0_xxxyyz_zz[k] = -g_y_0_x_0_xxyyz_zz[k] * ab_x + g_y_0_x_0_xxyyz_xzz[k];
            }

            /// Set up 552-558 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxyzz_xx = cbuffer.data(id_geom_1010_off + 552 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyzz_xy = cbuffer.data(id_geom_1010_off + 553 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyzz_xz = cbuffer.data(id_geom_1010_off + 554 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyzz_yy = cbuffer.data(id_geom_1010_off + 555 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyzz_yz = cbuffer.data(id_geom_1010_off + 556 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyzz_zz = cbuffer.data(id_geom_1010_off + 557 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxyzz_xx, g_y_0_x_0_xxxyzz_xy, g_y_0_x_0_xxxyzz_xz, g_y_0_x_0_xxxyzz_yy, g_y_0_x_0_xxxyzz_yz, g_y_0_x_0_xxxyzz_zz, g_y_0_x_0_xxyzz_xx, g_y_0_x_0_xxyzz_xxx, g_y_0_x_0_xxyzz_xxy, g_y_0_x_0_xxyzz_xxz, g_y_0_x_0_xxyzz_xy, g_y_0_x_0_xxyzz_xyy, g_y_0_x_0_xxyzz_xyz, g_y_0_x_0_xxyzz_xz, g_y_0_x_0_xxyzz_xzz, g_y_0_x_0_xxyzz_yy, g_y_0_x_0_xxyzz_yz, g_y_0_x_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxyzz_xx[k] = -g_y_0_x_0_xxyzz_xx[k] * ab_x + g_y_0_x_0_xxyzz_xxx[k];

                g_y_0_x_0_xxxyzz_xy[k] = -g_y_0_x_0_xxyzz_xy[k] * ab_x + g_y_0_x_0_xxyzz_xxy[k];

                g_y_0_x_0_xxxyzz_xz[k] = -g_y_0_x_0_xxyzz_xz[k] * ab_x + g_y_0_x_0_xxyzz_xxz[k];

                g_y_0_x_0_xxxyzz_yy[k] = -g_y_0_x_0_xxyzz_yy[k] * ab_x + g_y_0_x_0_xxyzz_xyy[k];

                g_y_0_x_0_xxxyzz_yz[k] = -g_y_0_x_0_xxyzz_yz[k] * ab_x + g_y_0_x_0_xxyzz_xyz[k];

                g_y_0_x_0_xxxyzz_zz[k] = -g_y_0_x_0_xxyzz_zz[k] * ab_x + g_y_0_x_0_xxyzz_xzz[k];
            }

            /// Set up 558-564 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxzzz_xx = cbuffer.data(id_geom_1010_off + 558 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzzz_xy = cbuffer.data(id_geom_1010_off + 559 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzzz_xz = cbuffer.data(id_geom_1010_off + 560 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzzz_yy = cbuffer.data(id_geom_1010_off + 561 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzzz_yz = cbuffer.data(id_geom_1010_off + 562 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzzz_zz = cbuffer.data(id_geom_1010_off + 563 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxzzz_xx, g_y_0_x_0_xxxzzz_xy, g_y_0_x_0_xxxzzz_xz, g_y_0_x_0_xxxzzz_yy, g_y_0_x_0_xxxzzz_yz, g_y_0_x_0_xxxzzz_zz, g_y_0_x_0_xxzzz_xx, g_y_0_x_0_xxzzz_xxx, g_y_0_x_0_xxzzz_xxy, g_y_0_x_0_xxzzz_xxz, g_y_0_x_0_xxzzz_xy, g_y_0_x_0_xxzzz_xyy, g_y_0_x_0_xxzzz_xyz, g_y_0_x_0_xxzzz_xz, g_y_0_x_0_xxzzz_xzz, g_y_0_x_0_xxzzz_yy, g_y_0_x_0_xxzzz_yz, g_y_0_x_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxzzz_xx[k] = -g_y_0_x_0_xxzzz_xx[k] * ab_x + g_y_0_x_0_xxzzz_xxx[k];

                g_y_0_x_0_xxxzzz_xy[k] = -g_y_0_x_0_xxzzz_xy[k] * ab_x + g_y_0_x_0_xxzzz_xxy[k];

                g_y_0_x_0_xxxzzz_xz[k] = -g_y_0_x_0_xxzzz_xz[k] * ab_x + g_y_0_x_0_xxzzz_xxz[k];

                g_y_0_x_0_xxxzzz_yy[k] = -g_y_0_x_0_xxzzz_yy[k] * ab_x + g_y_0_x_0_xxzzz_xyy[k];

                g_y_0_x_0_xxxzzz_yz[k] = -g_y_0_x_0_xxzzz_yz[k] * ab_x + g_y_0_x_0_xxzzz_xyz[k];

                g_y_0_x_0_xxxzzz_zz[k] = -g_y_0_x_0_xxzzz_zz[k] * ab_x + g_y_0_x_0_xxzzz_xzz[k];
            }

            /// Set up 564-570 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyyyy_xx = cbuffer.data(id_geom_1010_off + 564 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyy_xy = cbuffer.data(id_geom_1010_off + 565 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyy_xz = cbuffer.data(id_geom_1010_off + 566 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyy_yy = cbuffer.data(id_geom_1010_off + 567 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyy_yz = cbuffer.data(id_geom_1010_off + 568 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyy_zz = cbuffer.data(id_geom_1010_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyyyy_xx, g_y_0_x_0_xxyyyy_xy, g_y_0_x_0_xxyyyy_xz, g_y_0_x_0_xxyyyy_yy, g_y_0_x_0_xxyyyy_yz, g_y_0_x_0_xxyyyy_zz, g_y_0_x_0_xyyyy_xx, g_y_0_x_0_xyyyy_xxx, g_y_0_x_0_xyyyy_xxy, g_y_0_x_0_xyyyy_xxz, g_y_0_x_0_xyyyy_xy, g_y_0_x_0_xyyyy_xyy, g_y_0_x_0_xyyyy_xyz, g_y_0_x_0_xyyyy_xz, g_y_0_x_0_xyyyy_xzz, g_y_0_x_0_xyyyy_yy, g_y_0_x_0_xyyyy_yz, g_y_0_x_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyyyy_xx[k] = -g_y_0_x_0_xyyyy_xx[k] * ab_x + g_y_0_x_0_xyyyy_xxx[k];

                g_y_0_x_0_xxyyyy_xy[k] = -g_y_0_x_0_xyyyy_xy[k] * ab_x + g_y_0_x_0_xyyyy_xxy[k];

                g_y_0_x_0_xxyyyy_xz[k] = -g_y_0_x_0_xyyyy_xz[k] * ab_x + g_y_0_x_0_xyyyy_xxz[k];

                g_y_0_x_0_xxyyyy_yy[k] = -g_y_0_x_0_xyyyy_yy[k] * ab_x + g_y_0_x_0_xyyyy_xyy[k];

                g_y_0_x_0_xxyyyy_yz[k] = -g_y_0_x_0_xyyyy_yz[k] * ab_x + g_y_0_x_0_xyyyy_xyz[k];

                g_y_0_x_0_xxyyyy_zz[k] = -g_y_0_x_0_xyyyy_zz[k] * ab_x + g_y_0_x_0_xyyyy_xzz[k];
            }

            /// Set up 570-576 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyyyz_xx = cbuffer.data(id_geom_1010_off + 570 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyz_xy = cbuffer.data(id_geom_1010_off + 571 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyz_xz = cbuffer.data(id_geom_1010_off + 572 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyz_yy = cbuffer.data(id_geom_1010_off + 573 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyz_yz = cbuffer.data(id_geom_1010_off + 574 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyz_zz = cbuffer.data(id_geom_1010_off + 575 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyyyz_xx, g_y_0_x_0_xxyyyz_xy, g_y_0_x_0_xxyyyz_xz, g_y_0_x_0_xxyyyz_yy, g_y_0_x_0_xxyyyz_yz, g_y_0_x_0_xxyyyz_zz, g_y_0_x_0_xyyyz_xx, g_y_0_x_0_xyyyz_xxx, g_y_0_x_0_xyyyz_xxy, g_y_0_x_0_xyyyz_xxz, g_y_0_x_0_xyyyz_xy, g_y_0_x_0_xyyyz_xyy, g_y_0_x_0_xyyyz_xyz, g_y_0_x_0_xyyyz_xz, g_y_0_x_0_xyyyz_xzz, g_y_0_x_0_xyyyz_yy, g_y_0_x_0_xyyyz_yz, g_y_0_x_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyyyz_xx[k] = -g_y_0_x_0_xyyyz_xx[k] * ab_x + g_y_0_x_0_xyyyz_xxx[k];

                g_y_0_x_0_xxyyyz_xy[k] = -g_y_0_x_0_xyyyz_xy[k] * ab_x + g_y_0_x_0_xyyyz_xxy[k];

                g_y_0_x_0_xxyyyz_xz[k] = -g_y_0_x_0_xyyyz_xz[k] * ab_x + g_y_0_x_0_xyyyz_xxz[k];

                g_y_0_x_0_xxyyyz_yy[k] = -g_y_0_x_0_xyyyz_yy[k] * ab_x + g_y_0_x_0_xyyyz_xyy[k];

                g_y_0_x_0_xxyyyz_yz[k] = -g_y_0_x_0_xyyyz_yz[k] * ab_x + g_y_0_x_0_xyyyz_xyz[k];

                g_y_0_x_0_xxyyyz_zz[k] = -g_y_0_x_0_xyyyz_zz[k] * ab_x + g_y_0_x_0_xyyyz_xzz[k];
            }

            /// Set up 576-582 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyyzz_xx = cbuffer.data(id_geom_1010_off + 576 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyzz_xy = cbuffer.data(id_geom_1010_off + 577 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyzz_xz = cbuffer.data(id_geom_1010_off + 578 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyzz_yy = cbuffer.data(id_geom_1010_off + 579 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyzz_yz = cbuffer.data(id_geom_1010_off + 580 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyzz_zz = cbuffer.data(id_geom_1010_off + 581 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyyzz_xx, g_y_0_x_0_xxyyzz_xy, g_y_0_x_0_xxyyzz_xz, g_y_0_x_0_xxyyzz_yy, g_y_0_x_0_xxyyzz_yz, g_y_0_x_0_xxyyzz_zz, g_y_0_x_0_xyyzz_xx, g_y_0_x_0_xyyzz_xxx, g_y_0_x_0_xyyzz_xxy, g_y_0_x_0_xyyzz_xxz, g_y_0_x_0_xyyzz_xy, g_y_0_x_0_xyyzz_xyy, g_y_0_x_0_xyyzz_xyz, g_y_0_x_0_xyyzz_xz, g_y_0_x_0_xyyzz_xzz, g_y_0_x_0_xyyzz_yy, g_y_0_x_0_xyyzz_yz, g_y_0_x_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyyzz_xx[k] = -g_y_0_x_0_xyyzz_xx[k] * ab_x + g_y_0_x_0_xyyzz_xxx[k];

                g_y_0_x_0_xxyyzz_xy[k] = -g_y_0_x_0_xyyzz_xy[k] * ab_x + g_y_0_x_0_xyyzz_xxy[k];

                g_y_0_x_0_xxyyzz_xz[k] = -g_y_0_x_0_xyyzz_xz[k] * ab_x + g_y_0_x_0_xyyzz_xxz[k];

                g_y_0_x_0_xxyyzz_yy[k] = -g_y_0_x_0_xyyzz_yy[k] * ab_x + g_y_0_x_0_xyyzz_xyy[k];

                g_y_0_x_0_xxyyzz_yz[k] = -g_y_0_x_0_xyyzz_yz[k] * ab_x + g_y_0_x_0_xyyzz_xyz[k];

                g_y_0_x_0_xxyyzz_zz[k] = -g_y_0_x_0_xyyzz_zz[k] * ab_x + g_y_0_x_0_xyyzz_xzz[k];
            }

            /// Set up 582-588 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyzzz_xx = cbuffer.data(id_geom_1010_off + 582 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzzz_xy = cbuffer.data(id_geom_1010_off + 583 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzzz_xz = cbuffer.data(id_geom_1010_off + 584 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzzz_yy = cbuffer.data(id_geom_1010_off + 585 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzzz_yz = cbuffer.data(id_geom_1010_off + 586 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzzz_zz = cbuffer.data(id_geom_1010_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyzzz_xx, g_y_0_x_0_xxyzzz_xy, g_y_0_x_0_xxyzzz_xz, g_y_0_x_0_xxyzzz_yy, g_y_0_x_0_xxyzzz_yz, g_y_0_x_0_xxyzzz_zz, g_y_0_x_0_xyzzz_xx, g_y_0_x_0_xyzzz_xxx, g_y_0_x_0_xyzzz_xxy, g_y_0_x_0_xyzzz_xxz, g_y_0_x_0_xyzzz_xy, g_y_0_x_0_xyzzz_xyy, g_y_0_x_0_xyzzz_xyz, g_y_0_x_0_xyzzz_xz, g_y_0_x_0_xyzzz_xzz, g_y_0_x_0_xyzzz_yy, g_y_0_x_0_xyzzz_yz, g_y_0_x_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyzzz_xx[k] = -g_y_0_x_0_xyzzz_xx[k] * ab_x + g_y_0_x_0_xyzzz_xxx[k];

                g_y_0_x_0_xxyzzz_xy[k] = -g_y_0_x_0_xyzzz_xy[k] * ab_x + g_y_0_x_0_xyzzz_xxy[k];

                g_y_0_x_0_xxyzzz_xz[k] = -g_y_0_x_0_xyzzz_xz[k] * ab_x + g_y_0_x_0_xyzzz_xxz[k];

                g_y_0_x_0_xxyzzz_yy[k] = -g_y_0_x_0_xyzzz_yy[k] * ab_x + g_y_0_x_0_xyzzz_xyy[k];

                g_y_0_x_0_xxyzzz_yz[k] = -g_y_0_x_0_xyzzz_yz[k] * ab_x + g_y_0_x_0_xyzzz_xyz[k];

                g_y_0_x_0_xxyzzz_zz[k] = -g_y_0_x_0_xyzzz_zz[k] * ab_x + g_y_0_x_0_xyzzz_xzz[k];
            }

            /// Set up 588-594 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxzzzz_xx = cbuffer.data(id_geom_1010_off + 588 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzzz_xy = cbuffer.data(id_geom_1010_off + 589 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzzz_xz = cbuffer.data(id_geom_1010_off + 590 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzzz_yy = cbuffer.data(id_geom_1010_off + 591 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzzz_yz = cbuffer.data(id_geom_1010_off + 592 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzzz_zz = cbuffer.data(id_geom_1010_off + 593 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxzzzz_xx, g_y_0_x_0_xxzzzz_xy, g_y_0_x_0_xxzzzz_xz, g_y_0_x_0_xxzzzz_yy, g_y_0_x_0_xxzzzz_yz, g_y_0_x_0_xxzzzz_zz, g_y_0_x_0_xzzzz_xx, g_y_0_x_0_xzzzz_xxx, g_y_0_x_0_xzzzz_xxy, g_y_0_x_0_xzzzz_xxz, g_y_0_x_0_xzzzz_xy, g_y_0_x_0_xzzzz_xyy, g_y_0_x_0_xzzzz_xyz, g_y_0_x_0_xzzzz_xz, g_y_0_x_0_xzzzz_xzz, g_y_0_x_0_xzzzz_yy, g_y_0_x_0_xzzzz_yz, g_y_0_x_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxzzzz_xx[k] = -g_y_0_x_0_xzzzz_xx[k] * ab_x + g_y_0_x_0_xzzzz_xxx[k];

                g_y_0_x_0_xxzzzz_xy[k] = -g_y_0_x_0_xzzzz_xy[k] * ab_x + g_y_0_x_0_xzzzz_xxy[k];

                g_y_0_x_0_xxzzzz_xz[k] = -g_y_0_x_0_xzzzz_xz[k] * ab_x + g_y_0_x_0_xzzzz_xxz[k];

                g_y_0_x_0_xxzzzz_yy[k] = -g_y_0_x_0_xzzzz_yy[k] * ab_x + g_y_0_x_0_xzzzz_xyy[k];

                g_y_0_x_0_xxzzzz_yz[k] = -g_y_0_x_0_xzzzz_yz[k] * ab_x + g_y_0_x_0_xzzzz_xyz[k];

                g_y_0_x_0_xxzzzz_zz[k] = -g_y_0_x_0_xzzzz_zz[k] * ab_x + g_y_0_x_0_xzzzz_xzz[k];
            }

            /// Set up 594-600 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyyyy_xx = cbuffer.data(id_geom_1010_off + 594 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyy_xy = cbuffer.data(id_geom_1010_off + 595 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyy_xz = cbuffer.data(id_geom_1010_off + 596 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyy_yy = cbuffer.data(id_geom_1010_off + 597 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyy_yz = cbuffer.data(id_geom_1010_off + 598 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyy_zz = cbuffer.data(id_geom_1010_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyyyy_xx, g_y_0_x_0_xyyyyy_xy, g_y_0_x_0_xyyyyy_xz, g_y_0_x_0_xyyyyy_yy, g_y_0_x_0_xyyyyy_yz, g_y_0_x_0_xyyyyy_zz, g_y_0_x_0_yyyyy_xx, g_y_0_x_0_yyyyy_xxx, g_y_0_x_0_yyyyy_xxy, g_y_0_x_0_yyyyy_xxz, g_y_0_x_0_yyyyy_xy, g_y_0_x_0_yyyyy_xyy, g_y_0_x_0_yyyyy_xyz, g_y_0_x_0_yyyyy_xz, g_y_0_x_0_yyyyy_xzz, g_y_0_x_0_yyyyy_yy, g_y_0_x_0_yyyyy_yz, g_y_0_x_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyyyy_xx[k] = -g_y_0_x_0_yyyyy_xx[k] * ab_x + g_y_0_x_0_yyyyy_xxx[k];

                g_y_0_x_0_xyyyyy_xy[k] = -g_y_0_x_0_yyyyy_xy[k] * ab_x + g_y_0_x_0_yyyyy_xxy[k];

                g_y_0_x_0_xyyyyy_xz[k] = -g_y_0_x_0_yyyyy_xz[k] * ab_x + g_y_0_x_0_yyyyy_xxz[k];

                g_y_0_x_0_xyyyyy_yy[k] = -g_y_0_x_0_yyyyy_yy[k] * ab_x + g_y_0_x_0_yyyyy_xyy[k];

                g_y_0_x_0_xyyyyy_yz[k] = -g_y_0_x_0_yyyyy_yz[k] * ab_x + g_y_0_x_0_yyyyy_xyz[k];

                g_y_0_x_0_xyyyyy_zz[k] = -g_y_0_x_0_yyyyy_zz[k] * ab_x + g_y_0_x_0_yyyyy_xzz[k];
            }

            /// Set up 600-606 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyyyz_xx = cbuffer.data(id_geom_1010_off + 600 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyz_xy = cbuffer.data(id_geom_1010_off + 601 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyz_xz = cbuffer.data(id_geom_1010_off + 602 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyz_yy = cbuffer.data(id_geom_1010_off + 603 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyz_yz = cbuffer.data(id_geom_1010_off + 604 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyz_zz = cbuffer.data(id_geom_1010_off + 605 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyyyz_xx, g_y_0_x_0_xyyyyz_xy, g_y_0_x_0_xyyyyz_xz, g_y_0_x_0_xyyyyz_yy, g_y_0_x_0_xyyyyz_yz, g_y_0_x_0_xyyyyz_zz, g_y_0_x_0_yyyyz_xx, g_y_0_x_0_yyyyz_xxx, g_y_0_x_0_yyyyz_xxy, g_y_0_x_0_yyyyz_xxz, g_y_0_x_0_yyyyz_xy, g_y_0_x_0_yyyyz_xyy, g_y_0_x_0_yyyyz_xyz, g_y_0_x_0_yyyyz_xz, g_y_0_x_0_yyyyz_xzz, g_y_0_x_0_yyyyz_yy, g_y_0_x_0_yyyyz_yz, g_y_0_x_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyyyz_xx[k] = -g_y_0_x_0_yyyyz_xx[k] * ab_x + g_y_0_x_0_yyyyz_xxx[k];

                g_y_0_x_0_xyyyyz_xy[k] = -g_y_0_x_0_yyyyz_xy[k] * ab_x + g_y_0_x_0_yyyyz_xxy[k];

                g_y_0_x_0_xyyyyz_xz[k] = -g_y_0_x_0_yyyyz_xz[k] * ab_x + g_y_0_x_0_yyyyz_xxz[k];

                g_y_0_x_0_xyyyyz_yy[k] = -g_y_0_x_0_yyyyz_yy[k] * ab_x + g_y_0_x_0_yyyyz_xyy[k];

                g_y_0_x_0_xyyyyz_yz[k] = -g_y_0_x_0_yyyyz_yz[k] * ab_x + g_y_0_x_0_yyyyz_xyz[k];

                g_y_0_x_0_xyyyyz_zz[k] = -g_y_0_x_0_yyyyz_zz[k] * ab_x + g_y_0_x_0_yyyyz_xzz[k];
            }

            /// Set up 606-612 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyyzz_xx = cbuffer.data(id_geom_1010_off + 606 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyzz_xy = cbuffer.data(id_geom_1010_off + 607 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyzz_xz = cbuffer.data(id_geom_1010_off + 608 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyzz_yy = cbuffer.data(id_geom_1010_off + 609 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyzz_yz = cbuffer.data(id_geom_1010_off + 610 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyzz_zz = cbuffer.data(id_geom_1010_off + 611 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyyzz_xx, g_y_0_x_0_xyyyzz_xy, g_y_0_x_0_xyyyzz_xz, g_y_0_x_0_xyyyzz_yy, g_y_0_x_0_xyyyzz_yz, g_y_0_x_0_xyyyzz_zz, g_y_0_x_0_yyyzz_xx, g_y_0_x_0_yyyzz_xxx, g_y_0_x_0_yyyzz_xxy, g_y_0_x_0_yyyzz_xxz, g_y_0_x_0_yyyzz_xy, g_y_0_x_0_yyyzz_xyy, g_y_0_x_0_yyyzz_xyz, g_y_0_x_0_yyyzz_xz, g_y_0_x_0_yyyzz_xzz, g_y_0_x_0_yyyzz_yy, g_y_0_x_0_yyyzz_yz, g_y_0_x_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyyzz_xx[k] = -g_y_0_x_0_yyyzz_xx[k] * ab_x + g_y_0_x_0_yyyzz_xxx[k];

                g_y_0_x_0_xyyyzz_xy[k] = -g_y_0_x_0_yyyzz_xy[k] * ab_x + g_y_0_x_0_yyyzz_xxy[k];

                g_y_0_x_0_xyyyzz_xz[k] = -g_y_0_x_0_yyyzz_xz[k] * ab_x + g_y_0_x_0_yyyzz_xxz[k];

                g_y_0_x_0_xyyyzz_yy[k] = -g_y_0_x_0_yyyzz_yy[k] * ab_x + g_y_0_x_0_yyyzz_xyy[k];

                g_y_0_x_0_xyyyzz_yz[k] = -g_y_0_x_0_yyyzz_yz[k] * ab_x + g_y_0_x_0_yyyzz_xyz[k];

                g_y_0_x_0_xyyyzz_zz[k] = -g_y_0_x_0_yyyzz_zz[k] * ab_x + g_y_0_x_0_yyyzz_xzz[k];
            }

            /// Set up 612-618 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyzzz_xx = cbuffer.data(id_geom_1010_off + 612 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzzz_xy = cbuffer.data(id_geom_1010_off + 613 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzzz_xz = cbuffer.data(id_geom_1010_off + 614 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzzz_yy = cbuffer.data(id_geom_1010_off + 615 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzzz_yz = cbuffer.data(id_geom_1010_off + 616 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzzz_zz = cbuffer.data(id_geom_1010_off + 617 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyzzz_xx, g_y_0_x_0_xyyzzz_xy, g_y_0_x_0_xyyzzz_xz, g_y_0_x_0_xyyzzz_yy, g_y_0_x_0_xyyzzz_yz, g_y_0_x_0_xyyzzz_zz, g_y_0_x_0_yyzzz_xx, g_y_0_x_0_yyzzz_xxx, g_y_0_x_0_yyzzz_xxy, g_y_0_x_0_yyzzz_xxz, g_y_0_x_0_yyzzz_xy, g_y_0_x_0_yyzzz_xyy, g_y_0_x_0_yyzzz_xyz, g_y_0_x_0_yyzzz_xz, g_y_0_x_0_yyzzz_xzz, g_y_0_x_0_yyzzz_yy, g_y_0_x_0_yyzzz_yz, g_y_0_x_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyzzz_xx[k] = -g_y_0_x_0_yyzzz_xx[k] * ab_x + g_y_0_x_0_yyzzz_xxx[k];

                g_y_0_x_0_xyyzzz_xy[k] = -g_y_0_x_0_yyzzz_xy[k] * ab_x + g_y_0_x_0_yyzzz_xxy[k];

                g_y_0_x_0_xyyzzz_xz[k] = -g_y_0_x_0_yyzzz_xz[k] * ab_x + g_y_0_x_0_yyzzz_xxz[k];

                g_y_0_x_0_xyyzzz_yy[k] = -g_y_0_x_0_yyzzz_yy[k] * ab_x + g_y_0_x_0_yyzzz_xyy[k];

                g_y_0_x_0_xyyzzz_yz[k] = -g_y_0_x_0_yyzzz_yz[k] * ab_x + g_y_0_x_0_yyzzz_xyz[k];

                g_y_0_x_0_xyyzzz_zz[k] = -g_y_0_x_0_yyzzz_zz[k] * ab_x + g_y_0_x_0_yyzzz_xzz[k];
            }

            /// Set up 618-624 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyzzzz_xx = cbuffer.data(id_geom_1010_off + 618 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzzz_xy = cbuffer.data(id_geom_1010_off + 619 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzzz_xz = cbuffer.data(id_geom_1010_off + 620 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzzz_yy = cbuffer.data(id_geom_1010_off + 621 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzzz_yz = cbuffer.data(id_geom_1010_off + 622 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzzz_zz = cbuffer.data(id_geom_1010_off + 623 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyzzzz_xx, g_y_0_x_0_xyzzzz_xy, g_y_0_x_0_xyzzzz_xz, g_y_0_x_0_xyzzzz_yy, g_y_0_x_0_xyzzzz_yz, g_y_0_x_0_xyzzzz_zz, g_y_0_x_0_yzzzz_xx, g_y_0_x_0_yzzzz_xxx, g_y_0_x_0_yzzzz_xxy, g_y_0_x_0_yzzzz_xxz, g_y_0_x_0_yzzzz_xy, g_y_0_x_0_yzzzz_xyy, g_y_0_x_0_yzzzz_xyz, g_y_0_x_0_yzzzz_xz, g_y_0_x_0_yzzzz_xzz, g_y_0_x_0_yzzzz_yy, g_y_0_x_0_yzzzz_yz, g_y_0_x_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyzzzz_xx[k] = -g_y_0_x_0_yzzzz_xx[k] * ab_x + g_y_0_x_0_yzzzz_xxx[k];

                g_y_0_x_0_xyzzzz_xy[k] = -g_y_0_x_0_yzzzz_xy[k] * ab_x + g_y_0_x_0_yzzzz_xxy[k];

                g_y_0_x_0_xyzzzz_xz[k] = -g_y_0_x_0_yzzzz_xz[k] * ab_x + g_y_0_x_0_yzzzz_xxz[k];

                g_y_0_x_0_xyzzzz_yy[k] = -g_y_0_x_0_yzzzz_yy[k] * ab_x + g_y_0_x_0_yzzzz_xyy[k];

                g_y_0_x_0_xyzzzz_yz[k] = -g_y_0_x_0_yzzzz_yz[k] * ab_x + g_y_0_x_0_yzzzz_xyz[k];

                g_y_0_x_0_xyzzzz_zz[k] = -g_y_0_x_0_yzzzz_zz[k] * ab_x + g_y_0_x_0_yzzzz_xzz[k];
            }

            /// Set up 624-630 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xzzzzz_xx = cbuffer.data(id_geom_1010_off + 624 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzzz_xy = cbuffer.data(id_geom_1010_off + 625 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzzz_xz = cbuffer.data(id_geom_1010_off + 626 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzzz_yy = cbuffer.data(id_geom_1010_off + 627 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzzz_yz = cbuffer.data(id_geom_1010_off + 628 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzzz_zz = cbuffer.data(id_geom_1010_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xzzzzz_xx, g_y_0_x_0_xzzzzz_xy, g_y_0_x_0_xzzzzz_xz, g_y_0_x_0_xzzzzz_yy, g_y_0_x_0_xzzzzz_yz, g_y_0_x_0_xzzzzz_zz, g_y_0_x_0_zzzzz_xx, g_y_0_x_0_zzzzz_xxx, g_y_0_x_0_zzzzz_xxy, g_y_0_x_0_zzzzz_xxz, g_y_0_x_0_zzzzz_xy, g_y_0_x_0_zzzzz_xyy, g_y_0_x_0_zzzzz_xyz, g_y_0_x_0_zzzzz_xz, g_y_0_x_0_zzzzz_xzz, g_y_0_x_0_zzzzz_yy, g_y_0_x_0_zzzzz_yz, g_y_0_x_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xzzzzz_xx[k] = -g_y_0_x_0_zzzzz_xx[k] * ab_x + g_y_0_x_0_zzzzz_xxx[k];

                g_y_0_x_0_xzzzzz_xy[k] = -g_y_0_x_0_zzzzz_xy[k] * ab_x + g_y_0_x_0_zzzzz_xxy[k];

                g_y_0_x_0_xzzzzz_xz[k] = -g_y_0_x_0_zzzzz_xz[k] * ab_x + g_y_0_x_0_zzzzz_xxz[k];

                g_y_0_x_0_xzzzzz_yy[k] = -g_y_0_x_0_zzzzz_yy[k] * ab_x + g_y_0_x_0_zzzzz_xyy[k];

                g_y_0_x_0_xzzzzz_yz[k] = -g_y_0_x_0_zzzzz_yz[k] * ab_x + g_y_0_x_0_zzzzz_xyz[k];

                g_y_0_x_0_xzzzzz_zz[k] = -g_y_0_x_0_zzzzz_zz[k] * ab_x + g_y_0_x_0_zzzzz_xzz[k];
            }

            /// Set up 630-636 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyyyy_xx = cbuffer.data(id_geom_1010_off + 630 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyy_xy = cbuffer.data(id_geom_1010_off + 631 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyy_xz = cbuffer.data(id_geom_1010_off + 632 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyy_yy = cbuffer.data(id_geom_1010_off + 633 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyy_yz = cbuffer.data(id_geom_1010_off + 634 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyy_zz = cbuffer.data(id_geom_1010_off + 635 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_yyyyy_xx, g_0_0_x_0_yyyyy_xy, g_0_0_x_0_yyyyy_xz, g_0_0_x_0_yyyyy_yy, g_0_0_x_0_yyyyy_yz, g_0_0_x_0_yyyyy_zz, g_y_0_x_0_yyyyy_xx, g_y_0_x_0_yyyyy_xxy, g_y_0_x_0_yyyyy_xy, g_y_0_x_0_yyyyy_xyy, g_y_0_x_0_yyyyy_xyz, g_y_0_x_0_yyyyy_xz, g_y_0_x_0_yyyyy_yy, g_y_0_x_0_yyyyy_yyy, g_y_0_x_0_yyyyy_yyz, g_y_0_x_0_yyyyy_yz, g_y_0_x_0_yyyyy_yzz, g_y_0_x_0_yyyyy_zz, g_y_0_x_0_yyyyyy_xx, g_y_0_x_0_yyyyyy_xy, g_y_0_x_0_yyyyyy_xz, g_y_0_x_0_yyyyyy_yy, g_y_0_x_0_yyyyyy_yz, g_y_0_x_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyyyy_xx[k] = -g_0_0_x_0_yyyyy_xx[k] - g_y_0_x_0_yyyyy_xx[k] * ab_y + g_y_0_x_0_yyyyy_xxy[k];

                g_y_0_x_0_yyyyyy_xy[k] = -g_0_0_x_0_yyyyy_xy[k] - g_y_0_x_0_yyyyy_xy[k] * ab_y + g_y_0_x_0_yyyyy_xyy[k];

                g_y_0_x_0_yyyyyy_xz[k] = -g_0_0_x_0_yyyyy_xz[k] - g_y_0_x_0_yyyyy_xz[k] * ab_y + g_y_0_x_0_yyyyy_xyz[k];

                g_y_0_x_0_yyyyyy_yy[k] = -g_0_0_x_0_yyyyy_yy[k] - g_y_0_x_0_yyyyy_yy[k] * ab_y + g_y_0_x_0_yyyyy_yyy[k];

                g_y_0_x_0_yyyyyy_yz[k] = -g_0_0_x_0_yyyyy_yz[k] - g_y_0_x_0_yyyyy_yz[k] * ab_y + g_y_0_x_0_yyyyy_yyz[k];

                g_y_0_x_0_yyyyyy_zz[k] = -g_0_0_x_0_yyyyy_zz[k] - g_y_0_x_0_yyyyy_zz[k] * ab_y + g_y_0_x_0_yyyyy_yzz[k];
            }

            /// Set up 636-642 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyyyz_xx = cbuffer.data(id_geom_1010_off + 636 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyz_xy = cbuffer.data(id_geom_1010_off + 637 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyz_xz = cbuffer.data(id_geom_1010_off + 638 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyz_yy = cbuffer.data(id_geom_1010_off + 639 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyz_yz = cbuffer.data(id_geom_1010_off + 640 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyz_zz = cbuffer.data(id_geom_1010_off + 641 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyyyy_xx, g_y_0_x_0_yyyyy_xxz, g_y_0_x_0_yyyyy_xy, g_y_0_x_0_yyyyy_xyz, g_y_0_x_0_yyyyy_xz, g_y_0_x_0_yyyyy_xzz, g_y_0_x_0_yyyyy_yy, g_y_0_x_0_yyyyy_yyz, g_y_0_x_0_yyyyy_yz, g_y_0_x_0_yyyyy_yzz, g_y_0_x_0_yyyyy_zz, g_y_0_x_0_yyyyy_zzz, g_y_0_x_0_yyyyyz_xx, g_y_0_x_0_yyyyyz_xy, g_y_0_x_0_yyyyyz_xz, g_y_0_x_0_yyyyyz_yy, g_y_0_x_0_yyyyyz_yz, g_y_0_x_0_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyyyz_xx[k] = -g_y_0_x_0_yyyyy_xx[k] * ab_z + g_y_0_x_0_yyyyy_xxz[k];

                g_y_0_x_0_yyyyyz_xy[k] = -g_y_0_x_0_yyyyy_xy[k] * ab_z + g_y_0_x_0_yyyyy_xyz[k];

                g_y_0_x_0_yyyyyz_xz[k] = -g_y_0_x_0_yyyyy_xz[k] * ab_z + g_y_0_x_0_yyyyy_xzz[k];

                g_y_0_x_0_yyyyyz_yy[k] = -g_y_0_x_0_yyyyy_yy[k] * ab_z + g_y_0_x_0_yyyyy_yyz[k];

                g_y_0_x_0_yyyyyz_yz[k] = -g_y_0_x_0_yyyyy_yz[k] * ab_z + g_y_0_x_0_yyyyy_yzz[k];

                g_y_0_x_0_yyyyyz_zz[k] = -g_y_0_x_0_yyyyy_zz[k] * ab_z + g_y_0_x_0_yyyyy_zzz[k];
            }

            /// Set up 642-648 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyyzz_xx = cbuffer.data(id_geom_1010_off + 642 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyzz_xy = cbuffer.data(id_geom_1010_off + 643 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyzz_xz = cbuffer.data(id_geom_1010_off + 644 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyzz_yy = cbuffer.data(id_geom_1010_off + 645 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyzz_yz = cbuffer.data(id_geom_1010_off + 646 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyzz_zz = cbuffer.data(id_geom_1010_off + 647 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyyyz_xx, g_y_0_x_0_yyyyz_xxz, g_y_0_x_0_yyyyz_xy, g_y_0_x_0_yyyyz_xyz, g_y_0_x_0_yyyyz_xz, g_y_0_x_0_yyyyz_xzz, g_y_0_x_0_yyyyz_yy, g_y_0_x_0_yyyyz_yyz, g_y_0_x_0_yyyyz_yz, g_y_0_x_0_yyyyz_yzz, g_y_0_x_0_yyyyz_zz, g_y_0_x_0_yyyyz_zzz, g_y_0_x_0_yyyyzz_xx, g_y_0_x_0_yyyyzz_xy, g_y_0_x_0_yyyyzz_xz, g_y_0_x_0_yyyyzz_yy, g_y_0_x_0_yyyyzz_yz, g_y_0_x_0_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyyzz_xx[k] = -g_y_0_x_0_yyyyz_xx[k] * ab_z + g_y_0_x_0_yyyyz_xxz[k];

                g_y_0_x_0_yyyyzz_xy[k] = -g_y_0_x_0_yyyyz_xy[k] * ab_z + g_y_0_x_0_yyyyz_xyz[k];

                g_y_0_x_0_yyyyzz_xz[k] = -g_y_0_x_0_yyyyz_xz[k] * ab_z + g_y_0_x_0_yyyyz_xzz[k];

                g_y_0_x_0_yyyyzz_yy[k] = -g_y_0_x_0_yyyyz_yy[k] * ab_z + g_y_0_x_0_yyyyz_yyz[k];

                g_y_0_x_0_yyyyzz_yz[k] = -g_y_0_x_0_yyyyz_yz[k] * ab_z + g_y_0_x_0_yyyyz_yzz[k];

                g_y_0_x_0_yyyyzz_zz[k] = -g_y_0_x_0_yyyyz_zz[k] * ab_z + g_y_0_x_0_yyyyz_zzz[k];
            }

            /// Set up 648-654 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyzzz_xx = cbuffer.data(id_geom_1010_off + 648 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzzz_xy = cbuffer.data(id_geom_1010_off + 649 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzzz_xz = cbuffer.data(id_geom_1010_off + 650 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzzz_yy = cbuffer.data(id_geom_1010_off + 651 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzzz_yz = cbuffer.data(id_geom_1010_off + 652 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzzz_zz = cbuffer.data(id_geom_1010_off + 653 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyyzz_xx, g_y_0_x_0_yyyzz_xxz, g_y_0_x_0_yyyzz_xy, g_y_0_x_0_yyyzz_xyz, g_y_0_x_0_yyyzz_xz, g_y_0_x_0_yyyzz_xzz, g_y_0_x_0_yyyzz_yy, g_y_0_x_0_yyyzz_yyz, g_y_0_x_0_yyyzz_yz, g_y_0_x_0_yyyzz_yzz, g_y_0_x_0_yyyzz_zz, g_y_0_x_0_yyyzz_zzz, g_y_0_x_0_yyyzzz_xx, g_y_0_x_0_yyyzzz_xy, g_y_0_x_0_yyyzzz_xz, g_y_0_x_0_yyyzzz_yy, g_y_0_x_0_yyyzzz_yz, g_y_0_x_0_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyzzz_xx[k] = -g_y_0_x_0_yyyzz_xx[k] * ab_z + g_y_0_x_0_yyyzz_xxz[k];

                g_y_0_x_0_yyyzzz_xy[k] = -g_y_0_x_0_yyyzz_xy[k] * ab_z + g_y_0_x_0_yyyzz_xyz[k];

                g_y_0_x_0_yyyzzz_xz[k] = -g_y_0_x_0_yyyzz_xz[k] * ab_z + g_y_0_x_0_yyyzz_xzz[k];

                g_y_0_x_0_yyyzzz_yy[k] = -g_y_0_x_0_yyyzz_yy[k] * ab_z + g_y_0_x_0_yyyzz_yyz[k];

                g_y_0_x_0_yyyzzz_yz[k] = -g_y_0_x_0_yyyzz_yz[k] * ab_z + g_y_0_x_0_yyyzz_yzz[k];

                g_y_0_x_0_yyyzzz_zz[k] = -g_y_0_x_0_yyyzz_zz[k] * ab_z + g_y_0_x_0_yyyzz_zzz[k];
            }

            /// Set up 654-660 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyzzzz_xx = cbuffer.data(id_geom_1010_off + 654 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzzz_xy = cbuffer.data(id_geom_1010_off + 655 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzzz_xz = cbuffer.data(id_geom_1010_off + 656 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzzz_yy = cbuffer.data(id_geom_1010_off + 657 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzzz_yz = cbuffer.data(id_geom_1010_off + 658 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzzz_zz = cbuffer.data(id_geom_1010_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyzzz_xx, g_y_0_x_0_yyzzz_xxz, g_y_0_x_0_yyzzz_xy, g_y_0_x_0_yyzzz_xyz, g_y_0_x_0_yyzzz_xz, g_y_0_x_0_yyzzz_xzz, g_y_0_x_0_yyzzz_yy, g_y_0_x_0_yyzzz_yyz, g_y_0_x_0_yyzzz_yz, g_y_0_x_0_yyzzz_yzz, g_y_0_x_0_yyzzz_zz, g_y_0_x_0_yyzzz_zzz, g_y_0_x_0_yyzzzz_xx, g_y_0_x_0_yyzzzz_xy, g_y_0_x_0_yyzzzz_xz, g_y_0_x_0_yyzzzz_yy, g_y_0_x_0_yyzzzz_yz, g_y_0_x_0_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyzzzz_xx[k] = -g_y_0_x_0_yyzzz_xx[k] * ab_z + g_y_0_x_0_yyzzz_xxz[k];

                g_y_0_x_0_yyzzzz_xy[k] = -g_y_0_x_0_yyzzz_xy[k] * ab_z + g_y_0_x_0_yyzzz_xyz[k];

                g_y_0_x_0_yyzzzz_xz[k] = -g_y_0_x_0_yyzzz_xz[k] * ab_z + g_y_0_x_0_yyzzz_xzz[k];

                g_y_0_x_0_yyzzzz_yy[k] = -g_y_0_x_0_yyzzz_yy[k] * ab_z + g_y_0_x_0_yyzzz_yyz[k];

                g_y_0_x_0_yyzzzz_yz[k] = -g_y_0_x_0_yyzzz_yz[k] * ab_z + g_y_0_x_0_yyzzz_yzz[k];

                g_y_0_x_0_yyzzzz_zz[k] = -g_y_0_x_0_yyzzz_zz[k] * ab_z + g_y_0_x_0_yyzzz_zzz[k];
            }

            /// Set up 660-666 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yzzzzz_xx = cbuffer.data(id_geom_1010_off + 660 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzzz_xy = cbuffer.data(id_geom_1010_off + 661 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzzz_xz = cbuffer.data(id_geom_1010_off + 662 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzzz_yy = cbuffer.data(id_geom_1010_off + 663 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzzz_yz = cbuffer.data(id_geom_1010_off + 664 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzzz_zz = cbuffer.data(id_geom_1010_off + 665 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yzzzz_xx, g_y_0_x_0_yzzzz_xxz, g_y_0_x_0_yzzzz_xy, g_y_0_x_0_yzzzz_xyz, g_y_0_x_0_yzzzz_xz, g_y_0_x_0_yzzzz_xzz, g_y_0_x_0_yzzzz_yy, g_y_0_x_0_yzzzz_yyz, g_y_0_x_0_yzzzz_yz, g_y_0_x_0_yzzzz_yzz, g_y_0_x_0_yzzzz_zz, g_y_0_x_0_yzzzz_zzz, g_y_0_x_0_yzzzzz_xx, g_y_0_x_0_yzzzzz_xy, g_y_0_x_0_yzzzzz_xz, g_y_0_x_0_yzzzzz_yy, g_y_0_x_0_yzzzzz_yz, g_y_0_x_0_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yzzzzz_xx[k] = -g_y_0_x_0_yzzzz_xx[k] * ab_z + g_y_0_x_0_yzzzz_xxz[k];

                g_y_0_x_0_yzzzzz_xy[k] = -g_y_0_x_0_yzzzz_xy[k] * ab_z + g_y_0_x_0_yzzzz_xyz[k];

                g_y_0_x_0_yzzzzz_xz[k] = -g_y_0_x_0_yzzzz_xz[k] * ab_z + g_y_0_x_0_yzzzz_xzz[k];

                g_y_0_x_0_yzzzzz_yy[k] = -g_y_0_x_0_yzzzz_yy[k] * ab_z + g_y_0_x_0_yzzzz_yyz[k];

                g_y_0_x_0_yzzzzz_yz[k] = -g_y_0_x_0_yzzzz_yz[k] * ab_z + g_y_0_x_0_yzzzz_yzz[k];

                g_y_0_x_0_yzzzzz_zz[k] = -g_y_0_x_0_yzzzz_zz[k] * ab_z + g_y_0_x_0_yzzzz_zzz[k];
            }

            /// Set up 666-672 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_zzzzzz_xx = cbuffer.data(id_geom_1010_off + 666 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzzz_xy = cbuffer.data(id_geom_1010_off + 667 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzzz_xz = cbuffer.data(id_geom_1010_off + 668 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzzz_yy = cbuffer.data(id_geom_1010_off + 669 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzzz_yz = cbuffer.data(id_geom_1010_off + 670 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzzz_zz = cbuffer.data(id_geom_1010_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_zzzzz_xx, g_y_0_x_0_zzzzz_xxz, g_y_0_x_0_zzzzz_xy, g_y_0_x_0_zzzzz_xyz, g_y_0_x_0_zzzzz_xz, g_y_0_x_0_zzzzz_xzz, g_y_0_x_0_zzzzz_yy, g_y_0_x_0_zzzzz_yyz, g_y_0_x_0_zzzzz_yz, g_y_0_x_0_zzzzz_yzz, g_y_0_x_0_zzzzz_zz, g_y_0_x_0_zzzzz_zzz, g_y_0_x_0_zzzzzz_xx, g_y_0_x_0_zzzzzz_xy, g_y_0_x_0_zzzzzz_xz, g_y_0_x_0_zzzzzz_yy, g_y_0_x_0_zzzzzz_yz, g_y_0_x_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_zzzzzz_xx[k] = -g_y_0_x_0_zzzzz_xx[k] * ab_z + g_y_0_x_0_zzzzz_xxz[k];

                g_y_0_x_0_zzzzzz_xy[k] = -g_y_0_x_0_zzzzz_xy[k] * ab_z + g_y_0_x_0_zzzzz_xyz[k];

                g_y_0_x_0_zzzzzz_xz[k] = -g_y_0_x_0_zzzzz_xz[k] * ab_z + g_y_0_x_0_zzzzz_xzz[k];

                g_y_0_x_0_zzzzzz_yy[k] = -g_y_0_x_0_zzzzz_yy[k] * ab_z + g_y_0_x_0_zzzzz_yyz[k];

                g_y_0_x_0_zzzzzz_yz[k] = -g_y_0_x_0_zzzzz_yz[k] * ab_z + g_y_0_x_0_zzzzz_yzz[k];

                g_y_0_x_0_zzzzzz_zz[k] = -g_y_0_x_0_zzzzz_zz[k] * ab_z + g_y_0_x_0_zzzzz_zzz[k];
            }

            /// Set up 672-678 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxxx_xx = cbuffer.data(id_geom_1010_off + 672 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxx_xy = cbuffer.data(id_geom_1010_off + 673 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxx_xz = cbuffer.data(id_geom_1010_off + 674 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxx_yy = cbuffer.data(id_geom_1010_off + 675 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxx_yz = cbuffer.data(id_geom_1010_off + 676 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxx_zz = cbuffer.data(id_geom_1010_off + 677 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxx_xx, g_y_0_y_0_xxxxx_xxx, g_y_0_y_0_xxxxx_xxy, g_y_0_y_0_xxxxx_xxz, g_y_0_y_0_xxxxx_xy, g_y_0_y_0_xxxxx_xyy, g_y_0_y_0_xxxxx_xyz, g_y_0_y_0_xxxxx_xz, g_y_0_y_0_xxxxx_xzz, g_y_0_y_0_xxxxx_yy, g_y_0_y_0_xxxxx_yz, g_y_0_y_0_xxxxx_zz, g_y_0_y_0_xxxxxx_xx, g_y_0_y_0_xxxxxx_xy, g_y_0_y_0_xxxxxx_xz, g_y_0_y_0_xxxxxx_yy, g_y_0_y_0_xxxxxx_yz, g_y_0_y_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxxx_xx[k] = -g_y_0_y_0_xxxxx_xx[k] * ab_x + g_y_0_y_0_xxxxx_xxx[k];

                g_y_0_y_0_xxxxxx_xy[k] = -g_y_0_y_0_xxxxx_xy[k] * ab_x + g_y_0_y_0_xxxxx_xxy[k];

                g_y_0_y_0_xxxxxx_xz[k] = -g_y_0_y_0_xxxxx_xz[k] * ab_x + g_y_0_y_0_xxxxx_xxz[k];

                g_y_0_y_0_xxxxxx_yy[k] = -g_y_0_y_0_xxxxx_yy[k] * ab_x + g_y_0_y_0_xxxxx_xyy[k];

                g_y_0_y_0_xxxxxx_yz[k] = -g_y_0_y_0_xxxxx_yz[k] * ab_x + g_y_0_y_0_xxxxx_xyz[k];

                g_y_0_y_0_xxxxxx_zz[k] = -g_y_0_y_0_xxxxx_zz[k] * ab_x + g_y_0_y_0_xxxxx_xzz[k];
            }

            /// Set up 678-684 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxxy_xx = cbuffer.data(id_geom_1010_off + 678 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxy_xy = cbuffer.data(id_geom_1010_off + 679 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxy_xz = cbuffer.data(id_geom_1010_off + 680 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxy_yy = cbuffer.data(id_geom_1010_off + 681 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxy_yz = cbuffer.data(id_geom_1010_off + 682 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxy_zz = cbuffer.data(id_geom_1010_off + 683 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxxy_xx, g_y_0_y_0_xxxxxy_xy, g_y_0_y_0_xxxxxy_xz, g_y_0_y_0_xxxxxy_yy, g_y_0_y_0_xxxxxy_yz, g_y_0_y_0_xxxxxy_zz, g_y_0_y_0_xxxxy_xx, g_y_0_y_0_xxxxy_xxx, g_y_0_y_0_xxxxy_xxy, g_y_0_y_0_xxxxy_xxz, g_y_0_y_0_xxxxy_xy, g_y_0_y_0_xxxxy_xyy, g_y_0_y_0_xxxxy_xyz, g_y_0_y_0_xxxxy_xz, g_y_0_y_0_xxxxy_xzz, g_y_0_y_0_xxxxy_yy, g_y_0_y_0_xxxxy_yz, g_y_0_y_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxxy_xx[k] = -g_y_0_y_0_xxxxy_xx[k] * ab_x + g_y_0_y_0_xxxxy_xxx[k];

                g_y_0_y_0_xxxxxy_xy[k] = -g_y_0_y_0_xxxxy_xy[k] * ab_x + g_y_0_y_0_xxxxy_xxy[k];

                g_y_0_y_0_xxxxxy_xz[k] = -g_y_0_y_0_xxxxy_xz[k] * ab_x + g_y_0_y_0_xxxxy_xxz[k];

                g_y_0_y_0_xxxxxy_yy[k] = -g_y_0_y_0_xxxxy_yy[k] * ab_x + g_y_0_y_0_xxxxy_xyy[k];

                g_y_0_y_0_xxxxxy_yz[k] = -g_y_0_y_0_xxxxy_yz[k] * ab_x + g_y_0_y_0_xxxxy_xyz[k];

                g_y_0_y_0_xxxxxy_zz[k] = -g_y_0_y_0_xxxxy_zz[k] * ab_x + g_y_0_y_0_xxxxy_xzz[k];
            }

            /// Set up 684-690 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxxz_xx = cbuffer.data(id_geom_1010_off + 684 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxz_xy = cbuffer.data(id_geom_1010_off + 685 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxz_xz = cbuffer.data(id_geom_1010_off + 686 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxz_yy = cbuffer.data(id_geom_1010_off + 687 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxz_yz = cbuffer.data(id_geom_1010_off + 688 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxz_zz = cbuffer.data(id_geom_1010_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxxz_xx, g_y_0_y_0_xxxxxz_xy, g_y_0_y_0_xxxxxz_xz, g_y_0_y_0_xxxxxz_yy, g_y_0_y_0_xxxxxz_yz, g_y_0_y_0_xxxxxz_zz, g_y_0_y_0_xxxxz_xx, g_y_0_y_0_xxxxz_xxx, g_y_0_y_0_xxxxz_xxy, g_y_0_y_0_xxxxz_xxz, g_y_0_y_0_xxxxz_xy, g_y_0_y_0_xxxxz_xyy, g_y_0_y_0_xxxxz_xyz, g_y_0_y_0_xxxxz_xz, g_y_0_y_0_xxxxz_xzz, g_y_0_y_0_xxxxz_yy, g_y_0_y_0_xxxxz_yz, g_y_0_y_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxxz_xx[k] = -g_y_0_y_0_xxxxz_xx[k] * ab_x + g_y_0_y_0_xxxxz_xxx[k];

                g_y_0_y_0_xxxxxz_xy[k] = -g_y_0_y_0_xxxxz_xy[k] * ab_x + g_y_0_y_0_xxxxz_xxy[k];

                g_y_0_y_0_xxxxxz_xz[k] = -g_y_0_y_0_xxxxz_xz[k] * ab_x + g_y_0_y_0_xxxxz_xxz[k];

                g_y_0_y_0_xxxxxz_yy[k] = -g_y_0_y_0_xxxxz_yy[k] * ab_x + g_y_0_y_0_xxxxz_xyy[k];

                g_y_0_y_0_xxxxxz_yz[k] = -g_y_0_y_0_xxxxz_yz[k] * ab_x + g_y_0_y_0_xxxxz_xyz[k];

                g_y_0_y_0_xxxxxz_zz[k] = -g_y_0_y_0_xxxxz_zz[k] * ab_x + g_y_0_y_0_xxxxz_xzz[k];
            }

            /// Set up 690-696 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxyy_xx = cbuffer.data(id_geom_1010_off + 690 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyy_xy = cbuffer.data(id_geom_1010_off + 691 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyy_xz = cbuffer.data(id_geom_1010_off + 692 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyy_yy = cbuffer.data(id_geom_1010_off + 693 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyy_yz = cbuffer.data(id_geom_1010_off + 694 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyy_zz = cbuffer.data(id_geom_1010_off + 695 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxyy_xx, g_y_0_y_0_xxxxyy_xy, g_y_0_y_0_xxxxyy_xz, g_y_0_y_0_xxxxyy_yy, g_y_0_y_0_xxxxyy_yz, g_y_0_y_0_xxxxyy_zz, g_y_0_y_0_xxxyy_xx, g_y_0_y_0_xxxyy_xxx, g_y_0_y_0_xxxyy_xxy, g_y_0_y_0_xxxyy_xxz, g_y_0_y_0_xxxyy_xy, g_y_0_y_0_xxxyy_xyy, g_y_0_y_0_xxxyy_xyz, g_y_0_y_0_xxxyy_xz, g_y_0_y_0_xxxyy_xzz, g_y_0_y_0_xxxyy_yy, g_y_0_y_0_xxxyy_yz, g_y_0_y_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxyy_xx[k] = -g_y_0_y_0_xxxyy_xx[k] * ab_x + g_y_0_y_0_xxxyy_xxx[k];

                g_y_0_y_0_xxxxyy_xy[k] = -g_y_0_y_0_xxxyy_xy[k] * ab_x + g_y_0_y_0_xxxyy_xxy[k];

                g_y_0_y_0_xxxxyy_xz[k] = -g_y_0_y_0_xxxyy_xz[k] * ab_x + g_y_0_y_0_xxxyy_xxz[k];

                g_y_0_y_0_xxxxyy_yy[k] = -g_y_0_y_0_xxxyy_yy[k] * ab_x + g_y_0_y_0_xxxyy_xyy[k];

                g_y_0_y_0_xxxxyy_yz[k] = -g_y_0_y_0_xxxyy_yz[k] * ab_x + g_y_0_y_0_xxxyy_xyz[k];

                g_y_0_y_0_xxxxyy_zz[k] = -g_y_0_y_0_xxxyy_zz[k] * ab_x + g_y_0_y_0_xxxyy_xzz[k];
            }

            /// Set up 696-702 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxyz_xx = cbuffer.data(id_geom_1010_off + 696 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyz_xy = cbuffer.data(id_geom_1010_off + 697 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyz_xz = cbuffer.data(id_geom_1010_off + 698 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyz_yy = cbuffer.data(id_geom_1010_off + 699 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyz_yz = cbuffer.data(id_geom_1010_off + 700 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyz_zz = cbuffer.data(id_geom_1010_off + 701 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxyz_xx, g_y_0_y_0_xxxxyz_xy, g_y_0_y_0_xxxxyz_xz, g_y_0_y_0_xxxxyz_yy, g_y_0_y_0_xxxxyz_yz, g_y_0_y_0_xxxxyz_zz, g_y_0_y_0_xxxyz_xx, g_y_0_y_0_xxxyz_xxx, g_y_0_y_0_xxxyz_xxy, g_y_0_y_0_xxxyz_xxz, g_y_0_y_0_xxxyz_xy, g_y_0_y_0_xxxyz_xyy, g_y_0_y_0_xxxyz_xyz, g_y_0_y_0_xxxyz_xz, g_y_0_y_0_xxxyz_xzz, g_y_0_y_0_xxxyz_yy, g_y_0_y_0_xxxyz_yz, g_y_0_y_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxyz_xx[k] = -g_y_0_y_0_xxxyz_xx[k] * ab_x + g_y_0_y_0_xxxyz_xxx[k];

                g_y_0_y_0_xxxxyz_xy[k] = -g_y_0_y_0_xxxyz_xy[k] * ab_x + g_y_0_y_0_xxxyz_xxy[k];

                g_y_0_y_0_xxxxyz_xz[k] = -g_y_0_y_0_xxxyz_xz[k] * ab_x + g_y_0_y_0_xxxyz_xxz[k];

                g_y_0_y_0_xxxxyz_yy[k] = -g_y_0_y_0_xxxyz_yy[k] * ab_x + g_y_0_y_0_xxxyz_xyy[k];

                g_y_0_y_0_xxxxyz_yz[k] = -g_y_0_y_0_xxxyz_yz[k] * ab_x + g_y_0_y_0_xxxyz_xyz[k];

                g_y_0_y_0_xxxxyz_zz[k] = -g_y_0_y_0_xxxyz_zz[k] * ab_x + g_y_0_y_0_xxxyz_xzz[k];
            }

            /// Set up 702-708 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxzz_xx = cbuffer.data(id_geom_1010_off + 702 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxzz_xy = cbuffer.data(id_geom_1010_off + 703 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxzz_xz = cbuffer.data(id_geom_1010_off + 704 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxzz_yy = cbuffer.data(id_geom_1010_off + 705 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxzz_yz = cbuffer.data(id_geom_1010_off + 706 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxzz_zz = cbuffer.data(id_geom_1010_off + 707 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxzz_xx, g_y_0_y_0_xxxxzz_xy, g_y_0_y_0_xxxxzz_xz, g_y_0_y_0_xxxxzz_yy, g_y_0_y_0_xxxxzz_yz, g_y_0_y_0_xxxxzz_zz, g_y_0_y_0_xxxzz_xx, g_y_0_y_0_xxxzz_xxx, g_y_0_y_0_xxxzz_xxy, g_y_0_y_0_xxxzz_xxz, g_y_0_y_0_xxxzz_xy, g_y_0_y_0_xxxzz_xyy, g_y_0_y_0_xxxzz_xyz, g_y_0_y_0_xxxzz_xz, g_y_0_y_0_xxxzz_xzz, g_y_0_y_0_xxxzz_yy, g_y_0_y_0_xxxzz_yz, g_y_0_y_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxzz_xx[k] = -g_y_0_y_0_xxxzz_xx[k] * ab_x + g_y_0_y_0_xxxzz_xxx[k];

                g_y_0_y_0_xxxxzz_xy[k] = -g_y_0_y_0_xxxzz_xy[k] * ab_x + g_y_0_y_0_xxxzz_xxy[k];

                g_y_0_y_0_xxxxzz_xz[k] = -g_y_0_y_0_xxxzz_xz[k] * ab_x + g_y_0_y_0_xxxzz_xxz[k];

                g_y_0_y_0_xxxxzz_yy[k] = -g_y_0_y_0_xxxzz_yy[k] * ab_x + g_y_0_y_0_xxxzz_xyy[k];

                g_y_0_y_0_xxxxzz_yz[k] = -g_y_0_y_0_xxxzz_yz[k] * ab_x + g_y_0_y_0_xxxzz_xyz[k];

                g_y_0_y_0_xxxxzz_zz[k] = -g_y_0_y_0_xxxzz_zz[k] * ab_x + g_y_0_y_0_xxxzz_xzz[k];
            }

            /// Set up 708-714 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxyyy_xx = cbuffer.data(id_geom_1010_off + 708 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyy_xy = cbuffer.data(id_geom_1010_off + 709 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyy_xz = cbuffer.data(id_geom_1010_off + 710 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyy_yy = cbuffer.data(id_geom_1010_off + 711 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyy_yz = cbuffer.data(id_geom_1010_off + 712 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyy_zz = cbuffer.data(id_geom_1010_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxyyy_xx, g_y_0_y_0_xxxyyy_xy, g_y_0_y_0_xxxyyy_xz, g_y_0_y_0_xxxyyy_yy, g_y_0_y_0_xxxyyy_yz, g_y_0_y_0_xxxyyy_zz, g_y_0_y_0_xxyyy_xx, g_y_0_y_0_xxyyy_xxx, g_y_0_y_0_xxyyy_xxy, g_y_0_y_0_xxyyy_xxz, g_y_0_y_0_xxyyy_xy, g_y_0_y_0_xxyyy_xyy, g_y_0_y_0_xxyyy_xyz, g_y_0_y_0_xxyyy_xz, g_y_0_y_0_xxyyy_xzz, g_y_0_y_0_xxyyy_yy, g_y_0_y_0_xxyyy_yz, g_y_0_y_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxyyy_xx[k] = -g_y_0_y_0_xxyyy_xx[k] * ab_x + g_y_0_y_0_xxyyy_xxx[k];

                g_y_0_y_0_xxxyyy_xy[k] = -g_y_0_y_0_xxyyy_xy[k] * ab_x + g_y_0_y_0_xxyyy_xxy[k];

                g_y_0_y_0_xxxyyy_xz[k] = -g_y_0_y_0_xxyyy_xz[k] * ab_x + g_y_0_y_0_xxyyy_xxz[k];

                g_y_0_y_0_xxxyyy_yy[k] = -g_y_0_y_0_xxyyy_yy[k] * ab_x + g_y_0_y_0_xxyyy_xyy[k];

                g_y_0_y_0_xxxyyy_yz[k] = -g_y_0_y_0_xxyyy_yz[k] * ab_x + g_y_0_y_0_xxyyy_xyz[k];

                g_y_0_y_0_xxxyyy_zz[k] = -g_y_0_y_0_xxyyy_zz[k] * ab_x + g_y_0_y_0_xxyyy_xzz[k];
            }

            /// Set up 714-720 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxyyz_xx = cbuffer.data(id_geom_1010_off + 714 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyz_xy = cbuffer.data(id_geom_1010_off + 715 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyz_xz = cbuffer.data(id_geom_1010_off + 716 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyz_yy = cbuffer.data(id_geom_1010_off + 717 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyz_yz = cbuffer.data(id_geom_1010_off + 718 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyz_zz = cbuffer.data(id_geom_1010_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxyyz_xx, g_y_0_y_0_xxxyyz_xy, g_y_0_y_0_xxxyyz_xz, g_y_0_y_0_xxxyyz_yy, g_y_0_y_0_xxxyyz_yz, g_y_0_y_0_xxxyyz_zz, g_y_0_y_0_xxyyz_xx, g_y_0_y_0_xxyyz_xxx, g_y_0_y_0_xxyyz_xxy, g_y_0_y_0_xxyyz_xxz, g_y_0_y_0_xxyyz_xy, g_y_0_y_0_xxyyz_xyy, g_y_0_y_0_xxyyz_xyz, g_y_0_y_0_xxyyz_xz, g_y_0_y_0_xxyyz_xzz, g_y_0_y_0_xxyyz_yy, g_y_0_y_0_xxyyz_yz, g_y_0_y_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxyyz_xx[k] = -g_y_0_y_0_xxyyz_xx[k] * ab_x + g_y_0_y_0_xxyyz_xxx[k];

                g_y_0_y_0_xxxyyz_xy[k] = -g_y_0_y_0_xxyyz_xy[k] * ab_x + g_y_0_y_0_xxyyz_xxy[k];

                g_y_0_y_0_xxxyyz_xz[k] = -g_y_0_y_0_xxyyz_xz[k] * ab_x + g_y_0_y_0_xxyyz_xxz[k];

                g_y_0_y_0_xxxyyz_yy[k] = -g_y_0_y_0_xxyyz_yy[k] * ab_x + g_y_0_y_0_xxyyz_xyy[k];

                g_y_0_y_0_xxxyyz_yz[k] = -g_y_0_y_0_xxyyz_yz[k] * ab_x + g_y_0_y_0_xxyyz_xyz[k];

                g_y_0_y_0_xxxyyz_zz[k] = -g_y_0_y_0_xxyyz_zz[k] * ab_x + g_y_0_y_0_xxyyz_xzz[k];
            }

            /// Set up 720-726 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxyzz_xx = cbuffer.data(id_geom_1010_off + 720 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyzz_xy = cbuffer.data(id_geom_1010_off + 721 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyzz_xz = cbuffer.data(id_geom_1010_off + 722 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyzz_yy = cbuffer.data(id_geom_1010_off + 723 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyzz_yz = cbuffer.data(id_geom_1010_off + 724 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyzz_zz = cbuffer.data(id_geom_1010_off + 725 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxyzz_xx, g_y_0_y_0_xxxyzz_xy, g_y_0_y_0_xxxyzz_xz, g_y_0_y_0_xxxyzz_yy, g_y_0_y_0_xxxyzz_yz, g_y_0_y_0_xxxyzz_zz, g_y_0_y_0_xxyzz_xx, g_y_0_y_0_xxyzz_xxx, g_y_0_y_0_xxyzz_xxy, g_y_0_y_0_xxyzz_xxz, g_y_0_y_0_xxyzz_xy, g_y_0_y_0_xxyzz_xyy, g_y_0_y_0_xxyzz_xyz, g_y_0_y_0_xxyzz_xz, g_y_0_y_0_xxyzz_xzz, g_y_0_y_0_xxyzz_yy, g_y_0_y_0_xxyzz_yz, g_y_0_y_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxyzz_xx[k] = -g_y_0_y_0_xxyzz_xx[k] * ab_x + g_y_0_y_0_xxyzz_xxx[k];

                g_y_0_y_0_xxxyzz_xy[k] = -g_y_0_y_0_xxyzz_xy[k] * ab_x + g_y_0_y_0_xxyzz_xxy[k];

                g_y_0_y_0_xxxyzz_xz[k] = -g_y_0_y_0_xxyzz_xz[k] * ab_x + g_y_0_y_0_xxyzz_xxz[k];

                g_y_0_y_0_xxxyzz_yy[k] = -g_y_0_y_0_xxyzz_yy[k] * ab_x + g_y_0_y_0_xxyzz_xyy[k];

                g_y_0_y_0_xxxyzz_yz[k] = -g_y_0_y_0_xxyzz_yz[k] * ab_x + g_y_0_y_0_xxyzz_xyz[k];

                g_y_0_y_0_xxxyzz_zz[k] = -g_y_0_y_0_xxyzz_zz[k] * ab_x + g_y_0_y_0_xxyzz_xzz[k];
            }

            /// Set up 726-732 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxzzz_xx = cbuffer.data(id_geom_1010_off + 726 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzzz_xy = cbuffer.data(id_geom_1010_off + 727 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzzz_xz = cbuffer.data(id_geom_1010_off + 728 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzzz_yy = cbuffer.data(id_geom_1010_off + 729 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzzz_yz = cbuffer.data(id_geom_1010_off + 730 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzzz_zz = cbuffer.data(id_geom_1010_off + 731 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxzzz_xx, g_y_0_y_0_xxxzzz_xy, g_y_0_y_0_xxxzzz_xz, g_y_0_y_0_xxxzzz_yy, g_y_0_y_0_xxxzzz_yz, g_y_0_y_0_xxxzzz_zz, g_y_0_y_0_xxzzz_xx, g_y_0_y_0_xxzzz_xxx, g_y_0_y_0_xxzzz_xxy, g_y_0_y_0_xxzzz_xxz, g_y_0_y_0_xxzzz_xy, g_y_0_y_0_xxzzz_xyy, g_y_0_y_0_xxzzz_xyz, g_y_0_y_0_xxzzz_xz, g_y_0_y_0_xxzzz_xzz, g_y_0_y_0_xxzzz_yy, g_y_0_y_0_xxzzz_yz, g_y_0_y_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxzzz_xx[k] = -g_y_0_y_0_xxzzz_xx[k] * ab_x + g_y_0_y_0_xxzzz_xxx[k];

                g_y_0_y_0_xxxzzz_xy[k] = -g_y_0_y_0_xxzzz_xy[k] * ab_x + g_y_0_y_0_xxzzz_xxy[k];

                g_y_0_y_0_xxxzzz_xz[k] = -g_y_0_y_0_xxzzz_xz[k] * ab_x + g_y_0_y_0_xxzzz_xxz[k];

                g_y_0_y_0_xxxzzz_yy[k] = -g_y_0_y_0_xxzzz_yy[k] * ab_x + g_y_0_y_0_xxzzz_xyy[k];

                g_y_0_y_0_xxxzzz_yz[k] = -g_y_0_y_0_xxzzz_yz[k] * ab_x + g_y_0_y_0_xxzzz_xyz[k];

                g_y_0_y_0_xxxzzz_zz[k] = -g_y_0_y_0_xxzzz_zz[k] * ab_x + g_y_0_y_0_xxzzz_xzz[k];
            }

            /// Set up 732-738 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyyyy_xx = cbuffer.data(id_geom_1010_off + 732 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyy_xy = cbuffer.data(id_geom_1010_off + 733 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyy_xz = cbuffer.data(id_geom_1010_off + 734 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyy_yy = cbuffer.data(id_geom_1010_off + 735 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyy_yz = cbuffer.data(id_geom_1010_off + 736 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyy_zz = cbuffer.data(id_geom_1010_off + 737 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyyyy_xx, g_y_0_y_0_xxyyyy_xy, g_y_0_y_0_xxyyyy_xz, g_y_0_y_0_xxyyyy_yy, g_y_0_y_0_xxyyyy_yz, g_y_0_y_0_xxyyyy_zz, g_y_0_y_0_xyyyy_xx, g_y_0_y_0_xyyyy_xxx, g_y_0_y_0_xyyyy_xxy, g_y_0_y_0_xyyyy_xxz, g_y_0_y_0_xyyyy_xy, g_y_0_y_0_xyyyy_xyy, g_y_0_y_0_xyyyy_xyz, g_y_0_y_0_xyyyy_xz, g_y_0_y_0_xyyyy_xzz, g_y_0_y_0_xyyyy_yy, g_y_0_y_0_xyyyy_yz, g_y_0_y_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyyyy_xx[k] = -g_y_0_y_0_xyyyy_xx[k] * ab_x + g_y_0_y_0_xyyyy_xxx[k];

                g_y_0_y_0_xxyyyy_xy[k] = -g_y_0_y_0_xyyyy_xy[k] * ab_x + g_y_0_y_0_xyyyy_xxy[k];

                g_y_0_y_0_xxyyyy_xz[k] = -g_y_0_y_0_xyyyy_xz[k] * ab_x + g_y_0_y_0_xyyyy_xxz[k];

                g_y_0_y_0_xxyyyy_yy[k] = -g_y_0_y_0_xyyyy_yy[k] * ab_x + g_y_0_y_0_xyyyy_xyy[k];

                g_y_0_y_0_xxyyyy_yz[k] = -g_y_0_y_0_xyyyy_yz[k] * ab_x + g_y_0_y_0_xyyyy_xyz[k];

                g_y_0_y_0_xxyyyy_zz[k] = -g_y_0_y_0_xyyyy_zz[k] * ab_x + g_y_0_y_0_xyyyy_xzz[k];
            }

            /// Set up 738-744 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyyyz_xx = cbuffer.data(id_geom_1010_off + 738 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyz_xy = cbuffer.data(id_geom_1010_off + 739 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyz_xz = cbuffer.data(id_geom_1010_off + 740 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyz_yy = cbuffer.data(id_geom_1010_off + 741 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyz_yz = cbuffer.data(id_geom_1010_off + 742 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyz_zz = cbuffer.data(id_geom_1010_off + 743 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyyyz_xx, g_y_0_y_0_xxyyyz_xy, g_y_0_y_0_xxyyyz_xz, g_y_0_y_0_xxyyyz_yy, g_y_0_y_0_xxyyyz_yz, g_y_0_y_0_xxyyyz_zz, g_y_0_y_0_xyyyz_xx, g_y_0_y_0_xyyyz_xxx, g_y_0_y_0_xyyyz_xxy, g_y_0_y_0_xyyyz_xxz, g_y_0_y_0_xyyyz_xy, g_y_0_y_0_xyyyz_xyy, g_y_0_y_0_xyyyz_xyz, g_y_0_y_0_xyyyz_xz, g_y_0_y_0_xyyyz_xzz, g_y_0_y_0_xyyyz_yy, g_y_0_y_0_xyyyz_yz, g_y_0_y_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyyyz_xx[k] = -g_y_0_y_0_xyyyz_xx[k] * ab_x + g_y_0_y_0_xyyyz_xxx[k];

                g_y_0_y_0_xxyyyz_xy[k] = -g_y_0_y_0_xyyyz_xy[k] * ab_x + g_y_0_y_0_xyyyz_xxy[k];

                g_y_0_y_0_xxyyyz_xz[k] = -g_y_0_y_0_xyyyz_xz[k] * ab_x + g_y_0_y_0_xyyyz_xxz[k];

                g_y_0_y_0_xxyyyz_yy[k] = -g_y_0_y_0_xyyyz_yy[k] * ab_x + g_y_0_y_0_xyyyz_xyy[k];

                g_y_0_y_0_xxyyyz_yz[k] = -g_y_0_y_0_xyyyz_yz[k] * ab_x + g_y_0_y_0_xyyyz_xyz[k];

                g_y_0_y_0_xxyyyz_zz[k] = -g_y_0_y_0_xyyyz_zz[k] * ab_x + g_y_0_y_0_xyyyz_xzz[k];
            }

            /// Set up 744-750 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyyzz_xx = cbuffer.data(id_geom_1010_off + 744 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyzz_xy = cbuffer.data(id_geom_1010_off + 745 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyzz_xz = cbuffer.data(id_geom_1010_off + 746 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyzz_yy = cbuffer.data(id_geom_1010_off + 747 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyzz_yz = cbuffer.data(id_geom_1010_off + 748 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyzz_zz = cbuffer.data(id_geom_1010_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyyzz_xx, g_y_0_y_0_xxyyzz_xy, g_y_0_y_0_xxyyzz_xz, g_y_0_y_0_xxyyzz_yy, g_y_0_y_0_xxyyzz_yz, g_y_0_y_0_xxyyzz_zz, g_y_0_y_0_xyyzz_xx, g_y_0_y_0_xyyzz_xxx, g_y_0_y_0_xyyzz_xxy, g_y_0_y_0_xyyzz_xxz, g_y_0_y_0_xyyzz_xy, g_y_0_y_0_xyyzz_xyy, g_y_0_y_0_xyyzz_xyz, g_y_0_y_0_xyyzz_xz, g_y_0_y_0_xyyzz_xzz, g_y_0_y_0_xyyzz_yy, g_y_0_y_0_xyyzz_yz, g_y_0_y_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyyzz_xx[k] = -g_y_0_y_0_xyyzz_xx[k] * ab_x + g_y_0_y_0_xyyzz_xxx[k];

                g_y_0_y_0_xxyyzz_xy[k] = -g_y_0_y_0_xyyzz_xy[k] * ab_x + g_y_0_y_0_xyyzz_xxy[k];

                g_y_0_y_0_xxyyzz_xz[k] = -g_y_0_y_0_xyyzz_xz[k] * ab_x + g_y_0_y_0_xyyzz_xxz[k];

                g_y_0_y_0_xxyyzz_yy[k] = -g_y_0_y_0_xyyzz_yy[k] * ab_x + g_y_0_y_0_xyyzz_xyy[k];

                g_y_0_y_0_xxyyzz_yz[k] = -g_y_0_y_0_xyyzz_yz[k] * ab_x + g_y_0_y_0_xyyzz_xyz[k];

                g_y_0_y_0_xxyyzz_zz[k] = -g_y_0_y_0_xyyzz_zz[k] * ab_x + g_y_0_y_0_xyyzz_xzz[k];
            }

            /// Set up 750-756 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyzzz_xx = cbuffer.data(id_geom_1010_off + 750 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzzz_xy = cbuffer.data(id_geom_1010_off + 751 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzzz_xz = cbuffer.data(id_geom_1010_off + 752 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzzz_yy = cbuffer.data(id_geom_1010_off + 753 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzzz_yz = cbuffer.data(id_geom_1010_off + 754 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzzz_zz = cbuffer.data(id_geom_1010_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyzzz_xx, g_y_0_y_0_xxyzzz_xy, g_y_0_y_0_xxyzzz_xz, g_y_0_y_0_xxyzzz_yy, g_y_0_y_0_xxyzzz_yz, g_y_0_y_0_xxyzzz_zz, g_y_0_y_0_xyzzz_xx, g_y_0_y_0_xyzzz_xxx, g_y_0_y_0_xyzzz_xxy, g_y_0_y_0_xyzzz_xxz, g_y_0_y_0_xyzzz_xy, g_y_0_y_0_xyzzz_xyy, g_y_0_y_0_xyzzz_xyz, g_y_0_y_0_xyzzz_xz, g_y_0_y_0_xyzzz_xzz, g_y_0_y_0_xyzzz_yy, g_y_0_y_0_xyzzz_yz, g_y_0_y_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyzzz_xx[k] = -g_y_0_y_0_xyzzz_xx[k] * ab_x + g_y_0_y_0_xyzzz_xxx[k];

                g_y_0_y_0_xxyzzz_xy[k] = -g_y_0_y_0_xyzzz_xy[k] * ab_x + g_y_0_y_0_xyzzz_xxy[k];

                g_y_0_y_0_xxyzzz_xz[k] = -g_y_0_y_0_xyzzz_xz[k] * ab_x + g_y_0_y_0_xyzzz_xxz[k];

                g_y_0_y_0_xxyzzz_yy[k] = -g_y_0_y_0_xyzzz_yy[k] * ab_x + g_y_0_y_0_xyzzz_xyy[k];

                g_y_0_y_0_xxyzzz_yz[k] = -g_y_0_y_0_xyzzz_yz[k] * ab_x + g_y_0_y_0_xyzzz_xyz[k];

                g_y_0_y_0_xxyzzz_zz[k] = -g_y_0_y_0_xyzzz_zz[k] * ab_x + g_y_0_y_0_xyzzz_xzz[k];
            }

            /// Set up 756-762 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxzzzz_xx = cbuffer.data(id_geom_1010_off + 756 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzzz_xy = cbuffer.data(id_geom_1010_off + 757 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzzz_xz = cbuffer.data(id_geom_1010_off + 758 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzzz_yy = cbuffer.data(id_geom_1010_off + 759 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzzz_yz = cbuffer.data(id_geom_1010_off + 760 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzzz_zz = cbuffer.data(id_geom_1010_off + 761 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxzzzz_xx, g_y_0_y_0_xxzzzz_xy, g_y_0_y_0_xxzzzz_xz, g_y_0_y_0_xxzzzz_yy, g_y_0_y_0_xxzzzz_yz, g_y_0_y_0_xxzzzz_zz, g_y_0_y_0_xzzzz_xx, g_y_0_y_0_xzzzz_xxx, g_y_0_y_0_xzzzz_xxy, g_y_0_y_0_xzzzz_xxz, g_y_0_y_0_xzzzz_xy, g_y_0_y_0_xzzzz_xyy, g_y_0_y_0_xzzzz_xyz, g_y_0_y_0_xzzzz_xz, g_y_0_y_0_xzzzz_xzz, g_y_0_y_0_xzzzz_yy, g_y_0_y_0_xzzzz_yz, g_y_0_y_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxzzzz_xx[k] = -g_y_0_y_0_xzzzz_xx[k] * ab_x + g_y_0_y_0_xzzzz_xxx[k];

                g_y_0_y_0_xxzzzz_xy[k] = -g_y_0_y_0_xzzzz_xy[k] * ab_x + g_y_0_y_0_xzzzz_xxy[k];

                g_y_0_y_0_xxzzzz_xz[k] = -g_y_0_y_0_xzzzz_xz[k] * ab_x + g_y_0_y_0_xzzzz_xxz[k];

                g_y_0_y_0_xxzzzz_yy[k] = -g_y_0_y_0_xzzzz_yy[k] * ab_x + g_y_0_y_0_xzzzz_xyy[k];

                g_y_0_y_0_xxzzzz_yz[k] = -g_y_0_y_0_xzzzz_yz[k] * ab_x + g_y_0_y_0_xzzzz_xyz[k];

                g_y_0_y_0_xxzzzz_zz[k] = -g_y_0_y_0_xzzzz_zz[k] * ab_x + g_y_0_y_0_xzzzz_xzz[k];
            }

            /// Set up 762-768 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyyyy_xx = cbuffer.data(id_geom_1010_off + 762 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyy_xy = cbuffer.data(id_geom_1010_off + 763 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyy_xz = cbuffer.data(id_geom_1010_off + 764 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyy_yy = cbuffer.data(id_geom_1010_off + 765 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyy_yz = cbuffer.data(id_geom_1010_off + 766 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyy_zz = cbuffer.data(id_geom_1010_off + 767 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyyyy_xx, g_y_0_y_0_xyyyyy_xy, g_y_0_y_0_xyyyyy_xz, g_y_0_y_0_xyyyyy_yy, g_y_0_y_0_xyyyyy_yz, g_y_0_y_0_xyyyyy_zz, g_y_0_y_0_yyyyy_xx, g_y_0_y_0_yyyyy_xxx, g_y_0_y_0_yyyyy_xxy, g_y_0_y_0_yyyyy_xxz, g_y_0_y_0_yyyyy_xy, g_y_0_y_0_yyyyy_xyy, g_y_0_y_0_yyyyy_xyz, g_y_0_y_0_yyyyy_xz, g_y_0_y_0_yyyyy_xzz, g_y_0_y_0_yyyyy_yy, g_y_0_y_0_yyyyy_yz, g_y_0_y_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyyyy_xx[k] = -g_y_0_y_0_yyyyy_xx[k] * ab_x + g_y_0_y_0_yyyyy_xxx[k];

                g_y_0_y_0_xyyyyy_xy[k] = -g_y_0_y_0_yyyyy_xy[k] * ab_x + g_y_0_y_0_yyyyy_xxy[k];

                g_y_0_y_0_xyyyyy_xz[k] = -g_y_0_y_0_yyyyy_xz[k] * ab_x + g_y_0_y_0_yyyyy_xxz[k];

                g_y_0_y_0_xyyyyy_yy[k] = -g_y_0_y_0_yyyyy_yy[k] * ab_x + g_y_0_y_0_yyyyy_xyy[k];

                g_y_0_y_0_xyyyyy_yz[k] = -g_y_0_y_0_yyyyy_yz[k] * ab_x + g_y_0_y_0_yyyyy_xyz[k];

                g_y_0_y_0_xyyyyy_zz[k] = -g_y_0_y_0_yyyyy_zz[k] * ab_x + g_y_0_y_0_yyyyy_xzz[k];
            }

            /// Set up 768-774 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyyyz_xx = cbuffer.data(id_geom_1010_off + 768 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyz_xy = cbuffer.data(id_geom_1010_off + 769 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyz_xz = cbuffer.data(id_geom_1010_off + 770 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyz_yy = cbuffer.data(id_geom_1010_off + 771 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyz_yz = cbuffer.data(id_geom_1010_off + 772 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyz_zz = cbuffer.data(id_geom_1010_off + 773 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyyyz_xx, g_y_0_y_0_xyyyyz_xy, g_y_0_y_0_xyyyyz_xz, g_y_0_y_0_xyyyyz_yy, g_y_0_y_0_xyyyyz_yz, g_y_0_y_0_xyyyyz_zz, g_y_0_y_0_yyyyz_xx, g_y_0_y_0_yyyyz_xxx, g_y_0_y_0_yyyyz_xxy, g_y_0_y_0_yyyyz_xxz, g_y_0_y_0_yyyyz_xy, g_y_0_y_0_yyyyz_xyy, g_y_0_y_0_yyyyz_xyz, g_y_0_y_0_yyyyz_xz, g_y_0_y_0_yyyyz_xzz, g_y_0_y_0_yyyyz_yy, g_y_0_y_0_yyyyz_yz, g_y_0_y_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyyyz_xx[k] = -g_y_0_y_0_yyyyz_xx[k] * ab_x + g_y_0_y_0_yyyyz_xxx[k];

                g_y_0_y_0_xyyyyz_xy[k] = -g_y_0_y_0_yyyyz_xy[k] * ab_x + g_y_0_y_0_yyyyz_xxy[k];

                g_y_0_y_0_xyyyyz_xz[k] = -g_y_0_y_0_yyyyz_xz[k] * ab_x + g_y_0_y_0_yyyyz_xxz[k];

                g_y_0_y_0_xyyyyz_yy[k] = -g_y_0_y_0_yyyyz_yy[k] * ab_x + g_y_0_y_0_yyyyz_xyy[k];

                g_y_0_y_0_xyyyyz_yz[k] = -g_y_0_y_0_yyyyz_yz[k] * ab_x + g_y_0_y_0_yyyyz_xyz[k];

                g_y_0_y_0_xyyyyz_zz[k] = -g_y_0_y_0_yyyyz_zz[k] * ab_x + g_y_0_y_0_yyyyz_xzz[k];
            }

            /// Set up 774-780 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyyzz_xx = cbuffer.data(id_geom_1010_off + 774 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyzz_xy = cbuffer.data(id_geom_1010_off + 775 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyzz_xz = cbuffer.data(id_geom_1010_off + 776 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyzz_yy = cbuffer.data(id_geom_1010_off + 777 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyzz_yz = cbuffer.data(id_geom_1010_off + 778 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyzz_zz = cbuffer.data(id_geom_1010_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyyzz_xx, g_y_0_y_0_xyyyzz_xy, g_y_0_y_0_xyyyzz_xz, g_y_0_y_0_xyyyzz_yy, g_y_0_y_0_xyyyzz_yz, g_y_0_y_0_xyyyzz_zz, g_y_0_y_0_yyyzz_xx, g_y_0_y_0_yyyzz_xxx, g_y_0_y_0_yyyzz_xxy, g_y_0_y_0_yyyzz_xxz, g_y_0_y_0_yyyzz_xy, g_y_0_y_0_yyyzz_xyy, g_y_0_y_0_yyyzz_xyz, g_y_0_y_0_yyyzz_xz, g_y_0_y_0_yyyzz_xzz, g_y_0_y_0_yyyzz_yy, g_y_0_y_0_yyyzz_yz, g_y_0_y_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyyzz_xx[k] = -g_y_0_y_0_yyyzz_xx[k] * ab_x + g_y_0_y_0_yyyzz_xxx[k];

                g_y_0_y_0_xyyyzz_xy[k] = -g_y_0_y_0_yyyzz_xy[k] * ab_x + g_y_0_y_0_yyyzz_xxy[k];

                g_y_0_y_0_xyyyzz_xz[k] = -g_y_0_y_0_yyyzz_xz[k] * ab_x + g_y_0_y_0_yyyzz_xxz[k];

                g_y_0_y_0_xyyyzz_yy[k] = -g_y_0_y_0_yyyzz_yy[k] * ab_x + g_y_0_y_0_yyyzz_xyy[k];

                g_y_0_y_0_xyyyzz_yz[k] = -g_y_0_y_0_yyyzz_yz[k] * ab_x + g_y_0_y_0_yyyzz_xyz[k];

                g_y_0_y_0_xyyyzz_zz[k] = -g_y_0_y_0_yyyzz_zz[k] * ab_x + g_y_0_y_0_yyyzz_xzz[k];
            }

            /// Set up 780-786 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyzzz_xx = cbuffer.data(id_geom_1010_off + 780 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzzz_xy = cbuffer.data(id_geom_1010_off + 781 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzzz_xz = cbuffer.data(id_geom_1010_off + 782 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzzz_yy = cbuffer.data(id_geom_1010_off + 783 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzzz_yz = cbuffer.data(id_geom_1010_off + 784 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzzz_zz = cbuffer.data(id_geom_1010_off + 785 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyzzz_xx, g_y_0_y_0_xyyzzz_xy, g_y_0_y_0_xyyzzz_xz, g_y_0_y_0_xyyzzz_yy, g_y_0_y_0_xyyzzz_yz, g_y_0_y_0_xyyzzz_zz, g_y_0_y_0_yyzzz_xx, g_y_0_y_0_yyzzz_xxx, g_y_0_y_0_yyzzz_xxy, g_y_0_y_0_yyzzz_xxz, g_y_0_y_0_yyzzz_xy, g_y_0_y_0_yyzzz_xyy, g_y_0_y_0_yyzzz_xyz, g_y_0_y_0_yyzzz_xz, g_y_0_y_0_yyzzz_xzz, g_y_0_y_0_yyzzz_yy, g_y_0_y_0_yyzzz_yz, g_y_0_y_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyzzz_xx[k] = -g_y_0_y_0_yyzzz_xx[k] * ab_x + g_y_0_y_0_yyzzz_xxx[k];

                g_y_0_y_0_xyyzzz_xy[k] = -g_y_0_y_0_yyzzz_xy[k] * ab_x + g_y_0_y_0_yyzzz_xxy[k];

                g_y_0_y_0_xyyzzz_xz[k] = -g_y_0_y_0_yyzzz_xz[k] * ab_x + g_y_0_y_0_yyzzz_xxz[k];

                g_y_0_y_0_xyyzzz_yy[k] = -g_y_0_y_0_yyzzz_yy[k] * ab_x + g_y_0_y_0_yyzzz_xyy[k];

                g_y_0_y_0_xyyzzz_yz[k] = -g_y_0_y_0_yyzzz_yz[k] * ab_x + g_y_0_y_0_yyzzz_xyz[k];

                g_y_0_y_0_xyyzzz_zz[k] = -g_y_0_y_0_yyzzz_zz[k] * ab_x + g_y_0_y_0_yyzzz_xzz[k];
            }

            /// Set up 786-792 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyzzzz_xx = cbuffer.data(id_geom_1010_off + 786 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzzz_xy = cbuffer.data(id_geom_1010_off + 787 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzzz_xz = cbuffer.data(id_geom_1010_off + 788 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzzz_yy = cbuffer.data(id_geom_1010_off + 789 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzzz_yz = cbuffer.data(id_geom_1010_off + 790 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzzz_zz = cbuffer.data(id_geom_1010_off + 791 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyzzzz_xx, g_y_0_y_0_xyzzzz_xy, g_y_0_y_0_xyzzzz_xz, g_y_0_y_0_xyzzzz_yy, g_y_0_y_0_xyzzzz_yz, g_y_0_y_0_xyzzzz_zz, g_y_0_y_0_yzzzz_xx, g_y_0_y_0_yzzzz_xxx, g_y_0_y_0_yzzzz_xxy, g_y_0_y_0_yzzzz_xxz, g_y_0_y_0_yzzzz_xy, g_y_0_y_0_yzzzz_xyy, g_y_0_y_0_yzzzz_xyz, g_y_0_y_0_yzzzz_xz, g_y_0_y_0_yzzzz_xzz, g_y_0_y_0_yzzzz_yy, g_y_0_y_0_yzzzz_yz, g_y_0_y_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyzzzz_xx[k] = -g_y_0_y_0_yzzzz_xx[k] * ab_x + g_y_0_y_0_yzzzz_xxx[k];

                g_y_0_y_0_xyzzzz_xy[k] = -g_y_0_y_0_yzzzz_xy[k] * ab_x + g_y_0_y_0_yzzzz_xxy[k];

                g_y_0_y_0_xyzzzz_xz[k] = -g_y_0_y_0_yzzzz_xz[k] * ab_x + g_y_0_y_0_yzzzz_xxz[k];

                g_y_0_y_0_xyzzzz_yy[k] = -g_y_0_y_0_yzzzz_yy[k] * ab_x + g_y_0_y_0_yzzzz_xyy[k];

                g_y_0_y_0_xyzzzz_yz[k] = -g_y_0_y_0_yzzzz_yz[k] * ab_x + g_y_0_y_0_yzzzz_xyz[k];

                g_y_0_y_0_xyzzzz_zz[k] = -g_y_0_y_0_yzzzz_zz[k] * ab_x + g_y_0_y_0_yzzzz_xzz[k];
            }

            /// Set up 792-798 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xzzzzz_xx = cbuffer.data(id_geom_1010_off + 792 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzzz_xy = cbuffer.data(id_geom_1010_off + 793 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzzz_xz = cbuffer.data(id_geom_1010_off + 794 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzzz_yy = cbuffer.data(id_geom_1010_off + 795 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzzz_yz = cbuffer.data(id_geom_1010_off + 796 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzzz_zz = cbuffer.data(id_geom_1010_off + 797 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xzzzzz_xx, g_y_0_y_0_xzzzzz_xy, g_y_0_y_0_xzzzzz_xz, g_y_0_y_0_xzzzzz_yy, g_y_0_y_0_xzzzzz_yz, g_y_0_y_0_xzzzzz_zz, g_y_0_y_0_zzzzz_xx, g_y_0_y_0_zzzzz_xxx, g_y_0_y_0_zzzzz_xxy, g_y_0_y_0_zzzzz_xxz, g_y_0_y_0_zzzzz_xy, g_y_0_y_0_zzzzz_xyy, g_y_0_y_0_zzzzz_xyz, g_y_0_y_0_zzzzz_xz, g_y_0_y_0_zzzzz_xzz, g_y_0_y_0_zzzzz_yy, g_y_0_y_0_zzzzz_yz, g_y_0_y_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xzzzzz_xx[k] = -g_y_0_y_0_zzzzz_xx[k] * ab_x + g_y_0_y_0_zzzzz_xxx[k];

                g_y_0_y_0_xzzzzz_xy[k] = -g_y_0_y_0_zzzzz_xy[k] * ab_x + g_y_0_y_0_zzzzz_xxy[k];

                g_y_0_y_0_xzzzzz_xz[k] = -g_y_0_y_0_zzzzz_xz[k] * ab_x + g_y_0_y_0_zzzzz_xxz[k];

                g_y_0_y_0_xzzzzz_yy[k] = -g_y_0_y_0_zzzzz_yy[k] * ab_x + g_y_0_y_0_zzzzz_xyy[k];

                g_y_0_y_0_xzzzzz_yz[k] = -g_y_0_y_0_zzzzz_yz[k] * ab_x + g_y_0_y_0_zzzzz_xyz[k];

                g_y_0_y_0_xzzzzz_zz[k] = -g_y_0_y_0_zzzzz_zz[k] * ab_x + g_y_0_y_0_zzzzz_xzz[k];
            }

            /// Set up 798-804 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyyyy_xx = cbuffer.data(id_geom_1010_off + 798 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyy_xy = cbuffer.data(id_geom_1010_off + 799 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyy_xz = cbuffer.data(id_geom_1010_off + 800 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyy_yy = cbuffer.data(id_geom_1010_off + 801 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyy_yz = cbuffer.data(id_geom_1010_off + 802 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyy_zz = cbuffer.data(id_geom_1010_off + 803 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_yyyyy_xx, g_0_0_y_0_yyyyy_xy, g_0_0_y_0_yyyyy_xz, g_0_0_y_0_yyyyy_yy, g_0_0_y_0_yyyyy_yz, g_0_0_y_0_yyyyy_zz, g_y_0_y_0_yyyyy_xx, g_y_0_y_0_yyyyy_xxy, g_y_0_y_0_yyyyy_xy, g_y_0_y_0_yyyyy_xyy, g_y_0_y_0_yyyyy_xyz, g_y_0_y_0_yyyyy_xz, g_y_0_y_0_yyyyy_yy, g_y_0_y_0_yyyyy_yyy, g_y_0_y_0_yyyyy_yyz, g_y_0_y_0_yyyyy_yz, g_y_0_y_0_yyyyy_yzz, g_y_0_y_0_yyyyy_zz, g_y_0_y_0_yyyyyy_xx, g_y_0_y_0_yyyyyy_xy, g_y_0_y_0_yyyyyy_xz, g_y_0_y_0_yyyyyy_yy, g_y_0_y_0_yyyyyy_yz, g_y_0_y_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyyyy_xx[k] = -g_0_0_y_0_yyyyy_xx[k] - g_y_0_y_0_yyyyy_xx[k] * ab_y + g_y_0_y_0_yyyyy_xxy[k];

                g_y_0_y_0_yyyyyy_xy[k] = -g_0_0_y_0_yyyyy_xy[k] - g_y_0_y_0_yyyyy_xy[k] * ab_y + g_y_0_y_0_yyyyy_xyy[k];

                g_y_0_y_0_yyyyyy_xz[k] = -g_0_0_y_0_yyyyy_xz[k] - g_y_0_y_0_yyyyy_xz[k] * ab_y + g_y_0_y_0_yyyyy_xyz[k];

                g_y_0_y_0_yyyyyy_yy[k] = -g_0_0_y_0_yyyyy_yy[k] - g_y_0_y_0_yyyyy_yy[k] * ab_y + g_y_0_y_0_yyyyy_yyy[k];

                g_y_0_y_0_yyyyyy_yz[k] = -g_0_0_y_0_yyyyy_yz[k] - g_y_0_y_0_yyyyy_yz[k] * ab_y + g_y_0_y_0_yyyyy_yyz[k];

                g_y_0_y_0_yyyyyy_zz[k] = -g_0_0_y_0_yyyyy_zz[k] - g_y_0_y_0_yyyyy_zz[k] * ab_y + g_y_0_y_0_yyyyy_yzz[k];
            }

            /// Set up 804-810 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyyyz_xx = cbuffer.data(id_geom_1010_off + 804 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyz_xy = cbuffer.data(id_geom_1010_off + 805 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyz_xz = cbuffer.data(id_geom_1010_off + 806 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyz_yy = cbuffer.data(id_geom_1010_off + 807 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyz_yz = cbuffer.data(id_geom_1010_off + 808 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyz_zz = cbuffer.data(id_geom_1010_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyyyy_xx, g_y_0_y_0_yyyyy_xxz, g_y_0_y_0_yyyyy_xy, g_y_0_y_0_yyyyy_xyz, g_y_0_y_0_yyyyy_xz, g_y_0_y_0_yyyyy_xzz, g_y_0_y_0_yyyyy_yy, g_y_0_y_0_yyyyy_yyz, g_y_0_y_0_yyyyy_yz, g_y_0_y_0_yyyyy_yzz, g_y_0_y_0_yyyyy_zz, g_y_0_y_0_yyyyy_zzz, g_y_0_y_0_yyyyyz_xx, g_y_0_y_0_yyyyyz_xy, g_y_0_y_0_yyyyyz_xz, g_y_0_y_0_yyyyyz_yy, g_y_0_y_0_yyyyyz_yz, g_y_0_y_0_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyyyz_xx[k] = -g_y_0_y_0_yyyyy_xx[k] * ab_z + g_y_0_y_0_yyyyy_xxz[k];

                g_y_0_y_0_yyyyyz_xy[k] = -g_y_0_y_0_yyyyy_xy[k] * ab_z + g_y_0_y_0_yyyyy_xyz[k];

                g_y_0_y_0_yyyyyz_xz[k] = -g_y_0_y_0_yyyyy_xz[k] * ab_z + g_y_0_y_0_yyyyy_xzz[k];

                g_y_0_y_0_yyyyyz_yy[k] = -g_y_0_y_0_yyyyy_yy[k] * ab_z + g_y_0_y_0_yyyyy_yyz[k];

                g_y_0_y_0_yyyyyz_yz[k] = -g_y_0_y_0_yyyyy_yz[k] * ab_z + g_y_0_y_0_yyyyy_yzz[k];

                g_y_0_y_0_yyyyyz_zz[k] = -g_y_0_y_0_yyyyy_zz[k] * ab_z + g_y_0_y_0_yyyyy_zzz[k];
            }

            /// Set up 810-816 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyyzz_xx = cbuffer.data(id_geom_1010_off + 810 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyzz_xy = cbuffer.data(id_geom_1010_off + 811 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyzz_xz = cbuffer.data(id_geom_1010_off + 812 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyzz_yy = cbuffer.data(id_geom_1010_off + 813 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyzz_yz = cbuffer.data(id_geom_1010_off + 814 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyzz_zz = cbuffer.data(id_geom_1010_off + 815 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyyyz_xx, g_y_0_y_0_yyyyz_xxz, g_y_0_y_0_yyyyz_xy, g_y_0_y_0_yyyyz_xyz, g_y_0_y_0_yyyyz_xz, g_y_0_y_0_yyyyz_xzz, g_y_0_y_0_yyyyz_yy, g_y_0_y_0_yyyyz_yyz, g_y_0_y_0_yyyyz_yz, g_y_0_y_0_yyyyz_yzz, g_y_0_y_0_yyyyz_zz, g_y_0_y_0_yyyyz_zzz, g_y_0_y_0_yyyyzz_xx, g_y_0_y_0_yyyyzz_xy, g_y_0_y_0_yyyyzz_xz, g_y_0_y_0_yyyyzz_yy, g_y_0_y_0_yyyyzz_yz, g_y_0_y_0_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyyzz_xx[k] = -g_y_0_y_0_yyyyz_xx[k] * ab_z + g_y_0_y_0_yyyyz_xxz[k];

                g_y_0_y_0_yyyyzz_xy[k] = -g_y_0_y_0_yyyyz_xy[k] * ab_z + g_y_0_y_0_yyyyz_xyz[k];

                g_y_0_y_0_yyyyzz_xz[k] = -g_y_0_y_0_yyyyz_xz[k] * ab_z + g_y_0_y_0_yyyyz_xzz[k];

                g_y_0_y_0_yyyyzz_yy[k] = -g_y_0_y_0_yyyyz_yy[k] * ab_z + g_y_0_y_0_yyyyz_yyz[k];

                g_y_0_y_0_yyyyzz_yz[k] = -g_y_0_y_0_yyyyz_yz[k] * ab_z + g_y_0_y_0_yyyyz_yzz[k];

                g_y_0_y_0_yyyyzz_zz[k] = -g_y_0_y_0_yyyyz_zz[k] * ab_z + g_y_0_y_0_yyyyz_zzz[k];
            }

            /// Set up 816-822 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyzzz_xx = cbuffer.data(id_geom_1010_off + 816 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzzz_xy = cbuffer.data(id_geom_1010_off + 817 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzzz_xz = cbuffer.data(id_geom_1010_off + 818 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzzz_yy = cbuffer.data(id_geom_1010_off + 819 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzzz_yz = cbuffer.data(id_geom_1010_off + 820 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzzz_zz = cbuffer.data(id_geom_1010_off + 821 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyyzz_xx, g_y_0_y_0_yyyzz_xxz, g_y_0_y_0_yyyzz_xy, g_y_0_y_0_yyyzz_xyz, g_y_0_y_0_yyyzz_xz, g_y_0_y_0_yyyzz_xzz, g_y_0_y_0_yyyzz_yy, g_y_0_y_0_yyyzz_yyz, g_y_0_y_0_yyyzz_yz, g_y_0_y_0_yyyzz_yzz, g_y_0_y_0_yyyzz_zz, g_y_0_y_0_yyyzz_zzz, g_y_0_y_0_yyyzzz_xx, g_y_0_y_0_yyyzzz_xy, g_y_0_y_0_yyyzzz_xz, g_y_0_y_0_yyyzzz_yy, g_y_0_y_0_yyyzzz_yz, g_y_0_y_0_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyzzz_xx[k] = -g_y_0_y_0_yyyzz_xx[k] * ab_z + g_y_0_y_0_yyyzz_xxz[k];

                g_y_0_y_0_yyyzzz_xy[k] = -g_y_0_y_0_yyyzz_xy[k] * ab_z + g_y_0_y_0_yyyzz_xyz[k];

                g_y_0_y_0_yyyzzz_xz[k] = -g_y_0_y_0_yyyzz_xz[k] * ab_z + g_y_0_y_0_yyyzz_xzz[k];

                g_y_0_y_0_yyyzzz_yy[k] = -g_y_0_y_0_yyyzz_yy[k] * ab_z + g_y_0_y_0_yyyzz_yyz[k];

                g_y_0_y_0_yyyzzz_yz[k] = -g_y_0_y_0_yyyzz_yz[k] * ab_z + g_y_0_y_0_yyyzz_yzz[k];

                g_y_0_y_0_yyyzzz_zz[k] = -g_y_0_y_0_yyyzz_zz[k] * ab_z + g_y_0_y_0_yyyzz_zzz[k];
            }

            /// Set up 822-828 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyzzzz_xx = cbuffer.data(id_geom_1010_off + 822 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzzz_xy = cbuffer.data(id_geom_1010_off + 823 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzzz_xz = cbuffer.data(id_geom_1010_off + 824 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzzz_yy = cbuffer.data(id_geom_1010_off + 825 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzzz_yz = cbuffer.data(id_geom_1010_off + 826 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzzz_zz = cbuffer.data(id_geom_1010_off + 827 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyzzz_xx, g_y_0_y_0_yyzzz_xxz, g_y_0_y_0_yyzzz_xy, g_y_0_y_0_yyzzz_xyz, g_y_0_y_0_yyzzz_xz, g_y_0_y_0_yyzzz_xzz, g_y_0_y_0_yyzzz_yy, g_y_0_y_0_yyzzz_yyz, g_y_0_y_0_yyzzz_yz, g_y_0_y_0_yyzzz_yzz, g_y_0_y_0_yyzzz_zz, g_y_0_y_0_yyzzz_zzz, g_y_0_y_0_yyzzzz_xx, g_y_0_y_0_yyzzzz_xy, g_y_0_y_0_yyzzzz_xz, g_y_0_y_0_yyzzzz_yy, g_y_0_y_0_yyzzzz_yz, g_y_0_y_0_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyzzzz_xx[k] = -g_y_0_y_0_yyzzz_xx[k] * ab_z + g_y_0_y_0_yyzzz_xxz[k];

                g_y_0_y_0_yyzzzz_xy[k] = -g_y_0_y_0_yyzzz_xy[k] * ab_z + g_y_0_y_0_yyzzz_xyz[k];

                g_y_0_y_0_yyzzzz_xz[k] = -g_y_0_y_0_yyzzz_xz[k] * ab_z + g_y_0_y_0_yyzzz_xzz[k];

                g_y_0_y_0_yyzzzz_yy[k] = -g_y_0_y_0_yyzzz_yy[k] * ab_z + g_y_0_y_0_yyzzz_yyz[k];

                g_y_0_y_0_yyzzzz_yz[k] = -g_y_0_y_0_yyzzz_yz[k] * ab_z + g_y_0_y_0_yyzzz_yzz[k];

                g_y_0_y_0_yyzzzz_zz[k] = -g_y_0_y_0_yyzzz_zz[k] * ab_z + g_y_0_y_0_yyzzz_zzz[k];
            }

            /// Set up 828-834 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yzzzzz_xx = cbuffer.data(id_geom_1010_off + 828 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzzz_xy = cbuffer.data(id_geom_1010_off + 829 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzzz_xz = cbuffer.data(id_geom_1010_off + 830 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzzz_yy = cbuffer.data(id_geom_1010_off + 831 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzzz_yz = cbuffer.data(id_geom_1010_off + 832 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzzz_zz = cbuffer.data(id_geom_1010_off + 833 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yzzzz_xx, g_y_0_y_0_yzzzz_xxz, g_y_0_y_0_yzzzz_xy, g_y_0_y_0_yzzzz_xyz, g_y_0_y_0_yzzzz_xz, g_y_0_y_0_yzzzz_xzz, g_y_0_y_0_yzzzz_yy, g_y_0_y_0_yzzzz_yyz, g_y_0_y_0_yzzzz_yz, g_y_0_y_0_yzzzz_yzz, g_y_0_y_0_yzzzz_zz, g_y_0_y_0_yzzzz_zzz, g_y_0_y_0_yzzzzz_xx, g_y_0_y_0_yzzzzz_xy, g_y_0_y_0_yzzzzz_xz, g_y_0_y_0_yzzzzz_yy, g_y_0_y_0_yzzzzz_yz, g_y_0_y_0_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yzzzzz_xx[k] = -g_y_0_y_0_yzzzz_xx[k] * ab_z + g_y_0_y_0_yzzzz_xxz[k];

                g_y_0_y_0_yzzzzz_xy[k] = -g_y_0_y_0_yzzzz_xy[k] * ab_z + g_y_0_y_0_yzzzz_xyz[k];

                g_y_0_y_0_yzzzzz_xz[k] = -g_y_0_y_0_yzzzz_xz[k] * ab_z + g_y_0_y_0_yzzzz_xzz[k];

                g_y_0_y_0_yzzzzz_yy[k] = -g_y_0_y_0_yzzzz_yy[k] * ab_z + g_y_0_y_0_yzzzz_yyz[k];

                g_y_0_y_0_yzzzzz_yz[k] = -g_y_0_y_0_yzzzz_yz[k] * ab_z + g_y_0_y_0_yzzzz_yzz[k];

                g_y_0_y_0_yzzzzz_zz[k] = -g_y_0_y_0_yzzzz_zz[k] * ab_z + g_y_0_y_0_yzzzz_zzz[k];
            }

            /// Set up 834-840 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_zzzzzz_xx = cbuffer.data(id_geom_1010_off + 834 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzzz_xy = cbuffer.data(id_geom_1010_off + 835 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzzz_xz = cbuffer.data(id_geom_1010_off + 836 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzzz_yy = cbuffer.data(id_geom_1010_off + 837 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzzz_yz = cbuffer.data(id_geom_1010_off + 838 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzzz_zz = cbuffer.data(id_geom_1010_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_zzzzz_xx, g_y_0_y_0_zzzzz_xxz, g_y_0_y_0_zzzzz_xy, g_y_0_y_0_zzzzz_xyz, g_y_0_y_0_zzzzz_xz, g_y_0_y_0_zzzzz_xzz, g_y_0_y_0_zzzzz_yy, g_y_0_y_0_zzzzz_yyz, g_y_0_y_0_zzzzz_yz, g_y_0_y_0_zzzzz_yzz, g_y_0_y_0_zzzzz_zz, g_y_0_y_0_zzzzz_zzz, g_y_0_y_0_zzzzzz_xx, g_y_0_y_0_zzzzzz_xy, g_y_0_y_0_zzzzzz_xz, g_y_0_y_0_zzzzzz_yy, g_y_0_y_0_zzzzzz_yz, g_y_0_y_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_zzzzzz_xx[k] = -g_y_0_y_0_zzzzz_xx[k] * ab_z + g_y_0_y_0_zzzzz_xxz[k];

                g_y_0_y_0_zzzzzz_xy[k] = -g_y_0_y_0_zzzzz_xy[k] * ab_z + g_y_0_y_0_zzzzz_xyz[k];

                g_y_0_y_0_zzzzzz_xz[k] = -g_y_0_y_0_zzzzz_xz[k] * ab_z + g_y_0_y_0_zzzzz_xzz[k];

                g_y_0_y_0_zzzzzz_yy[k] = -g_y_0_y_0_zzzzz_yy[k] * ab_z + g_y_0_y_0_zzzzz_yyz[k];

                g_y_0_y_0_zzzzzz_yz[k] = -g_y_0_y_0_zzzzz_yz[k] * ab_z + g_y_0_y_0_zzzzz_yzz[k];

                g_y_0_y_0_zzzzzz_zz[k] = -g_y_0_y_0_zzzzz_zz[k] * ab_z + g_y_0_y_0_zzzzz_zzz[k];
            }

            /// Set up 840-846 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxxx_xx = cbuffer.data(id_geom_1010_off + 840 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxx_xy = cbuffer.data(id_geom_1010_off + 841 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxx_xz = cbuffer.data(id_geom_1010_off + 842 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxx_yy = cbuffer.data(id_geom_1010_off + 843 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxx_yz = cbuffer.data(id_geom_1010_off + 844 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxx_zz = cbuffer.data(id_geom_1010_off + 845 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxx_xx, g_y_0_z_0_xxxxx_xxx, g_y_0_z_0_xxxxx_xxy, g_y_0_z_0_xxxxx_xxz, g_y_0_z_0_xxxxx_xy, g_y_0_z_0_xxxxx_xyy, g_y_0_z_0_xxxxx_xyz, g_y_0_z_0_xxxxx_xz, g_y_0_z_0_xxxxx_xzz, g_y_0_z_0_xxxxx_yy, g_y_0_z_0_xxxxx_yz, g_y_0_z_0_xxxxx_zz, g_y_0_z_0_xxxxxx_xx, g_y_0_z_0_xxxxxx_xy, g_y_0_z_0_xxxxxx_xz, g_y_0_z_0_xxxxxx_yy, g_y_0_z_0_xxxxxx_yz, g_y_0_z_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxxx_xx[k] = -g_y_0_z_0_xxxxx_xx[k] * ab_x + g_y_0_z_0_xxxxx_xxx[k];

                g_y_0_z_0_xxxxxx_xy[k] = -g_y_0_z_0_xxxxx_xy[k] * ab_x + g_y_0_z_0_xxxxx_xxy[k];

                g_y_0_z_0_xxxxxx_xz[k] = -g_y_0_z_0_xxxxx_xz[k] * ab_x + g_y_0_z_0_xxxxx_xxz[k];

                g_y_0_z_0_xxxxxx_yy[k] = -g_y_0_z_0_xxxxx_yy[k] * ab_x + g_y_0_z_0_xxxxx_xyy[k];

                g_y_0_z_0_xxxxxx_yz[k] = -g_y_0_z_0_xxxxx_yz[k] * ab_x + g_y_0_z_0_xxxxx_xyz[k];

                g_y_0_z_0_xxxxxx_zz[k] = -g_y_0_z_0_xxxxx_zz[k] * ab_x + g_y_0_z_0_xxxxx_xzz[k];
            }

            /// Set up 846-852 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxxy_xx = cbuffer.data(id_geom_1010_off + 846 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxy_xy = cbuffer.data(id_geom_1010_off + 847 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxy_xz = cbuffer.data(id_geom_1010_off + 848 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxy_yy = cbuffer.data(id_geom_1010_off + 849 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxy_yz = cbuffer.data(id_geom_1010_off + 850 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxy_zz = cbuffer.data(id_geom_1010_off + 851 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxxy_xx, g_y_0_z_0_xxxxxy_xy, g_y_0_z_0_xxxxxy_xz, g_y_0_z_0_xxxxxy_yy, g_y_0_z_0_xxxxxy_yz, g_y_0_z_0_xxxxxy_zz, g_y_0_z_0_xxxxy_xx, g_y_0_z_0_xxxxy_xxx, g_y_0_z_0_xxxxy_xxy, g_y_0_z_0_xxxxy_xxz, g_y_0_z_0_xxxxy_xy, g_y_0_z_0_xxxxy_xyy, g_y_0_z_0_xxxxy_xyz, g_y_0_z_0_xxxxy_xz, g_y_0_z_0_xxxxy_xzz, g_y_0_z_0_xxxxy_yy, g_y_0_z_0_xxxxy_yz, g_y_0_z_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxxy_xx[k] = -g_y_0_z_0_xxxxy_xx[k] * ab_x + g_y_0_z_0_xxxxy_xxx[k];

                g_y_0_z_0_xxxxxy_xy[k] = -g_y_0_z_0_xxxxy_xy[k] * ab_x + g_y_0_z_0_xxxxy_xxy[k];

                g_y_0_z_0_xxxxxy_xz[k] = -g_y_0_z_0_xxxxy_xz[k] * ab_x + g_y_0_z_0_xxxxy_xxz[k];

                g_y_0_z_0_xxxxxy_yy[k] = -g_y_0_z_0_xxxxy_yy[k] * ab_x + g_y_0_z_0_xxxxy_xyy[k];

                g_y_0_z_0_xxxxxy_yz[k] = -g_y_0_z_0_xxxxy_yz[k] * ab_x + g_y_0_z_0_xxxxy_xyz[k];

                g_y_0_z_0_xxxxxy_zz[k] = -g_y_0_z_0_xxxxy_zz[k] * ab_x + g_y_0_z_0_xxxxy_xzz[k];
            }

            /// Set up 852-858 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxxz_xx = cbuffer.data(id_geom_1010_off + 852 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxz_xy = cbuffer.data(id_geom_1010_off + 853 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxz_xz = cbuffer.data(id_geom_1010_off + 854 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxz_yy = cbuffer.data(id_geom_1010_off + 855 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxz_yz = cbuffer.data(id_geom_1010_off + 856 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxz_zz = cbuffer.data(id_geom_1010_off + 857 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxxz_xx, g_y_0_z_0_xxxxxz_xy, g_y_0_z_0_xxxxxz_xz, g_y_0_z_0_xxxxxz_yy, g_y_0_z_0_xxxxxz_yz, g_y_0_z_0_xxxxxz_zz, g_y_0_z_0_xxxxz_xx, g_y_0_z_0_xxxxz_xxx, g_y_0_z_0_xxxxz_xxy, g_y_0_z_0_xxxxz_xxz, g_y_0_z_0_xxxxz_xy, g_y_0_z_0_xxxxz_xyy, g_y_0_z_0_xxxxz_xyz, g_y_0_z_0_xxxxz_xz, g_y_0_z_0_xxxxz_xzz, g_y_0_z_0_xxxxz_yy, g_y_0_z_0_xxxxz_yz, g_y_0_z_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxxz_xx[k] = -g_y_0_z_0_xxxxz_xx[k] * ab_x + g_y_0_z_0_xxxxz_xxx[k];

                g_y_0_z_0_xxxxxz_xy[k] = -g_y_0_z_0_xxxxz_xy[k] * ab_x + g_y_0_z_0_xxxxz_xxy[k];

                g_y_0_z_0_xxxxxz_xz[k] = -g_y_0_z_0_xxxxz_xz[k] * ab_x + g_y_0_z_0_xxxxz_xxz[k];

                g_y_0_z_0_xxxxxz_yy[k] = -g_y_0_z_0_xxxxz_yy[k] * ab_x + g_y_0_z_0_xxxxz_xyy[k];

                g_y_0_z_0_xxxxxz_yz[k] = -g_y_0_z_0_xxxxz_yz[k] * ab_x + g_y_0_z_0_xxxxz_xyz[k];

                g_y_0_z_0_xxxxxz_zz[k] = -g_y_0_z_0_xxxxz_zz[k] * ab_x + g_y_0_z_0_xxxxz_xzz[k];
            }

            /// Set up 858-864 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxyy_xx = cbuffer.data(id_geom_1010_off + 858 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyy_xy = cbuffer.data(id_geom_1010_off + 859 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyy_xz = cbuffer.data(id_geom_1010_off + 860 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyy_yy = cbuffer.data(id_geom_1010_off + 861 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyy_yz = cbuffer.data(id_geom_1010_off + 862 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyy_zz = cbuffer.data(id_geom_1010_off + 863 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxyy_xx, g_y_0_z_0_xxxxyy_xy, g_y_0_z_0_xxxxyy_xz, g_y_0_z_0_xxxxyy_yy, g_y_0_z_0_xxxxyy_yz, g_y_0_z_0_xxxxyy_zz, g_y_0_z_0_xxxyy_xx, g_y_0_z_0_xxxyy_xxx, g_y_0_z_0_xxxyy_xxy, g_y_0_z_0_xxxyy_xxz, g_y_0_z_0_xxxyy_xy, g_y_0_z_0_xxxyy_xyy, g_y_0_z_0_xxxyy_xyz, g_y_0_z_0_xxxyy_xz, g_y_0_z_0_xxxyy_xzz, g_y_0_z_0_xxxyy_yy, g_y_0_z_0_xxxyy_yz, g_y_0_z_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxyy_xx[k] = -g_y_0_z_0_xxxyy_xx[k] * ab_x + g_y_0_z_0_xxxyy_xxx[k];

                g_y_0_z_0_xxxxyy_xy[k] = -g_y_0_z_0_xxxyy_xy[k] * ab_x + g_y_0_z_0_xxxyy_xxy[k];

                g_y_0_z_0_xxxxyy_xz[k] = -g_y_0_z_0_xxxyy_xz[k] * ab_x + g_y_0_z_0_xxxyy_xxz[k];

                g_y_0_z_0_xxxxyy_yy[k] = -g_y_0_z_0_xxxyy_yy[k] * ab_x + g_y_0_z_0_xxxyy_xyy[k];

                g_y_0_z_0_xxxxyy_yz[k] = -g_y_0_z_0_xxxyy_yz[k] * ab_x + g_y_0_z_0_xxxyy_xyz[k];

                g_y_0_z_0_xxxxyy_zz[k] = -g_y_0_z_0_xxxyy_zz[k] * ab_x + g_y_0_z_0_xxxyy_xzz[k];
            }

            /// Set up 864-870 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxyz_xx = cbuffer.data(id_geom_1010_off + 864 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyz_xy = cbuffer.data(id_geom_1010_off + 865 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyz_xz = cbuffer.data(id_geom_1010_off + 866 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyz_yy = cbuffer.data(id_geom_1010_off + 867 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyz_yz = cbuffer.data(id_geom_1010_off + 868 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyz_zz = cbuffer.data(id_geom_1010_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxyz_xx, g_y_0_z_0_xxxxyz_xy, g_y_0_z_0_xxxxyz_xz, g_y_0_z_0_xxxxyz_yy, g_y_0_z_0_xxxxyz_yz, g_y_0_z_0_xxxxyz_zz, g_y_0_z_0_xxxyz_xx, g_y_0_z_0_xxxyz_xxx, g_y_0_z_0_xxxyz_xxy, g_y_0_z_0_xxxyz_xxz, g_y_0_z_0_xxxyz_xy, g_y_0_z_0_xxxyz_xyy, g_y_0_z_0_xxxyz_xyz, g_y_0_z_0_xxxyz_xz, g_y_0_z_0_xxxyz_xzz, g_y_0_z_0_xxxyz_yy, g_y_0_z_0_xxxyz_yz, g_y_0_z_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxyz_xx[k] = -g_y_0_z_0_xxxyz_xx[k] * ab_x + g_y_0_z_0_xxxyz_xxx[k];

                g_y_0_z_0_xxxxyz_xy[k] = -g_y_0_z_0_xxxyz_xy[k] * ab_x + g_y_0_z_0_xxxyz_xxy[k];

                g_y_0_z_0_xxxxyz_xz[k] = -g_y_0_z_0_xxxyz_xz[k] * ab_x + g_y_0_z_0_xxxyz_xxz[k];

                g_y_0_z_0_xxxxyz_yy[k] = -g_y_0_z_0_xxxyz_yy[k] * ab_x + g_y_0_z_0_xxxyz_xyy[k];

                g_y_0_z_0_xxxxyz_yz[k] = -g_y_0_z_0_xxxyz_yz[k] * ab_x + g_y_0_z_0_xxxyz_xyz[k];

                g_y_0_z_0_xxxxyz_zz[k] = -g_y_0_z_0_xxxyz_zz[k] * ab_x + g_y_0_z_0_xxxyz_xzz[k];
            }

            /// Set up 870-876 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxzz_xx = cbuffer.data(id_geom_1010_off + 870 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxzz_xy = cbuffer.data(id_geom_1010_off + 871 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxzz_xz = cbuffer.data(id_geom_1010_off + 872 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxzz_yy = cbuffer.data(id_geom_1010_off + 873 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxzz_yz = cbuffer.data(id_geom_1010_off + 874 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxzz_zz = cbuffer.data(id_geom_1010_off + 875 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxzz_xx, g_y_0_z_0_xxxxzz_xy, g_y_0_z_0_xxxxzz_xz, g_y_0_z_0_xxxxzz_yy, g_y_0_z_0_xxxxzz_yz, g_y_0_z_0_xxxxzz_zz, g_y_0_z_0_xxxzz_xx, g_y_0_z_0_xxxzz_xxx, g_y_0_z_0_xxxzz_xxy, g_y_0_z_0_xxxzz_xxz, g_y_0_z_0_xxxzz_xy, g_y_0_z_0_xxxzz_xyy, g_y_0_z_0_xxxzz_xyz, g_y_0_z_0_xxxzz_xz, g_y_0_z_0_xxxzz_xzz, g_y_0_z_0_xxxzz_yy, g_y_0_z_0_xxxzz_yz, g_y_0_z_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxzz_xx[k] = -g_y_0_z_0_xxxzz_xx[k] * ab_x + g_y_0_z_0_xxxzz_xxx[k];

                g_y_0_z_0_xxxxzz_xy[k] = -g_y_0_z_0_xxxzz_xy[k] * ab_x + g_y_0_z_0_xxxzz_xxy[k];

                g_y_0_z_0_xxxxzz_xz[k] = -g_y_0_z_0_xxxzz_xz[k] * ab_x + g_y_0_z_0_xxxzz_xxz[k];

                g_y_0_z_0_xxxxzz_yy[k] = -g_y_0_z_0_xxxzz_yy[k] * ab_x + g_y_0_z_0_xxxzz_xyy[k];

                g_y_0_z_0_xxxxzz_yz[k] = -g_y_0_z_0_xxxzz_yz[k] * ab_x + g_y_0_z_0_xxxzz_xyz[k];

                g_y_0_z_0_xxxxzz_zz[k] = -g_y_0_z_0_xxxzz_zz[k] * ab_x + g_y_0_z_0_xxxzz_xzz[k];
            }

            /// Set up 876-882 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxyyy_xx = cbuffer.data(id_geom_1010_off + 876 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyy_xy = cbuffer.data(id_geom_1010_off + 877 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyy_xz = cbuffer.data(id_geom_1010_off + 878 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyy_yy = cbuffer.data(id_geom_1010_off + 879 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyy_yz = cbuffer.data(id_geom_1010_off + 880 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyy_zz = cbuffer.data(id_geom_1010_off + 881 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxyyy_xx, g_y_0_z_0_xxxyyy_xy, g_y_0_z_0_xxxyyy_xz, g_y_0_z_0_xxxyyy_yy, g_y_0_z_0_xxxyyy_yz, g_y_0_z_0_xxxyyy_zz, g_y_0_z_0_xxyyy_xx, g_y_0_z_0_xxyyy_xxx, g_y_0_z_0_xxyyy_xxy, g_y_0_z_0_xxyyy_xxz, g_y_0_z_0_xxyyy_xy, g_y_0_z_0_xxyyy_xyy, g_y_0_z_0_xxyyy_xyz, g_y_0_z_0_xxyyy_xz, g_y_0_z_0_xxyyy_xzz, g_y_0_z_0_xxyyy_yy, g_y_0_z_0_xxyyy_yz, g_y_0_z_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxyyy_xx[k] = -g_y_0_z_0_xxyyy_xx[k] * ab_x + g_y_0_z_0_xxyyy_xxx[k];

                g_y_0_z_0_xxxyyy_xy[k] = -g_y_0_z_0_xxyyy_xy[k] * ab_x + g_y_0_z_0_xxyyy_xxy[k];

                g_y_0_z_0_xxxyyy_xz[k] = -g_y_0_z_0_xxyyy_xz[k] * ab_x + g_y_0_z_0_xxyyy_xxz[k];

                g_y_0_z_0_xxxyyy_yy[k] = -g_y_0_z_0_xxyyy_yy[k] * ab_x + g_y_0_z_0_xxyyy_xyy[k];

                g_y_0_z_0_xxxyyy_yz[k] = -g_y_0_z_0_xxyyy_yz[k] * ab_x + g_y_0_z_0_xxyyy_xyz[k];

                g_y_0_z_0_xxxyyy_zz[k] = -g_y_0_z_0_xxyyy_zz[k] * ab_x + g_y_0_z_0_xxyyy_xzz[k];
            }

            /// Set up 882-888 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxyyz_xx = cbuffer.data(id_geom_1010_off + 882 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyz_xy = cbuffer.data(id_geom_1010_off + 883 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyz_xz = cbuffer.data(id_geom_1010_off + 884 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyz_yy = cbuffer.data(id_geom_1010_off + 885 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyz_yz = cbuffer.data(id_geom_1010_off + 886 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyz_zz = cbuffer.data(id_geom_1010_off + 887 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxyyz_xx, g_y_0_z_0_xxxyyz_xy, g_y_0_z_0_xxxyyz_xz, g_y_0_z_0_xxxyyz_yy, g_y_0_z_0_xxxyyz_yz, g_y_0_z_0_xxxyyz_zz, g_y_0_z_0_xxyyz_xx, g_y_0_z_0_xxyyz_xxx, g_y_0_z_0_xxyyz_xxy, g_y_0_z_0_xxyyz_xxz, g_y_0_z_0_xxyyz_xy, g_y_0_z_0_xxyyz_xyy, g_y_0_z_0_xxyyz_xyz, g_y_0_z_0_xxyyz_xz, g_y_0_z_0_xxyyz_xzz, g_y_0_z_0_xxyyz_yy, g_y_0_z_0_xxyyz_yz, g_y_0_z_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxyyz_xx[k] = -g_y_0_z_0_xxyyz_xx[k] * ab_x + g_y_0_z_0_xxyyz_xxx[k];

                g_y_0_z_0_xxxyyz_xy[k] = -g_y_0_z_0_xxyyz_xy[k] * ab_x + g_y_0_z_0_xxyyz_xxy[k];

                g_y_0_z_0_xxxyyz_xz[k] = -g_y_0_z_0_xxyyz_xz[k] * ab_x + g_y_0_z_0_xxyyz_xxz[k];

                g_y_0_z_0_xxxyyz_yy[k] = -g_y_0_z_0_xxyyz_yy[k] * ab_x + g_y_0_z_0_xxyyz_xyy[k];

                g_y_0_z_0_xxxyyz_yz[k] = -g_y_0_z_0_xxyyz_yz[k] * ab_x + g_y_0_z_0_xxyyz_xyz[k];

                g_y_0_z_0_xxxyyz_zz[k] = -g_y_0_z_0_xxyyz_zz[k] * ab_x + g_y_0_z_0_xxyyz_xzz[k];
            }

            /// Set up 888-894 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxyzz_xx = cbuffer.data(id_geom_1010_off + 888 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyzz_xy = cbuffer.data(id_geom_1010_off + 889 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyzz_xz = cbuffer.data(id_geom_1010_off + 890 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyzz_yy = cbuffer.data(id_geom_1010_off + 891 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyzz_yz = cbuffer.data(id_geom_1010_off + 892 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyzz_zz = cbuffer.data(id_geom_1010_off + 893 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxyzz_xx, g_y_0_z_0_xxxyzz_xy, g_y_0_z_0_xxxyzz_xz, g_y_0_z_0_xxxyzz_yy, g_y_0_z_0_xxxyzz_yz, g_y_0_z_0_xxxyzz_zz, g_y_0_z_0_xxyzz_xx, g_y_0_z_0_xxyzz_xxx, g_y_0_z_0_xxyzz_xxy, g_y_0_z_0_xxyzz_xxz, g_y_0_z_0_xxyzz_xy, g_y_0_z_0_xxyzz_xyy, g_y_0_z_0_xxyzz_xyz, g_y_0_z_0_xxyzz_xz, g_y_0_z_0_xxyzz_xzz, g_y_0_z_0_xxyzz_yy, g_y_0_z_0_xxyzz_yz, g_y_0_z_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxyzz_xx[k] = -g_y_0_z_0_xxyzz_xx[k] * ab_x + g_y_0_z_0_xxyzz_xxx[k];

                g_y_0_z_0_xxxyzz_xy[k] = -g_y_0_z_0_xxyzz_xy[k] * ab_x + g_y_0_z_0_xxyzz_xxy[k];

                g_y_0_z_0_xxxyzz_xz[k] = -g_y_0_z_0_xxyzz_xz[k] * ab_x + g_y_0_z_0_xxyzz_xxz[k];

                g_y_0_z_0_xxxyzz_yy[k] = -g_y_0_z_0_xxyzz_yy[k] * ab_x + g_y_0_z_0_xxyzz_xyy[k];

                g_y_0_z_0_xxxyzz_yz[k] = -g_y_0_z_0_xxyzz_yz[k] * ab_x + g_y_0_z_0_xxyzz_xyz[k];

                g_y_0_z_0_xxxyzz_zz[k] = -g_y_0_z_0_xxyzz_zz[k] * ab_x + g_y_0_z_0_xxyzz_xzz[k];
            }

            /// Set up 894-900 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxzzz_xx = cbuffer.data(id_geom_1010_off + 894 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzzz_xy = cbuffer.data(id_geom_1010_off + 895 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzzz_xz = cbuffer.data(id_geom_1010_off + 896 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzzz_yy = cbuffer.data(id_geom_1010_off + 897 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzzz_yz = cbuffer.data(id_geom_1010_off + 898 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzzz_zz = cbuffer.data(id_geom_1010_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxzzz_xx, g_y_0_z_0_xxxzzz_xy, g_y_0_z_0_xxxzzz_xz, g_y_0_z_0_xxxzzz_yy, g_y_0_z_0_xxxzzz_yz, g_y_0_z_0_xxxzzz_zz, g_y_0_z_0_xxzzz_xx, g_y_0_z_0_xxzzz_xxx, g_y_0_z_0_xxzzz_xxy, g_y_0_z_0_xxzzz_xxz, g_y_0_z_0_xxzzz_xy, g_y_0_z_0_xxzzz_xyy, g_y_0_z_0_xxzzz_xyz, g_y_0_z_0_xxzzz_xz, g_y_0_z_0_xxzzz_xzz, g_y_0_z_0_xxzzz_yy, g_y_0_z_0_xxzzz_yz, g_y_0_z_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxzzz_xx[k] = -g_y_0_z_0_xxzzz_xx[k] * ab_x + g_y_0_z_0_xxzzz_xxx[k];

                g_y_0_z_0_xxxzzz_xy[k] = -g_y_0_z_0_xxzzz_xy[k] * ab_x + g_y_0_z_0_xxzzz_xxy[k];

                g_y_0_z_0_xxxzzz_xz[k] = -g_y_0_z_0_xxzzz_xz[k] * ab_x + g_y_0_z_0_xxzzz_xxz[k];

                g_y_0_z_0_xxxzzz_yy[k] = -g_y_0_z_0_xxzzz_yy[k] * ab_x + g_y_0_z_0_xxzzz_xyy[k];

                g_y_0_z_0_xxxzzz_yz[k] = -g_y_0_z_0_xxzzz_yz[k] * ab_x + g_y_0_z_0_xxzzz_xyz[k];

                g_y_0_z_0_xxxzzz_zz[k] = -g_y_0_z_0_xxzzz_zz[k] * ab_x + g_y_0_z_0_xxzzz_xzz[k];
            }

            /// Set up 900-906 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyyyy_xx = cbuffer.data(id_geom_1010_off + 900 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyy_xy = cbuffer.data(id_geom_1010_off + 901 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyy_xz = cbuffer.data(id_geom_1010_off + 902 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyy_yy = cbuffer.data(id_geom_1010_off + 903 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyy_yz = cbuffer.data(id_geom_1010_off + 904 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyy_zz = cbuffer.data(id_geom_1010_off + 905 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyyyy_xx, g_y_0_z_0_xxyyyy_xy, g_y_0_z_0_xxyyyy_xz, g_y_0_z_0_xxyyyy_yy, g_y_0_z_0_xxyyyy_yz, g_y_0_z_0_xxyyyy_zz, g_y_0_z_0_xyyyy_xx, g_y_0_z_0_xyyyy_xxx, g_y_0_z_0_xyyyy_xxy, g_y_0_z_0_xyyyy_xxz, g_y_0_z_0_xyyyy_xy, g_y_0_z_0_xyyyy_xyy, g_y_0_z_0_xyyyy_xyz, g_y_0_z_0_xyyyy_xz, g_y_0_z_0_xyyyy_xzz, g_y_0_z_0_xyyyy_yy, g_y_0_z_0_xyyyy_yz, g_y_0_z_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyyyy_xx[k] = -g_y_0_z_0_xyyyy_xx[k] * ab_x + g_y_0_z_0_xyyyy_xxx[k];

                g_y_0_z_0_xxyyyy_xy[k] = -g_y_0_z_0_xyyyy_xy[k] * ab_x + g_y_0_z_0_xyyyy_xxy[k];

                g_y_0_z_0_xxyyyy_xz[k] = -g_y_0_z_0_xyyyy_xz[k] * ab_x + g_y_0_z_0_xyyyy_xxz[k];

                g_y_0_z_0_xxyyyy_yy[k] = -g_y_0_z_0_xyyyy_yy[k] * ab_x + g_y_0_z_0_xyyyy_xyy[k];

                g_y_0_z_0_xxyyyy_yz[k] = -g_y_0_z_0_xyyyy_yz[k] * ab_x + g_y_0_z_0_xyyyy_xyz[k];

                g_y_0_z_0_xxyyyy_zz[k] = -g_y_0_z_0_xyyyy_zz[k] * ab_x + g_y_0_z_0_xyyyy_xzz[k];
            }

            /// Set up 906-912 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyyyz_xx = cbuffer.data(id_geom_1010_off + 906 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyz_xy = cbuffer.data(id_geom_1010_off + 907 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyz_xz = cbuffer.data(id_geom_1010_off + 908 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyz_yy = cbuffer.data(id_geom_1010_off + 909 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyz_yz = cbuffer.data(id_geom_1010_off + 910 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyz_zz = cbuffer.data(id_geom_1010_off + 911 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyyyz_xx, g_y_0_z_0_xxyyyz_xy, g_y_0_z_0_xxyyyz_xz, g_y_0_z_0_xxyyyz_yy, g_y_0_z_0_xxyyyz_yz, g_y_0_z_0_xxyyyz_zz, g_y_0_z_0_xyyyz_xx, g_y_0_z_0_xyyyz_xxx, g_y_0_z_0_xyyyz_xxy, g_y_0_z_0_xyyyz_xxz, g_y_0_z_0_xyyyz_xy, g_y_0_z_0_xyyyz_xyy, g_y_0_z_0_xyyyz_xyz, g_y_0_z_0_xyyyz_xz, g_y_0_z_0_xyyyz_xzz, g_y_0_z_0_xyyyz_yy, g_y_0_z_0_xyyyz_yz, g_y_0_z_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyyyz_xx[k] = -g_y_0_z_0_xyyyz_xx[k] * ab_x + g_y_0_z_0_xyyyz_xxx[k];

                g_y_0_z_0_xxyyyz_xy[k] = -g_y_0_z_0_xyyyz_xy[k] * ab_x + g_y_0_z_0_xyyyz_xxy[k];

                g_y_0_z_0_xxyyyz_xz[k] = -g_y_0_z_0_xyyyz_xz[k] * ab_x + g_y_0_z_0_xyyyz_xxz[k];

                g_y_0_z_0_xxyyyz_yy[k] = -g_y_0_z_0_xyyyz_yy[k] * ab_x + g_y_0_z_0_xyyyz_xyy[k];

                g_y_0_z_0_xxyyyz_yz[k] = -g_y_0_z_0_xyyyz_yz[k] * ab_x + g_y_0_z_0_xyyyz_xyz[k];

                g_y_0_z_0_xxyyyz_zz[k] = -g_y_0_z_0_xyyyz_zz[k] * ab_x + g_y_0_z_0_xyyyz_xzz[k];
            }

            /// Set up 912-918 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyyzz_xx = cbuffer.data(id_geom_1010_off + 912 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyzz_xy = cbuffer.data(id_geom_1010_off + 913 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyzz_xz = cbuffer.data(id_geom_1010_off + 914 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyzz_yy = cbuffer.data(id_geom_1010_off + 915 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyzz_yz = cbuffer.data(id_geom_1010_off + 916 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyzz_zz = cbuffer.data(id_geom_1010_off + 917 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyyzz_xx, g_y_0_z_0_xxyyzz_xy, g_y_0_z_0_xxyyzz_xz, g_y_0_z_0_xxyyzz_yy, g_y_0_z_0_xxyyzz_yz, g_y_0_z_0_xxyyzz_zz, g_y_0_z_0_xyyzz_xx, g_y_0_z_0_xyyzz_xxx, g_y_0_z_0_xyyzz_xxy, g_y_0_z_0_xyyzz_xxz, g_y_0_z_0_xyyzz_xy, g_y_0_z_0_xyyzz_xyy, g_y_0_z_0_xyyzz_xyz, g_y_0_z_0_xyyzz_xz, g_y_0_z_0_xyyzz_xzz, g_y_0_z_0_xyyzz_yy, g_y_0_z_0_xyyzz_yz, g_y_0_z_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyyzz_xx[k] = -g_y_0_z_0_xyyzz_xx[k] * ab_x + g_y_0_z_0_xyyzz_xxx[k];

                g_y_0_z_0_xxyyzz_xy[k] = -g_y_0_z_0_xyyzz_xy[k] * ab_x + g_y_0_z_0_xyyzz_xxy[k];

                g_y_0_z_0_xxyyzz_xz[k] = -g_y_0_z_0_xyyzz_xz[k] * ab_x + g_y_0_z_0_xyyzz_xxz[k];

                g_y_0_z_0_xxyyzz_yy[k] = -g_y_0_z_0_xyyzz_yy[k] * ab_x + g_y_0_z_0_xyyzz_xyy[k];

                g_y_0_z_0_xxyyzz_yz[k] = -g_y_0_z_0_xyyzz_yz[k] * ab_x + g_y_0_z_0_xyyzz_xyz[k];

                g_y_0_z_0_xxyyzz_zz[k] = -g_y_0_z_0_xyyzz_zz[k] * ab_x + g_y_0_z_0_xyyzz_xzz[k];
            }

            /// Set up 918-924 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyzzz_xx = cbuffer.data(id_geom_1010_off + 918 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzzz_xy = cbuffer.data(id_geom_1010_off + 919 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzzz_xz = cbuffer.data(id_geom_1010_off + 920 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzzz_yy = cbuffer.data(id_geom_1010_off + 921 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzzz_yz = cbuffer.data(id_geom_1010_off + 922 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzzz_zz = cbuffer.data(id_geom_1010_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyzzz_xx, g_y_0_z_0_xxyzzz_xy, g_y_0_z_0_xxyzzz_xz, g_y_0_z_0_xxyzzz_yy, g_y_0_z_0_xxyzzz_yz, g_y_0_z_0_xxyzzz_zz, g_y_0_z_0_xyzzz_xx, g_y_0_z_0_xyzzz_xxx, g_y_0_z_0_xyzzz_xxy, g_y_0_z_0_xyzzz_xxz, g_y_0_z_0_xyzzz_xy, g_y_0_z_0_xyzzz_xyy, g_y_0_z_0_xyzzz_xyz, g_y_0_z_0_xyzzz_xz, g_y_0_z_0_xyzzz_xzz, g_y_0_z_0_xyzzz_yy, g_y_0_z_0_xyzzz_yz, g_y_0_z_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyzzz_xx[k] = -g_y_0_z_0_xyzzz_xx[k] * ab_x + g_y_0_z_0_xyzzz_xxx[k];

                g_y_0_z_0_xxyzzz_xy[k] = -g_y_0_z_0_xyzzz_xy[k] * ab_x + g_y_0_z_0_xyzzz_xxy[k];

                g_y_0_z_0_xxyzzz_xz[k] = -g_y_0_z_0_xyzzz_xz[k] * ab_x + g_y_0_z_0_xyzzz_xxz[k];

                g_y_0_z_0_xxyzzz_yy[k] = -g_y_0_z_0_xyzzz_yy[k] * ab_x + g_y_0_z_0_xyzzz_xyy[k];

                g_y_0_z_0_xxyzzz_yz[k] = -g_y_0_z_0_xyzzz_yz[k] * ab_x + g_y_0_z_0_xyzzz_xyz[k];

                g_y_0_z_0_xxyzzz_zz[k] = -g_y_0_z_0_xyzzz_zz[k] * ab_x + g_y_0_z_0_xyzzz_xzz[k];
            }

            /// Set up 924-930 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxzzzz_xx = cbuffer.data(id_geom_1010_off + 924 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzzz_xy = cbuffer.data(id_geom_1010_off + 925 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzzz_xz = cbuffer.data(id_geom_1010_off + 926 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzzz_yy = cbuffer.data(id_geom_1010_off + 927 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzzz_yz = cbuffer.data(id_geom_1010_off + 928 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzzz_zz = cbuffer.data(id_geom_1010_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxzzzz_xx, g_y_0_z_0_xxzzzz_xy, g_y_0_z_0_xxzzzz_xz, g_y_0_z_0_xxzzzz_yy, g_y_0_z_0_xxzzzz_yz, g_y_0_z_0_xxzzzz_zz, g_y_0_z_0_xzzzz_xx, g_y_0_z_0_xzzzz_xxx, g_y_0_z_0_xzzzz_xxy, g_y_0_z_0_xzzzz_xxz, g_y_0_z_0_xzzzz_xy, g_y_0_z_0_xzzzz_xyy, g_y_0_z_0_xzzzz_xyz, g_y_0_z_0_xzzzz_xz, g_y_0_z_0_xzzzz_xzz, g_y_0_z_0_xzzzz_yy, g_y_0_z_0_xzzzz_yz, g_y_0_z_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxzzzz_xx[k] = -g_y_0_z_0_xzzzz_xx[k] * ab_x + g_y_0_z_0_xzzzz_xxx[k];

                g_y_0_z_0_xxzzzz_xy[k] = -g_y_0_z_0_xzzzz_xy[k] * ab_x + g_y_0_z_0_xzzzz_xxy[k];

                g_y_0_z_0_xxzzzz_xz[k] = -g_y_0_z_0_xzzzz_xz[k] * ab_x + g_y_0_z_0_xzzzz_xxz[k];

                g_y_0_z_0_xxzzzz_yy[k] = -g_y_0_z_0_xzzzz_yy[k] * ab_x + g_y_0_z_0_xzzzz_xyy[k];

                g_y_0_z_0_xxzzzz_yz[k] = -g_y_0_z_0_xzzzz_yz[k] * ab_x + g_y_0_z_0_xzzzz_xyz[k];

                g_y_0_z_0_xxzzzz_zz[k] = -g_y_0_z_0_xzzzz_zz[k] * ab_x + g_y_0_z_0_xzzzz_xzz[k];
            }

            /// Set up 930-936 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyyyy_xx = cbuffer.data(id_geom_1010_off + 930 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyy_xy = cbuffer.data(id_geom_1010_off + 931 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyy_xz = cbuffer.data(id_geom_1010_off + 932 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyy_yy = cbuffer.data(id_geom_1010_off + 933 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyy_yz = cbuffer.data(id_geom_1010_off + 934 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyy_zz = cbuffer.data(id_geom_1010_off + 935 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyyyy_xx, g_y_0_z_0_xyyyyy_xy, g_y_0_z_0_xyyyyy_xz, g_y_0_z_0_xyyyyy_yy, g_y_0_z_0_xyyyyy_yz, g_y_0_z_0_xyyyyy_zz, g_y_0_z_0_yyyyy_xx, g_y_0_z_0_yyyyy_xxx, g_y_0_z_0_yyyyy_xxy, g_y_0_z_0_yyyyy_xxz, g_y_0_z_0_yyyyy_xy, g_y_0_z_0_yyyyy_xyy, g_y_0_z_0_yyyyy_xyz, g_y_0_z_0_yyyyy_xz, g_y_0_z_0_yyyyy_xzz, g_y_0_z_0_yyyyy_yy, g_y_0_z_0_yyyyy_yz, g_y_0_z_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyyyy_xx[k] = -g_y_0_z_0_yyyyy_xx[k] * ab_x + g_y_0_z_0_yyyyy_xxx[k];

                g_y_0_z_0_xyyyyy_xy[k] = -g_y_0_z_0_yyyyy_xy[k] * ab_x + g_y_0_z_0_yyyyy_xxy[k];

                g_y_0_z_0_xyyyyy_xz[k] = -g_y_0_z_0_yyyyy_xz[k] * ab_x + g_y_0_z_0_yyyyy_xxz[k];

                g_y_0_z_0_xyyyyy_yy[k] = -g_y_0_z_0_yyyyy_yy[k] * ab_x + g_y_0_z_0_yyyyy_xyy[k];

                g_y_0_z_0_xyyyyy_yz[k] = -g_y_0_z_0_yyyyy_yz[k] * ab_x + g_y_0_z_0_yyyyy_xyz[k];

                g_y_0_z_0_xyyyyy_zz[k] = -g_y_0_z_0_yyyyy_zz[k] * ab_x + g_y_0_z_0_yyyyy_xzz[k];
            }

            /// Set up 936-942 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyyyz_xx = cbuffer.data(id_geom_1010_off + 936 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyz_xy = cbuffer.data(id_geom_1010_off + 937 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyz_xz = cbuffer.data(id_geom_1010_off + 938 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyz_yy = cbuffer.data(id_geom_1010_off + 939 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyz_yz = cbuffer.data(id_geom_1010_off + 940 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyz_zz = cbuffer.data(id_geom_1010_off + 941 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyyyz_xx, g_y_0_z_0_xyyyyz_xy, g_y_0_z_0_xyyyyz_xz, g_y_0_z_0_xyyyyz_yy, g_y_0_z_0_xyyyyz_yz, g_y_0_z_0_xyyyyz_zz, g_y_0_z_0_yyyyz_xx, g_y_0_z_0_yyyyz_xxx, g_y_0_z_0_yyyyz_xxy, g_y_0_z_0_yyyyz_xxz, g_y_0_z_0_yyyyz_xy, g_y_0_z_0_yyyyz_xyy, g_y_0_z_0_yyyyz_xyz, g_y_0_z_0_yyyyz_xz, g_y_0_z_0_yyyyz_xzz, g_y_0_z_0_yyyyz_yy, g_y_0_z_0_yyyyz_yz, g_y_0_z_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyyyz_xx[k] = -g_y_0_z_0_yyyyz_xx[k] * ab_x + g_y_0_z_0_yyyyz_xxx[k];

                g_y_0_z_0_xyyyyz_xy[k] = -g_y_0_z_0_yyyyz_xy[k] * ab_x + g_y_0_z_0_yyyyz_xxy[k];

                g_y_0_z_0_xyyyyz_xz[k] = -g_y_0_z_0_yyyyz_xz[k] * ab_x + g_y_0_z_0_yyyyz_xxz[k];

                g_y_0_z_0_xyyyyz_yy[k] = -g_y_0_z_0_yyyyz_yy[k] * ab_x + g_y_0_z_0_yyyyz_xyy[k];

                g_y_0_z_0_xyyyyz_yz[k] = -g_y_0_z_0_yyyyz_yz[k] * ab_x + g_y_0_z_0_yyyyz_xyz[k];

                g_y_0_z_0_xyyyyz_zz[k] = -g_y_0_z_0_yyyyz_zz[k] * ab_x + g_y_0_z_0_yyyyz_xzz[k];
            }

            /// Set up 942-948 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyyzz_xx = cbuffer.data(id_geom_1010_off + 942 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyzz_xy = cbuffer.data(id_geom_1010_off + 943 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyzz_xz = cbuffer.data(id_geom_1010_off + 944 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyzz_yy = cbuffer.data(id_geom_1010_off + 945 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyzz_yz = cbuffer.data(id_geom_1010_off + 946 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyzz_zz = cbuffer.data(id_geom_1010_off + 947 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyyzz_xx, g_y_0_z_0_xyyyzz_xy, g_y_0_z_0_xyyyzz_xz, g_y_0_z_0_xyyyzz_yy, g_y_0_z_0_xyyyzz_yz, g_y_0_z_0_xyyyzz_zz, g_y_0_z_0_yyyzz_xx, g_y_0_z_0_yyyzz_xxx, g_y_0_z_0_yyyzz_xxy, g_y_0_z_0_yyyzz_xxz, g_y_0_z_0_yyyzz_xy, g_y_0_z_0_yyyzz_xyy, g_y_0_z_0_yyyzz_xyz, g_y_0_z_0_yyyzz_xz, g_y_0_z_0_yyyzz_xzz, g_y_0_z_0_yyyzz_yy, g_y_0_z_0_yyyzz_yz, g_y_0_z_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyyzz_xx[k] = -g_y_0_z_0_yyyzz_xx[k] * ab_x + g_y_0_z_0_yyyzz_xxx[k];

                g_y_0_z_0_xyyyzz_xy[k] = -g_y_0_z_0_yyyzz_xy[k] * ab_x + g_y_0_z_0_yyyzz_xxy[k];

                g_y_0_z_0_xyyyzz_xz[k] = -g_y_0_z_0_yyyzz_xz[k] * ab_x + g_y_0_z_0_yyyzz_xxz[k];

                g_y_0_z_0_xyyyzz_yy[k] = -g_y_0_z_0_yyyzz_yy[k] * ab_x + g_y_0_z_0_yyyzz_xyy[k];

                g_y_0_z_0_xyyyzz_yz[k] = -g_y_0_z_0_yyyzz_yz[k] * ab_x + g_y_0_z_0_yyyzz_xyz[k];

                g_y_0_z_0_xyyyzz_zz[k] = -g_y_0_z_0_yyyzz_zz[k] * ab_x + g_y_0_z_0_yyyzz_xzz[k];
            }

            /// Set up 948-954 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyzzz_xx = cbuffer.data(id_geom_1010_off + 948 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzzz_xy = cbuffer.data(id_geom_1010_off + 949 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzzz_xz = cbuffer.data(id_geom_1010_off + 950 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzzz_yy = cbuffer.data(id_geom_1010_off + 951 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzzz_yz = cbuffer.data(id_geom_1010_off + 952 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzzz_zz = cbuffer.data(id_geom_1010_off + 953 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyzzz_xx, g_y_0_z_0_xyyzzz_xy, g_y_0_z_0_xyyzzz_xz, g_y_0_z_0_xyyzzz_yy, g_y_0_z_0_xyyzzz_yz, g_y_0_z_0_xyyzzz_zz, g_y_0_z_0_yyzzz_xx, g_y_0_z_0_yyzzz_xxx, g_y_0_z_0_yyzzz_xxy, g_y_0_z_0_yyzzz_xxz, g_y_0_z_0_yyzzz_xy, g_y_0_z_0_yyzzz_xyy, g_y_0_z_0_yyzzz_xyz, g_y_0_z_0_yyzzz_xz, g_y_0_z_0_yyzzz_xzz, g_y_0_z_0_yyzzz_yy, g_y_0_z_0_yyzzz_yz, g_y_0_z_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyzzz_xx[k] = -g_y_0_z_0_yyzzz_xx[k] * ab_x + g_y_0_z_0_yyzzz_xxx[k];

                g_y_0_z_0_xyyzzz_xy[k] = -g_y_0_z_0_yyzzz_xy[k] * ab_x + g_y_0_z_0_yyzzz_xxy[k];

                g_y_0_z_0_xyyzzz_xz[k] = -g_y_0_z_0_yyzzz_xz[k] * ab_x + g_y_0_z_0_yyzzz_xxz[k];

                g_y_0_z_0_xyyzzz_yy[k] = -g_y_0_z_0_yyzzz_yy[k] * ab_x + g_y_0_z_0_yyzzz_xyy[k];

                g_y_0_z_0_xyyzzz_yz[k] = -g_y_0_z_0_yyzzz_yz[k] * ab_x + g_y_0_z_0_yyzzz_xyz[k];

                g_y_0_z_0_xyyzzz_zz[k] = -g_y_0_z_0_yyzzz_zz[k] * ab_x + g_y_0_z_0_yyzzz_xzz[k];
            }

            /// Set up 954-960 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyzzzz_xx = cbuffer.data(id_geom_1010_off + 954 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzzz_xy = cbuffer.data(id_geom_1010_off + 955 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzzz_xz = cbuffer.data(id_geom_1010_off + 956 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzzz_yy = cbuffer.data(id_geom_1010_off + 957 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzzz_yz = cbuffer.data(id_geom_1010_off + 958 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzzz_zz = cbuffer.data(id_geom_1010_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyzzzz_xx, g_y_0_z_0_xyzzzz_xy, g_y_0_z_0_xyzzzz_xz, g_y_0_z_0_xyzzzz_yy, g_y_0_z_0_xyzzzz_yz, g_y_0_z_0_xyzzzz_zz, g_y_0_z_0_yzzzz_xx, g_y_0_z_0_yzzzz_xxx, g_y_0_z_0_yzzzz_xxy, g_y_0_z_0_yzzzz_xxz, g_y_0_z_0_yzzzz_xy, g_y_0_z_0_yzzzz_xyy, g_y_0_z_0_yzzzz_xyz, g_y_0_z_0_yzzzz_xz, g_y_0_z_0_yzzzz_xzz, g_y_0_z_0_yzzzz_yy, g_y_0_z_0_yzzzz_yz, g_y_0_z_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyzzzz_xx[k] = -g_y_0_z_0_yzzzz_xx[k] * ab_x + g_y_0_z_0_yzzzz_xxx[k];

                g_y_0_z_0_xyzzzz_xy[k] = -g_y_0_z_0_yzzzz_xy[k] * ab_x + g_y_0_z_0_yzzzz_xxy[k];

                g_y_0_z_0_xyzzzz_xz[k] = -g_y_0_z_0_yzzzz_xz[k] * ab_x + g_y_0_z_0_yzzzz_xxz[k];

                g_y_0_z_0_xyzzzz_yy[k] = -g_y_0_z_0_yzzzz_yy[k] * ab_x + g_y_0_z_0_yzzzz_xyy[k];

                g_y_0_z_0_xyzzzz_yz[k] = -g_y_0_z_0_yzzzz_yz[k] * ab_x + g_y_0_z_0_yzzzz_xyz[k];

                g_y_0_z_0_xyzzzz_zz[k] = -g_y_0_z_0_yzzzz_zz[k] * ab_x + g_y_0_z_0_yzzzz_xzz[k];
            }

            /// Set up 960-966 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xzzzzz_xx = cbuffer.data(id_geom_1010_off + 960 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzzz_xy = cbuffer.data(id_geom_1010_off + 961 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzzz_xz = cbuffer.data(id_geom_1010_off + 962 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzzz_yy = cbuffer.data(id_geom_1010_off + 963 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzzz_yz = cbuffer.data(id_geom_1010_off + 964 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzzz_zz = cbuffer.data(id_geom_1010_off + 965 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xzzzzz_xx, g_y_0_z_0_xzzzzz_xy, g_y_0_z_0_xzzzzz_xz, g_y_0_z_0_xzzzzz_yy, g_y_0_z_0_xzzzzz_yz, g_y_0_z_0_xzzzzz_zz, g_y_0_z_0_zzzzz_xx, g_y_0_z_0_zzzzz_xxx, g_y_0_z_0_zzzzz_xxy, g_y_0_z_0_zzzzz_xxz, g_y_0_z_0_zzzzz_xy, g_y_0_z_0_zzzzz_xyy, g_y_0_z_0_zzzzz_xyz, g_y_0_z_0_zzzzz_xz, g_y_0_z_0_zzzzz_xzz, g_y_0_z_0_zzzzz_yy, g_y_0_z_0_zzzzz_yz, g_y_0_z_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xzzzzz_xx[k] = -g_y_0_z_0_zzzzz_xx[k] * ab_x + g_y_0_z_0_zzzzz_xxx[k];

                g_y_0_z_0_xzzzzz_xy[k] = -g_y_0_z_0_zzzzz_xy[k] * ab_x + g_y_0_z_0_zzzzz_xxy[k];

                g_y_0_z_0_xzzzzz_xz[k] = -g_y_0_z_0_zzzzz_xz[k] * ab_x + g_y_0_z_0_zzzzz_xxz[k];

                g_y_0_z_0_xzzzzz_yy[k] = -g_y_0_z_0_zzzzz_yy[k] * ab_x + g_y_0_z_0_zzzzz_xyy[k];

                g_y_0_z_0_xzzzzz_yz[k] = -g_y_0_z_0_zzzzz_yz[k] * ab_x + g_y_0_z_0_zzzzz_xyz[k];

                g_y_0_z_0_xzzzzz_zz[k] = -g_y_0_z_0_zzzzz_zz[k] * ab_x + g_y_0_z_0_zzzzz_xzz[k];
            }

            /// Set up 966-972 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyyyy_xx = cbuffer.data(id_geom_1010_off + 966 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyy_xy = cbuffer.data(id_geom_1010_off + 967 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyy_xz = cbuffer.data(id_geom_1010_off + 968 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyy_yy = cbuffer.data(id_geom_1010_off + 969 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyy_yz = cbuffer.data(id_geom_1010_off + 970 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyy_zz = cbuffer.data(id_geom_1010_off + 971 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_yyyyy_xx, g_0_0_z_0_yyyyy_xy, g_0_0_z_0_yyyyy_xz, g_0_0_z_0_yyyyy_yy, g_0_0_z_0_yyyyy_yz, g_0_0_z_0_yyyyy_zz, g_y_0_z_0_yyyyy_xx, g_y_0_z_0_yyyyy_xxy, g_y_0_z_0_yyyyy_xy, g_y_0_z_0_yyyyy_xyy, g_y_0_z_0_yyyyy_xyz, g_y_0_z_0_yyyyy_xz, g_y_0_z_0_yyyyy_yy, g_y_0_z_0_yyyyy_yyy, g_y_0_z_0_yyyyy_yyz, g_y_0_z_0_yyyyy_yz, g_y_0_z_0_yyyyy_yzz, g_y_0_z_0_yyyyy_zz, g_y_0_z_0_yyyyyy_xx, g_y_0_z_0_yyyyyy_xy, g_y_0_z_0_yyyyyy_xz, g_y_0_z_0_yyyyyy_yy, g_y_0_z_0_yyyyyy_yz, g_y_0_z_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyyyy_xx[k] = -g_0_0_z_0_yyyyy_xx[k] - g_y_0_z_0_yyyyy_xx[k] * ab_y + g_y_0_z_0_yyyyy_xxy[k];

                g_y_0_z_0_yyyyyy_xy[k] = -g_0_0_z_0_yyyyy_xy[k] - g_y_0_z_0_yyyyy_xy[k] * ab_y + g_y_0_z_0_yyyyy_xyy[k];

                g_y_0_z_0_yyyyyy_xz[k] = -g_0_0_z_0_yyyyy_xz[k] - g_y_0_z_0_yyyyy_xz[k] * ab_y + g_y_0_z_0_yyyyy_xyz[k];

                g_y_0_z_0_yyyyyy_yy[k] = -g_0_0_z_0_yyyyy_yy[k] - g_y_0_z_0_yyyyy_yy[k] * ab_y + g_y_0_z_0_yyyyy_yyy[k];

                g_y_0_z_0_yyyyyy_yz[k] = -g_0_0_z_0_yyyyy_yz[k] - g_y_0_z_0_yyyyy_yz[k] * ab_y + g_y_0_z_0_yyyyy_yyz[k];

                g_y_0_z_0_yyyyyy_zz[k] = -g_0_0_z_0_yyyyy_zz[k] - g_y_0_z_0_yyyyy_zz[k] * ab_y + g_y_0_z_0_yyyyy_yzz[k];
            }

            /// Set up 972-978 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyyyz_xx = cbuffer.data(id_geom_1010_off + 972 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyz_xy = cbuffer.data(id_geom_1010_off + 973 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyz_xz = cbuffer.data(id_geom_1010_off + 974 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyz_yy = cbuffer.data(id_geom_1010_off + 975 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyz_yz = cbuffer.data(id_geom_1010_off + 976 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyz_zz = cbuffer.data(id_geom_1010_off + 977 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyyyy_xx, g_y_0_z_0_yyyyy_xxz, g_y_0_z_0_yyyyy_xy, g_y_0_z_0_yyyyy_xyz, g_y_0_z_0_yyyyy_xz, g_y_0_z_0_yyyyy_xzz, g_y_0_z_0_yyyyy_yy, g_y_0_z_0_yyyyy_yyz, g_y_0_z_0_yyyyy_yz, g_y_0_z_0_yyyyy_yzz, g_y_0_z_0_yyyyy_zz, g_y_0_z_0_yyyyy_zzz, g_y_0_z_0_yyyyyz_xx, g_y_0_z_0_yyyyyz_xy, g_y_0_z_0_yyyyyz_xz, g_y_0_z_0_yyyyyz_yy, g_y_0_z_0_yyyyyz_yz, g_y_0_z_0_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyyyz_xx[k] = -g_y_0_z_0_yyyyy_xx[k] * ab_z + g_y_0_z_0_yyyyy_xxz[k];

                g_y_0_z_0_yyyyyz_xy[k] = -g_y_0_z_0_yyyyy_xy[k] * ab_z + g_y_0_z_0_yyyyy_xyz[k];

                g_y_0_z_0_yyyyyz_xz[k] = -g_y_0_z_0_yyyyy_xz[k] * ab_z + g_y_0_z_0_yyyyy_xzz[k];

                g_y_0_z_0_yyyyyz_yy[k] = -g_y_0_z_0_yyyyy_yy[k] * ab_z + g_y_0_z_0_yyyyy_yyz[k];

                g_y_0_z_0_yyyyyz_yz[k] = -g_y_0_z_0_yyyyy_yz[k] * ab_z + g_y_0_z_0_yyyyy_yzz[k];

                g_y_0_z_0_yyyyyz_zz[k] = -g_y_0_z_0_yyyyy_zz[k] * ab_z + g_y_0_z_0_yyyyy_zzz[k];
            }

            /// Set up 978-984 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyyzz_xx = cbuffer.data(id_geom_1010_off + 978 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyzz_xy = cbuffer.data(id_geom_1010_off + 979 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyzz_xz = cbuffer.data(id_geom_1010_off + 980 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyzz_yy = cbuffer.data(id_geom_1010_off + 981 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyzz_yz = cbuffer.data(id_geom_1010_off + 982 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyzz_zz = cbuffer.data(id_geom_1010_off + 983 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyyyz_xx, g_y_0_z_0_yyyyz_xxz, g_y_0_z_0_yyyyz_xy, g_y_0_z_0_yyyyz_xyz, g_y_0_z_0_yyyyz_xz, g_y_0_z_0_yyyyz_xzz, g_y_0_z_0_yyyyz_yy, g_y_0_z_0_yyyyz_yyz, g_y_0_z_0_yyyyz_yz, g_y_0_z_0_yyyyz_yzz, g_y_0_z_0_yyyyz_zz, g_y_0_z_0_yyyyz_zzz, g_y_0_z_0_yyyyzz_xx, g_y_0_z_0_yyyyzz_xy, g_y_0_z_0_yyyyzz_xz, g_y_0_z_0_yyyyzz_yy, g_y_0_z_0_yyyyzz_yz, g_y_0_z_0_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyyzz_xx[k] = -g_y_0_z_0_yyyyz_xx[k] * ab_z + g_y_0_z_0_yyyyz_xxz[k];

                g_y_0_z_0_yyyyzz_xy[k] = -g_y_0_z_0_yyyyz_xy[k] * ab_z + g_y_0_z_0_yyyyz_xyz[k];

                g_y_0_z_0_yyyyzz_xz[k] = -g_y_0_z_0_yyyyz_xz[k] * ab_z + g_y_0_z_0_yyyyz_xzz[k];

                g_y_0_z_0_yyyyzz_yy[k] = -g_y_0_z_0_yyyyz_yy[k] * ab_z + g_y_0_z_0_yyyyz_yyz[k];

                g_y_0_z_0_yyyyzz_yz[k] = -g_y_0_z_0_yyyyz_yz[k] * ab_z + g_y_0_z_0_yyyyz_yzz[k];

                g_y_0_z_0_yyyyzz_zz[k] = -g_y_0_z_0_yyyyz_zz[k] * ab_z + g_y_0_z_0_yyyyz_zzz[k];
            }

            /// Set up 984-990 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyzzz_xx = cbuffer.data(id_geom_1010_off + 984 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzzz_xy = cbuffer.data(id_geom_1010_off + 985 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzzz_xz = cbuffer.data(id_geom_1010_off + 986 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzzz_yy = cbuffer.data(id_geom_1010_off + 987 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzzz_yz = cbuffer.data(id_geom_1010_off + 988 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzzz_zz = cbuffer.data(id_geom_1010_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyyzz_xx, g_y_0_z_0_yyyzz_xxz, g_y_0_z_0_yyyzz_xy, g_y_0_z_0_yyyzz_xyz, g_y_0_z_0_yyyzz_xz, g_y_0_z_0_yyyzz_xzz, g_y_0_z_0_yyyzz_yy, g_y_0_z_0_yyyzz_yyz, g_y_0_z_0_yyyzz_yz, g_y_0_z_0_yyyzz_yzz, g_y_0_z_0_yyyzz_zz, g_y_0_z_0_yyyzz_zzz, g_y_0_z_0_yyyzzz_xx, g_y_0_z_0_yyyzzz_xy, g_y_0_z_0_yyyzzz_xz, g_y_0_z_0_yyyzzz_yy, g_y_0_z_0_yyyzzz_yz, g_y_0_z_0_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyzzz_xx[k] = -g_y_0_z_0_yyyzz_xx[k] * ab_z + g_y_0_z_0_yyyzz_xxz[k];

                g_y_0_z_0_yyyzzz_xy[k] = -g_y_0_z_0_yyyzz_xy[k] * ab_z + g_y_0_z_0_yyyzz_xyz[k];

                g_y_0_z_0_yyyzzz_xz[k] = -g_y_0_z_0_yyyzz_xz[k] * ab_z + g_y_0_z_0_yyyzz_xzz[k];

                g_y_0_z_0_yyyzzz_yy[k] = -g_y_0_z_0_yyyzz_yy[k] * ab_z + g_y_0_z_0_yyyzz_yyz[k];

                g_y_0_z_0_yyyzzz_yz[k] = -g_y_0_z_0_yyyzz_yz[k] * ab_z + g_y_0_z_0_yyyzz_yzz[k];

                g_y_0_z_0_yyyzzz_zz[k] = -g_y_0_z_0_yyyzz_zz[k] * ab_z + g_y_0_z_0_yyyzz_zzz[k];
            }

            /// Set up 990-996 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyzzzz_xx = cbuffer.data(id_geom_1010_off + 990 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzzz_xy = cbuffer.data(id_geom_1010_off + 991 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzzz_xz = cbuffer.data(id_geom_1010_off + 992 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzzz_yy = cbuffer.data(id_geom_1010_off + 993 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzzz_yz = cbuffer.data(id_geom_1010_off + 994 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzzz_zz = cbuffer.data(id_geom_1010_off + 995 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyzzz_xx, g_y_0_z_0_yyzzz_xxz, g_y_0_z_0_yyzzz_xy, g_y_0_z_0_yyzzz_xyz, g_y_0_z_0_yyzzz_xz, g_y_0_z_0_yyzzz_xzz, g_y_0_z_0_yyzzz_yy, g_y_0_z_0_yyzzz_yyz, g_y_0_z_0_yyzzz_yz, g_y_0_z_0_yyzzz_yzz, g_y_0_z_0_yyzzz_zz, g_y_0_z_0_yyzzz_zzz, g_y_0_z_0_yyzzzz_xx, g_y_0_z_0_yyzzzz_xy, g_y_0_z_0_yyzzzz_xz, g_y_0_z_0_yyzzzz_yy, g_y_0_z_0_yyzzzz_yz, g_y_0_z_0_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyzzzz_xx[k] = -g_y_0_z_0_yyzzz_xx[k] * ab_z + g_y_0_z_0_yyzzz_xxz[k];

                g_y_0_z_0_yyzzzz_xy[k] = -g_y_0_z_0_yyzzz_xy[k] * ab_z + g_y_0_z_0_yyzzz_xyz[k];

                g_y_0_z_0_yyzzzz_xz[k] = -g_y_0_z_0_yyzzz_xz[k] * ab_z + g_y_0_z_0_yyzzz_xzz[k];

                g_y_0_z_0_yyzzzz_yy[k] = -g_y_0_z_0_yyzzz_yy[k] * ab_z + g_y_0_z_0_yyzzz_yyz[k];

                g_y_0_z_0_yyzzzz_yz[k] = -g_y_0_z_0_yyzzz_yz[k] * ab_z + g_y_0_z_0_yyzzz_yzz[k];

                g_y_0_z_0_yyzzzz_zz[k] = -g_y_0_z_0_yyzzz_zz[k] * ab_z + g_y_0_z_0_yyzzz_zzz[k];
            }

            /// Set up 996-1002 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yzzzzz_xx = cbuffer.data(id_geom_1010_off + 996 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzzz_xy = cbuffer.data(id_geom_1010_off + 997 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzzz_xz = cbuffer.data(id_geom_1010_off + 998 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzzz_yy = cbuffer.data(id_geom_1010_off + 999 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzzz_yz = cbuffer.data(id_geom_1010_off + 1000 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzzz_zz = cbuffer.data(id_geom_1010_off + 1001 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yzzzz_xx, g_y_0_z_0_yzzzz_xxz, g_y_0_z_0_yzzzz_xy, g_y_0_z_0_yzzzz_xyz, g_y_0_z_0_yzzzz_xz, g_y_0_z_0_yzzzz_xzz, g_y_0_z_0_yzzzz_yy, g_y_0_z_0_yzzzz_yyz, g_y_0_z_0_yzzzz_yz, g_y_0_z_0_yzzzz_yzz, g_y_0_z_0_yzzzz_zz, g_y_0_z_0_yzzzz_zzz, g_y_0_z_0_yzzzzz_xx, g_y_0_z_0_yzzzzz_xy, g_y_0_z_0_yzzzzz_xz, g_y_0_z_0_yzzzzz_yy, g_y_0_z_0_yzzzzz_yz, g_y_0_z_0_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yzzzzz_xx[k] = -g_y_0_z_0_yzzzz_xx[k] * ab_z + g_y_0_z_0_yzzzz_xxz[k];

                g_y_0_z_0_yzzzzz_xy[k] = -g_y_0_z_0_yzzzz_xy[k] * ab_z + g_y_0_z_0_yzzzz_xyz[k];

                g_y_0_z_0_yzzzzz_xz[k] = -g_y_0_z_0_yzzzz_xz[k] * ab_z + g_y_0_z_0_yzzzz_xzz[k];

                g_y_0_z_0_yzzzzz_yy[k] = -g_y_0_z_0_yzzzz_yy[k] * ab_z + g_y_0_z_0_yzzzz_yyz[k];

                g_y_0_z_0_yzzzzz_yz[k] = -g_y_0_z_0_yzzzz_yz[k] * ab_z + g_y_0_z_0_yzzzz_yzz[k];

                g_y_0_z_0_yzzzzz_zz[k] = -g_y_0_z_0_yzzzz_zz[k] * ab_z + g_y_0_z_0_yzzzz_zzz[k];
            }

            /// Set up 1002-1008 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_zzzzzz_xx = cbuffer.data(id_geom_1010_off + 1002 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzzz_xy = cbuffer.data(id_geom_1010_off + 1003 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzzz_xz = cbuffer.data(id_geom_1010_off + 1004 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzzz_yy = cbuffer.data(id_geom_1010_off + 1005 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzzz_yz = cbuffer.data(id_geom_1010_off + 1006 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzzz_zz = cbuffer.data(id_geom_1010_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_zzzzz_xx, g_y_0_z_0_zzzzz_xxz, g_y_0_z_0_zzzzz_xy, g_y_0_z_0_zzzzz_xyz, g_y_0_z_0_zzzzz_xz, g_y_0_z_0_zzzzz_xzz, g_y_0_z_0_zzzzz_yy, g_y_0_z_0_zzzzz_yyz, g_y_0_z_0_zzzzz_yz, g_y_0_z_0_zzzzz_yzz, g_y_0_z_0_zzzzz_zz, g_y_0_z_0_zzzzz_zzz, g_y_0_z_0_zzzzzz_xx, g_y_0_z_0_zzzzzz_xy, g_y_0_z_0_zzzzzz_xz, g_y_0_z_0_zzzzzz_yy, g_y_0_z_0_zzzzzz_yz, g_y_0_z_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_zzzzzz_xx[k] = -g_y_0_z_0_zzzzz_xx[k] * ab_z + g_y_0_z_0_zzzzz_xxz[k];

                g_y_0_z_0_zzzzzz_xy[k] = -g_y_0_z_0_zzzzz_xy[k] * ab_z + g_y_0_z_0_zzzzz_xyz[k];

                g_y_0_z_0_zzzzzz_xz[k] = -g_y_0_z_0_zzzzz_xz[k] * ab_z + g_y_0_z_0_zzzzz_xzz[k];

                g_y_0_z_0_zzzzzz_yy[k] = -g_y_0_z_0_zzzzz_yy[k] * ab_z + g_y_0_z_0_zzzzz_yyz[k];

                g_y_0_z_0_zzzzzz_yz[k] = -g_y_0_z_0_zzzzz_yz[k] * ab_z + g_y_0_z_0_zzzzz_yzz[k];

                g_y_0_z_0_zzzzzz_zz[k] = -g_y_0_z_0_zzzzz_zz[k] * ab_z + g_y_0_z_0_zzzzz_zzz[k];
            }

            /// Set up 1008-1014 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxxx_xx = cbuffer.data(id_geom_1010_off + 1008 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxx_xy = cbuffer.data(id_geom_1010_off + 1009 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxx_xz = cbuffer.data(id_geom_1010_off + 1010 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxx_yy = cbuffer.data(id_geom_1010_off + 1011 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxx_yz = cbuffer.data(id_geom_1010_off + 1012 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxx_zz = cbuffer.data(id_geom_1010_off + 1013 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxx_xx, g_z_0_x_0_xxxxx_xxx, g_z_0_x_0_xxxxx_xxy, g_z_0_x_0_xxxxx_xxz, g_z_0_x_0_xxxxx_xy, g_z_0_x_0_xxxxx_xyy, g_z_0_x_0_xxxxx_xyz, g_z_0_x_0_xxxxx_xz, g_z_0_x_0_xxxxx_xzz, g_z_0_x_0_xxxxx_yy, g_z_0_x_0_xxxxx_yz, g_z_0_x_0_xxxxx_zz, g_z_0_x_0_xxxxxx_xx, g_z_0_x_0_xxxxxx_xy, g_z_0_x_0_xxxxxx_xz, g_z_0_x_0_xxxxxx_yy, g_z_0_x_0_xxxxxx_yz, g_z_0_x_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxxx_xx[k] = -g_z_0_x_0_xxxxx_xx[k] * ab_x + g_z_0_x_0_xxxxx_xxx[k];

                g_z_0_x_0_xxxxxx_xy[k] = -g_z_0_x_0_xxxxx_xy[k] * ab_x + g_z_0_x_0_xxxxx_xxy[k];

                g_z_0_x_0_xxxxxx_xz[k] = -g_z_0_x_0_xxxxx_xz[k] * ab_x + g_z_0_x_0_xxxxx_xxz[k];

                g_z_0_x_0_xxxxxx_yy[k] = -g_z_0_x_0_xxxxx_yy[k] * ab_x + g_z_0_x_0_xxxxx_xyy[k];

                g_z_0_x_0_xxxxxx_yz[k] = -g_z_0_x_0_xxxxx_yz[k] * ab_x + g_z_0_x_0_xxxxx_xyz[k];

                g_z_0_x_0_xxxxxx_zz[k] = -g_z_0_x_0_xxxxx_zz[k] * ab_x + g_z_0_x_0_xxxxx_xzz[k];
            }

            /// Set up 1014-1020 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxxy_xx = cbuffer.data(id_geom_1010_off + 1014 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxy_xy = cbuffer.data(id_geom_1010_off + 1015 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxy_xz = cbuffer.data(id_geom_1010_off + 1016 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxy_yy = cbuffer.data(id_geom_1010_off + 1017 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxy_yz = cbuffer.data(id_geom_1010_off + 1018 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxy_zz = cbuffer.data(id_geom_1010_off + 1019 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxxy_xx, g_z_0_x_0_xxxxxy_xy, g_z_0_x_0_xxxxxy_xz, g_z_0_x_0_xxxxxy_yy, g_z_0_x_0_xxxxxy_yz, g_z_0_x_0_xxxxxy_zz, g_z_0_x_0_xxxxy_xx, g_z_0_x_0_xxxxy_xxx, g_z_0_x_0_xxxxy_xxy, g_z_0_x_0_xxxxy_xxz, g_z_0_x_0_xxxxy_xy, g_z_0_x_0_xxxxy_xyy, g_z_0_x_0_xxxxy_xyz, g_z_0_x_0_xxxxy_xz, g_z_0_x_0_xxxxy_xzz, g_z_0_x_0_xxxxy_yy, g_z_0_x_0_xxxxy_yz, g_z_0_x_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxxy_xx[k] = -g_z_0_x_0_xxxxy_xx[k] * ab_x + g_z_0_x_0_xxxxy_xxx[k];

                g_z_0_x_0_xxxxxy_xy[k] = -g_z_0_x_0_xxxxy_xy[k] * ab_x + g_z_0_x_0_xxxxy_xxy[k];

                g_z_0_x_0_xxxxxy_xz[k] = -g_z_0_x_0_xxxxy_xz[k] * ab_x + g_z_0_x_0_xxxxy_xxz[k];

                g_z_0_x_0_xxxxxy_yy[k] = -g_z_0_x_0_xxxxy_yy[k] * ab_x + g_z_0_x_0_xxxxy_xyy[k];

                g_z_0_x_0_xxxxxy_yz[k] = -g_z_0_x_0_xxxxy_yz[k] * ab_x + g_z_0_x_0_xxxxy_xyz[k];

                g_z_0_x_0_xxxxxy_zz[k] = -g_z_0_x_0_xxxxy_zz[k] * ab_x + g_z_0_x_0_xxxxy_xzz[k];
            }

            /// Set up 1020-1026 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxxz_xx = cbuffer.data(id_geom_1010_off + 1020 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxz_xy = cbuffer.data(id_geom_1010_off + 1021 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxz_xz = cbuffer.data(id_geom_1010_off + 1022 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxz_yy = cbuffer.data(id_geom_1010_off + 1023 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxz_yz = cbuffer.data(id_geom_1010_off + 1024 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxz_zz = cbuffer.data(id_geom_1010_off + 1025 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxxz_xx, g_z_0_x_0_xxxxxz_xy, g_z_0_x_0_xxxxxz_xz, g_z_0_x_0_xxxxxz_yy, g_z_0_x_0_xxxxxz_yz, g_z_0_x_0_xxxxxz_zz, g_z_0_x_0_xxxxz_xx, g_z_0_x_0_xxxxz_xxx, g_z_0_x_0_xxxxz_xxy, g_z_0_x_0_xxxxz_xxz, g_z_0_x_0_xxxxz_xy, g_z_0_x_0_xxxxz_xyy, g_z_0_x_0_xxxxz_xyz, g_z_0_x_0_xxxxz_xz, g_z_0_x_0_xxxxz_xzz, g_z_0_x_0_xxxxz_yy, g_z_0_x_0_xxxxz_yz, g_z_0_x_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxxz_xx[k] = -g_z_0_x_0_xxxxz_xx[k] * ab_x + g_z_0_x_0_xxxxz_xxx[k];

                g_z_0_x_0_xxxxxz_xy[k] = -g_z_0_x_0_xxxxz_xy[k] * ab_x + g_z_0_x_0_xxxxz_xxy[k];

                g_z_0_x_0_xxxxxz_xz[k] = -g_z_0_x_0_xxxxz_xz[k] * ab_x + g_z_0_x_0_xxxxz_xxz[k];

                g_z_0_x_0_xxxxxz_yy[k] = -g_z_0_x_0_xxxxz_yy[k] * ab_x + g_z_0_x_0_xxxxz_xyy[k];

                g_z_0_x_0_xxxxxz_yz[k] = -g_z_0_x_0_xxxxz_yz[k] * ab_x + g_z_0_x_0_xxxxz_xyz[k];

                g_z_0_x_0_xxxxxz_zz[k] = -g_z_0_x_0_xxxxz_zz[k] * ab_x + g_z_0_x_0_xxxxz_xzz[k];
            }

            /// Set up 1026-1032 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxyy_xx = cbuffer.data(id_geom_1010_off + 1026 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyy_xy = cbuffer.data(id_geom_1010_off + 1027 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyy_xz = cbuffer.data(id_geom_1010_off + 1028 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyy_yy = cbuffer.data(id_geom_1010_off + 1029 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyy_yz = cbuffer.data(id_geom_1010_off + 1030 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyy_zz = cbuffer.data(id_geom_1010_off + 1031 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxyy_xx, g_z_0_x_0_xxxxyy_xy, g_z_0_x_0_xxxxyy_xz, g_z_0_x_0_xxxxyy_yy, g_z_0_x_0_xxxxyy_yz, g_z_0_x_0_xxxxyy_zz, g_z_0_x_0_xxxyy_xx, g_z_0_x_0_xxxyy_xxx, g_z_0_x_0_xxxyy_xxy, g_z_0_x_0_xxxyy_xxz, g_z_0_x_0_xxxyy_xy, g_z_0_x_0_xxxyy_xyy, g_z_0_x_0_xxxyy_xyz, g_z_0_x_0_xxxyy_xz, g_z_0_x_0_xxxyy_xzz, g_z_0_x_0_xxxyy_yy, g_z_0_x_0_xxxyy_yz, g_z_0_x_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxyy_xx[k] = -g_z_0_x_0_xxxyy_xx[k] * ab_x + g_z_0_x_0_xxxyy_xxx[k];

                g_z_0_x_0_xxxxyy_xy[k] = -g_z_0_x_0_xxxyy_xy[k] * ab_x + g_z_0_x_0_xxxyy_xxy[k];

                g_z_0_x_0_xxxxyy_xz[k] = -g_z_0_x_0_xxxyy_xz[k] * ab_x + g_z_0_x_0_xxxyy_xxz[k];

                g_z_0_x_0_xxxxyy_yy[k] = -g_z_0_x_0_xxxyy_yy[k] * ab_x + g_z_0_x_0_xxxyy_xyy[k];

                g_z_0_x_0_xxxxyy_yz[k] = -g_z_0_x_0_xxxyy_yz[k] * ab_x + g_z_0_x_0_xxxyy_xyz[k];

                g_z_0_x_0_xxxxyy_zz[k] = -g_z_0_x_0_xxxyy_zz[k] * ab_x + g_z_0_x_0_xxxyy_xzz[k];
            }

            /// Set up 1032-1038 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxyz_xx = cbuffer.data(id_geom_1010_off + 1032 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyz_xy = cbuffer.data(id_geom_1010_off + 1033 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyz_xz = cbuffer.data(id_geom_1010_off + 1034 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyz_yy = cbuffer.data(id_geom_1010_off + 1035 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyz_yz = cbuffer.data(id_geom_1010_off + 1036 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyz_zz = cbuffer.data(id_geom_1010_off + 1037 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxyz_xx, g_z_0_x_0_xxxxyz_xy, g_z_0_x_0_xxxxyz_xz, g_z_0_x_0_xxxxyz_yy, g_z_0_x_0_xxxxyz_yz, g_z_0_x_0_xxxxyz_zz, g_z_0_x_0_xxxyz_xx, g_z_0_x_0_xxxyz_xxx, g_z_0_x_0_xxxyz_xxy, g_z_0_x_0_xxxyz_xxz, g_z_0_x_0_xxxyz_xy, g_z_0_x_0_xxxyz_xyy, g_z_0_x_0_xxxyz_xyz, g_z_0_x_0_xxxyz_xz, g_z_0_x_0_xxxyz_xzz, g_z_0_x_0_xxxyz_yy, g_z_0_x_0_xxxyz_yz, g_z_0_x_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxyz_xx[k] = -g_z_0_x_0_xxxyz_xx[k] * ab_x + g_z_0_x_0_xxxyz_xxx[k];

                g_z_0_x_0_xxxxyz_xy[k] = -g_z_0_x_0_xxxyz_xy[k] * ab_x + g_z_0_x_0_xxxyz_xxy[k];

                g_z_0_x_0_xxxxyz_xz[k] = -g_z_0_x_0_xxxyz_xz[k] * ab_x + g_z_0_x_0_xxxyz_xxz[k];

                g_z_0_x_0_xxxxyz_yy[k] = -g_z_0_x_0_xxxyz_yy[k] * ab_x + g_z_0_x_0_xxxyz_xyy[k];

                g_z_0_x_0_xxxxyz_yz[k] = -g_z_0_x_0_xxxyz_yz[k] * ab_x + g_z_0_x_0_xxxyz_xyz[k];

                g_z_0_x_0_xxxxyz_zz[k] = -g_z_0_x_0_xxxyz_zz[k] * ab_x + g_z_0_x_0_xxxyz_xzz[k];
            }

            /// Set up 1038-1044 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxzz_xx = cbuffer.data(id_geom_1010_off + 1038 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxzz_xy = cbuffer.data(id_geom_1010_off + 1039 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxzz_xz = cbuffer.data(id_geom_1010_off + 1040 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxzz_yy = cbuffer.data(id_geom_1010_off + 1041 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxzz_yz = cbuffer.data(id_geom_1010_off + 1042 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxzz_zz = cbuffer.data(id_geom_1010_off + 1043 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxzz_xx, g_z_0_x_0_xxxxzz_xy, g_z_0_x_0_xxxxzz_xz, g_z_0_x_0_xxxxzz_yy, g_z_0_x_0_xxxxzz_yz, g_z_0_x_0_xxxxzz_zz, g_z_0_x_0_xxxzz_xx, g_z_0_x_0_xxxzz_xxx, g_z_0_x_0_xxxzz_xxy, g_z_0_x_0_xxxzz_xxz, g_z_0_x_0_xxxzz_xy, g_z_0_x_0_xxxzz_xyy, g_z_0_x_0_xxxzz_xyz, g_z_0_x_0_xxxzz_xz, g_z_0_x_0_xxxzz_xzz, g_z_0_x_0_xxxzz_yy, g_z_0_x_0_xxxzz_yz, g_z_0_x_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxzz_xx[k] = -g_z_0_x_0_xxxzz_xx[k] * ab_x + g_z_0_x_0_xxxzz_xxx[k];

                g_z_0_x_0_xxxxzz_xy[k] = -g_z_0_x_0_xxxzz_xy[k] * ab_x + g_z_0_x_0_xxxzz_xxy[k];

                g_z_0_x_0_xxxxzz_xz[k] = -g_z_0_x_0_xxxzz_xz[k] * ab_x + g_z_0_x_0_xxxzz_xxz[k];

                g_z_0_x_0_xxxxzz_yy[k] = -g_z_0_x_0_xxxzz_yy[k] * ab_x + g_z_0_x_0_xxxzz_xyy[k];

                g_z_0_x_0_xxxxzz_yz[k] = -g_z_0_x_0_xxxzz_yz[k] * ab_x + g_z_0_x_0_xxxzz_xyz[k];

                g_z_0_x_0_xxxxzz_zz[k] = -g_z_0_x_0_xxxzz_zz[k] * ab_x + g_z_0_x_0_xxxzz_xzz[k];
            }

            /// Set up 1044-1050 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxyyy_xx = cbuffer.data(id_geom_1010_off + 1044 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyy_xy = cbuffer.data(id_geom_1010_off + 1045 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyy_xz = cbuffer.data(id_geom_1010_off + 1046 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyy_yy = cbuffer.data(id_geom_1010_off + 1047 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyy_yz = cbuffer.data(id_geom_1010_off + 1048 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyy_zz = cbuffer.data(id_geom_1010_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxyyy_xx, g_z_0_x_0_xxxyyy_xy, g_z_0_x_0_xxxyyy_xz, g_z_0_x_0_xxxyyy_yy, g_z_0_x_0_xxxyyy_yz, g_z_0_x_0_xxxyyy_zz, g_z_0_x_0_xxyyy_xx, g_z_0_x_0_xxyyy_xxx, g_z_0_x_0_xxyyy_xxy, g_z_0_x_0_xxyyy_xxz, g_z_0_x_0_xxyyy_xy, g_z_0_x_0_xxyyy_xyy, g_z_0_x_0_xxyyy_xyz, g_z_0_x_0_xxyyy_xz, g_z_0_x_0_xxyyy_xzz, g_z_0_x_0_xxyyy_yy, g_z_0_x_0_xxyyy_yz, g_z_0_x_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxyyy_xx[k] = -g_z_0_x_0_xxyyy_xx[k] * ab_x + g_z_0_x_0_xxyyy_xxx[k];

                g_z_0_x_0_xxxyyy_xy[k] = -g_z_0_x_0_xxyyy_xy[k] * ab_x + g_z_0_x_0_xxyyy_xxy[k];

                g_z_0_x_0_xxxyyy_xz[k] = -g_z_0_x_0_xxyyy_xz[k] * ab_x + g_z_0_x_0_xxyyy_xxz[k];

                g_z_0_x_0_xxxyyy_yy[k] = -g_z_0_x_0_xxyyy_yy[k] * ab_x + g_z_0_x_0_xxyyy_xyy[k];

                g_z_0_x_0_xxxyyy_yz[k] = -g_z_0_x_0_xxyyy_yz[k] * ab_x + g_z_0_x_0_xxyyy_xyz[k];

                g_z_0_x_0_xxxyyy_zz[k] = -g_z_0_x_0_xxyyy_zz[k] * ab_x + g_z_0_x_0_xxyyy_xzz[k];
            }

            /// Set up 1050-1056 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxyyz_xx = cbuffer.data(id_geom_1010_off + 1050 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyz_xy = cbuffer.data(id_geom_1010_off + 1051 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyz_xz = cbuffer.data(id_geom_1010_off + 1052 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyz_yy = cbuffer.data(id_geom_1010_off + 1053 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyz_yz = cbuffer.data(id_geom_1010_off + 1054 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyz_zz = cbuffer.data(id_geom_1010_off + 1055 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxyyz_xx, g_z_0_x_0_xxxyyz_xy, g_z_0_x_0_xxxyyz_xz, g_z_0_x_0_xxxyyz_yy, g_z_0_x_0_xxxyyz_yz, g_z_0_x_0_xxxyyz_zz, g_z_0_x_0_xxyyz_xx, g_z_0_x_0_xxyyz_xxx, g_z_0_x_0_xxyyz_xxy, g_z_0_x_0_xxyyz_xxz, g_z_0_x_0_xxyyz_xy, g_z_0_x_0_xxyyz_xyy, g_z_0_x_0_xxyyz_xyz, g_z_0_x_0_xxyyz_xz, g_z_0_x_0_xxyyz_xzz, g_z_0_x_0_xxyyz_yy, g_z_0_x_0_xxyyz_yz, g_z_0_x_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxyyz_xx[k] = -g_z_0_x_0_xxyyz_xx[k] * ab_x + g_z_0_x_0_xxyyz_xxx[k];

                g_z_0_x_0_xxxyyz_xy[k] = -g_z_0_x_0_xxyyz_xy[k] * ab_x + g_z_0_x_0_xxyyz_xxy[k];

                g_z_0_x_0_xxxyyz_xz[k] = -g_z_0_x_0_xxyyz_xz[k] * ab_x + g_z_0_x_0_xxyyz_xxz[k];

                g_z_0_x_0_xxxyyz_yy[k] = -g_z_0_x_0_xxyyz_yy[k] * ab_x + g_z_0_x_0_xxyyz_xyy[k];

                g_z_0_x_0_xxxyyz_yz[k] = -g_z_0_x_0_xxyyz_yz[k] * ab_x + g_z_0_x_0_xxyyz_xyz[k];

                g_z_0_x_0_xxxyyz_zz[k] = -g_z_0_x_0_xxyyz_zz[k] * ab_x + g_z_0_x_0_xxyyz_xzz[k];
            }

            /// Set up 1056-1062 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxyzz_xx = cbuffer.data(id_geom_1010_off + 1056 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyzz_xy = cbuffer.data(id_geom_1010_off + 1057 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyzz_xz = cbuffer.data(id_geom_1010_off + 1058 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyzz_yy = cbuffer.data(id_geom_1010_off + 1059 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyzz_yz = cbuffer.data(id_geom_1010_off + 1060 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyzz_zz = cbuffer.data(id_geom_1010_off + 1061 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxyzz_xx, g_z_0_x_0_xxxyzz_xy, g_z_0_x_0_xxxyzz_xz, g_z_0_x_0_xxxyzz_yy, g_z_0_x_0_xxxyzz_yz, g_z_0_x_0_xxxyzz_zz, g_z_0_x_0_xxyzz_xx, g_z_0_x_0_xxyzz_xxx, g_z_0_x_0_xxyzz_xxy, g_z_0_x_0_xxyzz_xxz, g_z_0_x_0_xxyzz_xy, g_z_0_x_0_xxyzz_xyy, g_z_0_x_0_xxyzz_xyz, g_z_0_x_0_xxyzz_xz, g_z_0_x_0_xxyzz_xzz, g_z_0_x_0_xxyzz_yy, g_z_0_x_0_xxyzz_yz, g_z_0_x_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxyzz_xx[k] = -g_z_0_x_0_xxyzz_xx[k] * ab_x + g_z_0_x_0_xxyzz_xxx[k];

                g_z_0_x_0_xxxyzz_xy[k] = -g_z_0_x_0_xxyzz_xy[k] * ab_x + g_z_0_x_0_xxyzz_xxy[k];

                g_z_0_x_0_xxxyzz_xz[k] = -g_z_0_x_0_xxyzz_xz[k] * ab_x + g_z_0_x_0_xxyzz_xxz[k];

                g_z_0_x_0_xxxyzz_yy[k] = -g_z_0_x_0_xxyzz_yy[k] * ab_x + g_z_0_x_0_xxyzz_xyy[k];

                g_z_0_x_0_xxxyzz_yz[k] = -g_z_0_x_0_xxyzz_yz[k] * ab_x + g_z_0_x_0_xxyzz_xyz[k];

                g_z_0_x_0_xxxyzz_zz[k] = -g_z_0_x_0_xxyzz_zz[k] * ab_x + g_z_0_x_0_xxyzz_xzz[k];
            }

            /// Set up 1062-1068 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxzzz_xx = cbuffer.data(id_geom_1010_off + 1062 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzzz_xy = cbuffer.data(id_geom_1010_off + 1063 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzzz_xz = cbuffer.data(id_geom_1010_off + 1064 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzzz_yy = cbuffer.data(id_geom_1010_off + 1065 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzzz_yz = cbuffer.data(id_geom_1010_off + 1066 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzzz_zz = cbuffer.data(id_geom_1010_off + 1067 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxzzz_xx, g_z_0_x_0_xxxzzz_xy, g_z_0_x_0_xxxzzz_xz, g_z_0_x_0_xxxzzz_yy, g_z_0_x_0_xxxzzz_yz, g_z_0_x_0_xxxzzz_zz, g_z_0_x_0_xxzzz_xx, g_z_0_x_0_xxzzz_xxx, g_z_0_x_0_xxzzz_xxy, g_z_0_x_0_xxzzz_xxz, g_z_0_x_0_xxzzz_xy, g_z_0_x_0_xxzzz_xyy, g_z_0_x_0_xxzzz_xyz, g_z_0_x_0_xxzzz_xz, g_z_0_x_0_xxzzz_xzz, g_z_0_x_0_xxzzz_yy, g_z_0_x_0_xxzzz_yz, g_z_0_x_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxzzz_xx[k] = -g_z_0_x_0_xxzzz_xx[k] * ab_x + g_z_0_x_0_xxzzz_xxx[k];

                g_z_0_x_0_xxxzzz_xy[k] = -g_z_0_x_0_xxzzz_xy[k] * ab_x + g_z_0_x_0_xxzzz_xxy[k];

                g_z_0_x_0_xxxzzz_xz[k] = -g_z_0_x_0_xxzzz_xz[k] * ab_x + g_z_0_x_0_xxzzz_xxz[k];

                g_z_0_x_0_xxxzzz_yy[k] = -g_z_0_x_0_xxzzz_yy[k] * ab_x + g_z_0_x_0_xxzzz_xyy[k];

                g_z_0_x_0_xxxzzz_yz[k] = -g_z_0_x_0_xxzzz_yz[k] * ab_x + g_z_0_x_0_xxzzz_xyz[k];

                g_z_0_x_0_xxxzzz_zz[k] = -g_z_0_x_0_xxzzz_zz[k] * ab_x + g_z_0_x_0_xxzzz_xzz[k];
            }

            /// Set up 1068-1074 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyyyy_xx = cbuffer.data(id_geom_1010_off + 1068 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyy_xy = cbuffer.data(id_geom_1010_off + 1069 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyy_xz = cbuffer.data(id_geom_1010_off + 1070 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyy_yy = cbuffer.data(id_geom_1010_off + 1071 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyy_yz = cbuffer.data(id_geom_1010_off + 1072 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyy_zz = cbuffer.data(id_geom_1010_off + 1073 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyyyy_xx, g_z_0_x_0_xxyyyy_xy, g_z_0_x_0_xxyyyy_xz, g_z_0_x_0_xxyyyy_yy, g_z_0_x_0_xxyyyy_yz, g_z_0_x_0_xxyyyy_zz, g_z_0_x_0_xyyyy_xx, g_z_0_x_0_xyyyy_xxx, g_z_0_x_0_xyyyy_xxy, g_z_0_x_0_xyyyy_xxz, g_z_0_x_0_xyyyy_xy, g_z_0_x_0_xyyyy_xyy, g_z_0_x_0_xyyyy_xyz, g_z_0_x_0_xyyyy_xz, g_z_0_x_0_xyyyy_xzz, g_z_0_x_0_xyyyy_yy, g_z_0_x_0_xyyyy_yz, g_z_0_x_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyyyy_xx[k] = -g_z_0_x_0_xyyyy_xx[k] * ab_x + g_z_0_x_0_xyyyy_xxx[k];

                g_z_0_x_0_xxyyyy_xy[k] = -g_z_0_x_0_xyyyy_xy[k] * ab_x + g_z_0_x_0_xyyyy_xxy[k];

                g_z_0_x_0_xxyyyy_xz[k] = -g_z_0_x_0_xyyyy_xz[k] * ab_x + g_z_0_x_0_xyyyy_xxz[k];

                g_z_0_x_0_xxyyyy_yy[k] = -g_z_0_x_0_xyyyy_yy[k] * ab_x + g_z_0_x_0_xyyyy_xyy[k];

                g_z_0_x_0_xxyyyy_yz[k] = -g_z_0_x_0_xyyyy_yz[k] * ab_x + g_z_0_x_0_xyyyy_xyz[k];

                g_z_0_x_0_xxyyyy_zz[k] = -g_z_0_x_0_xyyyy_zz[k] * ab_x + g_z_0_x_0_xyyyy_xzz[k];
            }

            /// Set up 1074-1080 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyyyz_xx = cbuffer.data(id_geom_1010_off + 1074 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyz_xy = cbuffer.data(id_geom_1010_off + 1075 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyz_xz = cbuffer.data(id_geom_1010_off + 1076 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyz_yy = cbuffer.data(id_geom_1010_off + 1077 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyz_yz = cbuffer.data(id_geom_1010_off + 1078 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyz_zz = cbuffer.data(id_geom_1010_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyyyz_xx, g_z_0_x_0_xxyyyz_xy, g_z_0_x_0_xxyyyz_xz, g_z_0_x_0_xxyyyz_yy, g_z_0_x_0_xxyyyz_yz, g_z_0_x_0_xxyyyz_zz, g_z_0_x_0_xyyyz_xx, g_z_0_x_0_xyyyz_xxx, g_z_0_x_0_xyyyz_xxy, g_z_0_x_0_xyyyz_xxz, g_z_0_x_0_xyyyz_xy, g_z_0_x_0_xyyyz_xyy, g_z_0_x_0_xyyyz_xyz, g_z_0_x_0_xyyyz_xz, g_z_0_x_0_xyyyz_xzz, g_z_0_x_0_xyyyz_yy, g_z_0_x_0_xyyyz_yz, g_z_0_x_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyyyz_xx[k] = -g_z_0_x_0_xyyyz_xx[k] * ab_x + g_z_0_x_0_xyyyz_xxx[k];

                g_z_0_x_0_xxyyyz_xy[k] = -g_z_0_x_0_xyyyz_xy[k] * ab_x + g_z_0_x_0_xyyyz_xxy[k];

                g_z_0_x_0_xxyyyz_xz[k] = -g_z_0_x_0_xyyyz_xz[k] * ab_x + g_z_0_x_0_xyyyz_xxz[k];

                g_z_0_x_0_xxyyyz_yy[k] = -g_z_0_x_0_xyyyz_yy[k] * ab_x + g_z_0_x_0_xyyyz_xyy[k];

                g_z_0_x_0_xxyyyz_yz[k] = -g_z_0_x_0_xyyyz_yz[k] * ab_x + g_z_0_x_0_xyyyz_xyz[k];

                g_z_0_x_0_xxyyyz_zz[k] = -g_z_0_x_0_xyyyz_zz[k] * ab_x + g_z_0_x_0_xyyyz_xzz[k];
            }

            /// Set up 1080-1086 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyyzz_xx = cbuffer.data(id_geom_1010_off + 1080 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyzz_xy = cbuffer.data(id_geom_1010_off + 1081 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyzz_xz = cbuffer.data(id_geom_1010_off + 1082 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyzz_yy = cbuffer.data(id_geom_1010_off + 1083 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyzz_yz = cbuffer.data(id_geom_1010_off + 1084 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyzz_zz = cbuffer.data(id_geom_1010_off + 1085 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyyzz_xx, g_z_0_x_0_xxyyzz_xy, g_z_0_x_0_xxyyzz_xz, g_z_0_x_0_xxyyzz_yy, g_z_0_x_0_xxyyzz_yz, g_z_0_x_0_xxyyzz_zz, g_z_0_x_0_xyyzz_xx, g_z_0_x_0_xyyzz_xxx, g_z_0_x_0_xyyzz_xxy, g_z_0_x_0_xyyzz_xxz, g_z_0_x_0_xyyzz_xy, g_z_0_x_0_xyyzz_xyy, g_z_0_x_0_xyyzz_xyz, g_z_0_x_0_xyyzz_xz, g_z_0_x_0_xyyzz_xzz, g_z_0_x_0_xyyzz_yy, g_z_0_x_0_xyyzz_yz, g_z_0_x_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyyzz_xx[k] = -g_z_0_x_0_xyyzz_xx[k] * ab_x + g_z_0_x_0_xyyzz_xxx[k];

                g_z_0_x_0_xxyyzz_xy[k] = -g_z_0_x_0_xyyzz_xy[k] * ab_x + g_z_0_x_0_xyyzz_xxy[k];

                g_z_0_x_0_xxyyzz_xz[k] = -g_z_0_x_0_xyyzz_xz[k] * ab_x + g_z_0_x_0_xyyzz_xxz[k];

                g_z_0_x_0_xxyyzz_yy[k] = -g_z_0_x_0_xyyzz_yy[k] * ab_x + g_z_0_x_0_xyyzz_xyy[k];

                g_z_0_x_0_xxyyzz_yz[k] = -g_z_0_x_0_xyyzz_yz[k] * ab_x + g_z_0_x_0_xyyzz_xyz[k];

                g_z_0_x_0_xxyyzz_zz[k] = -g_z_0_x_0_xyyzz_zz[k] * ab_x + g_z_0_x_0_xyyzz_xzz[k];
            }

            /// Set up 1086-1092 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyzzz_xx = cbuffer.data(id_geom_1010_off + 1086 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzzz_xy = cbuffer.data(id_geom_1010_off + 1087 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzzz_xz = cbuffer.data(id_geom_1010_off + 1088 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzzz_yy = cbuffer.data(id_geom_1010_off + 1089 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzzz_yz = cbuffer.data(id_geom_1010_off + 1090 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzzz_zz = cbuffer.data(id_geom_1010_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyzzz_xx, g_z_0_x_0_xxyzzz_xy, g_z_0_x_0_xxyzzz_xz, g_z_0_x_0_xxyzzz_yy, g_z_0_x_0_xxyzzz_yz, g_z_0_x_0_xxyzzz_zz, g_z_0_x_0_xyzzz_xx, g_z_0_x_0_xyzzz_xxx, g_z_0_x_0_xyzzz_xxy, g_z_0_x_0_xyzzz_xxz, g_z_0_x_0_xyzzz_xy, g_z_0_x_0_xyzzz_xyy, g_z_0_x_0_xyzzz_xyz, g_z_0_x_0_xyzzz_xz, g_z_0_x_0_xyzzz_xzz, g_z_0_x_0_xyzzz_yy, g_z_0_x_0_xyzzz_yz, g_z_0_x_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyzzz_xx[k] = -g_z_0_x_0_xyzzz_xx[k] * ab_x + g_z_0_x_0_xyzzz_xxx[k];

                g_z_0_x_0_xxyzzz_xy[k] = -g_z_0_x_0_xyzzz_xy[k] * ab_x + g_z_0_x_0_xyzzz_xxy[k];

                g_z_0_x_0_xxyzzz_xz[k] = -g_z_0_x_0_xyzzz_xz[k] * ab_x + g_z_0_x_0_xyzzz_xxz[k];

                g_z_0_x_0_xxyzzz_yy[k] = -g_z_0_x_0_xyzzz_yy[k] * ab_x + g_z_0_x_0_xyzzz_xyy[k];

                g_z_0_x_0_xxyzzz_yz[k] = -g_z_0_x_0_xyzzz_yz[k] * ab_x + g_z_0_x_0_xyzzz_xyz[k];

                g_z_0_x_0_xxyzzz_zz[k] = -g_z_0_x_0_xyzzz_zz[k] * ab_x + g_z_0_x_0_xyzzz_xzz[k];
            }

            /// Set up 1092-1098 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxzzzz_xx = cbuffer.data(id_geom_1010_off + 1092 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzzz_xy = cbuffer.data(id_geom_1010_off + 1093 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzzz_xz = cbuffer.data(id_geom_1010_off + 1094 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzzz_yy = cbuffer.data(id_geom_1010_off + 1095 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzzz_yz = cbuffer.data(id_geom_1010_off + 1096 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzzz_zz = cbuffer.data(id_geom_1010_off + 1097 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxzzzz_xx, g_z_0_x_0_xxzzzz_xy, g_z_0_x_0_xxzzzz_xz, g_z_0_x_0_xxzzzz_yy, g_z_0_x_0_xxzzzz_yz, g_z_0_x_0_xxzzzz_zz, g_z_0_x_0_xzzzz_xx, g_z_0_x_0_xzzzz_xxx, g_z_0_x_0_xzzzz_xxy, g_z_0_x_0_xzzzz_xxz, g_z_0_x_0_xzzzz_xy, g_z_0_x_0_xzzzz_xyy, g_z_0_x_0_xzzzz_xyz, g_z_0_x_0_xzzzz_xz, g_z_0_x_0_xzzzz_xzz, g_z_0_x_0_xzzzz_yy, g_z_0_x_0_xzzzz_yz, g_z_0_x_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxzzzz_xx[k] = -g_z_0_x_0_xzzzz_xx[k] * ab_x + g_z_0_x_0_xzzzz_xxx[k];

                g_z_0_x_0_xxzzzz_xy[k] = -g_z_0_x_0_xzzzz_xy[k] * ab_x + g_z_0_x_0_xzzzz_xxy[k];

                g_z_0_x_0_xxzzzz_xz[k] = -g_z_0_x_0_xzzzz_xz[k] * ab_x + g_z_0_x_0_xzzzz_xxz[k];

                g_z_0_x_0_xxzzzz_yy[k] = -g_z_0_x_0_xzzzz_yy[k] * ab_x + g_z_0_x_0_xzzzz_xyy[k];

                g_z_0_x_0_xxzzzz_yz[k] = -g_z_0_x_0_xzzzz_yz[k] * ab_x + g_z_0_x_0_xzzzz_xyz[k];

                g_z_0_x_0_xxzzzz_zz[k] = -g_z_0_x_0_xzzzz_zz[k] * ab_x + g_z_0_x_0_xzzzz_xzz[k];
            }

            /// Set up 1098-1104 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyyyy_xx = cbuffer.data(id_geom_1010_off + 1098 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyy_xy = cbuffer.data(id_geom_1010_off + 1099 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyy_xz = cbuffer.data(id_geom_1010_off + 1100 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyy_yy = cbuffer.data(id_geom_1010_off + 1101 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyy_yz = cbuffer.data(id_geom_1010_off + 1102 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyy_zz = cbuffer.data(id_geom_1010_off + 1103 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyyyy_xx, g_z_0_x_0_xyyyyy_xy, g_z_0_x_0_xyyyyy_xz, g_z_0_x_0_xyyyyy_yy, g_z_0_x_0_xyyyyy_yz, g_z_0_x_0_xyyyyy_zz, g_z_0_x_0_yyyyy_xx, g_z_0_x_0_yyyyy_xxx, g_z_0_x_0_yyyyy_xxy, g_z_0_x_0_yyyyy_xxz, g_z_0_x_0_yyyyy_xy, g_z_0_x_0_yyyyy_xyy, g_z_0_x_0_yyyyy_xyz, g_z_0_x_0_yyyyy_xz, g_z_0_x_0_yyyyy_xzz, g_z_0_x_0_yyyyy_yy, g_z_0_x_0_yyyyy_yz, g_z_0_x_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyyyy_xx[k] = -g_z_0_x_0_yyyyy_xx[k] * ab_x + g_z_0_x_0_yyyyy_xxx[k];

                g_z_0_x_0_xyyyyy_xy[k] = -g_z_0_x_0_yyyyy_xy[k] * ab_x + g_z_0_x_0_yyyyy_xxy[k];

                g_z_0_x_0_xyyyyy_xz[k] = -g_z_0_x_0_yyyyy_xz[k] * ab_x + g_z_0_x_0_yyyyy_xxz[k];

                g_z_0_x_0_xyyyyy_yy[k] = -g_z_0_x_0_yyyyy_yy[k] * ab_x + g_z_0_x_0_yyyyy_xyy[k];

                g_z_0_x_0_xyyyyy_yz[k] = -g_z_0_x_0_yyyyy_yz[k] * ab_x + g_z_0_x_0_yyyyy_xyz[k];

                g_z_0_x_0_xyyyyy_zz[k] = -g_z_0_x_0_yyyyy_zz[k] * ab_x + g_z_0_x_0_yyyyy_xzz[k];
            }

            /// Set up 1104-1110 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyyyz_xx = cbuffer.data(id_geom_1010_off + 1104 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyz_xy = cbuffer.data(id_geom_1010_off + 1105 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyz_xz = cbuffer.data(id_geom_1010_off + 1106 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyz_yy = cbuffer.data(id_geom_1010_off + 1107 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyz_yz = cbuffer.data(id_geom_1010_off + 1108 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyz_zz = cbuffer.data(id_geom_1010_off + 1109 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyyyz_xx, g_z_0_x_0_xyyyyz_xy, g_z_0_x_0_xyyyyz_xz, g_z_0_x_0_xyyyyz_yy, g_z_0_x_0_xyyyyz_yz, g_z_0_x_0_xyyyyz_zz, g_z_0_x_0_yyyyz_xx, g_z_0_x_0_yyyyz_xxx, g_z_0_x_0_yyyyz_xxy, g_z_0_x_0_yyyyz_xxz, g_z_0_x_0_yyyyz_xy, g_z_0_x_0_yyyyz_xyy, g_z_0_x_0_yyyyz_xyz, g_z_0_x_0_yyyyz_xz, g_z_0_x_0_yyyyz_xzz, g_z_0_x_0_yyyyz_yy, g_z_0_x_0_yyyyz_yz, g_z_0_x_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyyyz_xx[k] = -g_z_0_x_0_yyyyz_xx[k] * ab_x + g_z_0_x_0_yyyyz_xxx[k];

                g_z_0_x_0_xyyyyz_xy[k] = -g_z_0_x_0_yyyyz_xy[k] * ab_x + g_z_0_x_0_yyyyz_xxy[k];

                g_z_0_x_0_xyyyyz_xz[k] = -g_z_0_x_0_yyyyz_xz[k] * ab_x + g_z_0_x_0_yyyyz_xxz[k];

                g_z_0_x_0_xyyyyz_yy[k] = -g_z_0_x_0_yyyyz_yy[k] * ab_x + g_z_0_x_0_yyyyz_xyy[k];

                g_z_0_x_0_xyyyyz_yz[k] = -g_z_0_x_0_yyyyz_yz[k] * ab_x + g_z_0_x_0_yyyyz_xyz[k];

                g_z_0_x_0_xyyyyz_zz[k] = -g_z_0_x_0_yyyyz_zz[k] * ab_x + g_z_0_x_0_yyyyz_xzz[k];
            }

            /// Set up 1110-1116 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyyzz_xx = cbuffer.data(id_geom_1010_off + 1110 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyzz_xy = cbuffer.data(id_geom_1010_off + 1111 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyzz_xz = cbuffer.data(id_geom_1010_off + 1112 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyzz_yy = cbuffer.data(id_geom_1010_off + 1113 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyzz_yz = cbuffer.data(id_geom_1010_off + 1114 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyzz_zz = cbuffer.data(id_geom_1010_off + 1115 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyyzz_xx, g_z_0_x_0_xyyyzz_xy, g_z_0_x_0_xyyyzz_xz, g_z_0_x_0_xyyyzz_yy, g_z_0_x_0_xyyyzz_yz, g_z_0_x_0_xyyyzz_zz, g_z_0_x_0_yyyzz_xx, g_z_0_x_0_yyyzz_xxx, g_z_0_x_0_yyyzz_xxy, g_z_0_x_0_yyyzz_xxz, g_z_0_x_0_yyyzz_xy, g_z_0_x_0_yyyzz_xyy, g_z_0_x_0_yyyzz_xyz, g_z_0_x_0_yyyzz_xz, g_z_0_x_0_yyyzz_xzz, g_z_0_x_0_yyyzz_yy, g_z_0_x_0_yyyzz_yz, g_z_0_x_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyyzz_xx[k] = -g_z_0_x_0_yyyzz_xx[k] * ab_x + g_z_0_x_0_yyyzz_xxx[k];

                g_z_0_x_0_xyyyzz_xy[k] = -g_z_0_x_0_yyyzz_xy[k] * ab_x + g_z_0_x_0_yyyzz_xxy[k];

                g_z_0_x_0_xyyyzz_xz[k] = -g_z_0_x_0_yyyzz_xz[k] * ab_x + g_z_0_x_0_yyyzz_xxz[k];

                g_z_0_x_0_xyyyzz_yy[k] = -g_z_0_x_0_yyyzz_yy[k] * ab_x + g_z_0_x_0_yyyzz_xyy[k];

                g_z_0_x_0_xyyyzz_yz[k] = -g_z_0_x_0_yyyzz_yz[k] * ab_x + g_z_0_x_0_yyyzz_xyz[k];

                g_z_0_x_0_xyyyzz_zz[k] = -g_z_0_x_0_yyyzz_zz[k] * ab_x + g_z_0_x_0_yyyzz_xzz[k];
            }

            /// Set up 1116-1122 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyzzz_xx = cbuffer.data(id_geom_1010_off + 1116 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzzz_xy = cbuffer.data(id_geom_1010_off + 1117 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzzz_xz = cbuffer.data(id_geom_1010_off + 1118 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzzz_yy = cbuffer.data(id_geom_1010_off + 1119 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzzz_yz = cbuffer.data(id_geom_1010_off + 1120 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzzz_zz = cbuffer.data(id_geom_1010_off + 1121 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyzzz_xx, g_z_0_x_0_xyyzzz_xy, g_z_0_x_0_xyyzzz_xz, g_z_0_x_0_xyyzzz_yy, g_z_0_x_0_xyyzzz_yz, g_z_0_x_0_xyyzzz_zz, g_z_0_x_0_yyzzz_xx, g_z_0_x_0_yyzzz_xxx, g_z_0_x_0_yyzzz_xxy, g_z_0_x_0_yyzzz_xxz, g_z_0_x_0_yyzzz_xy, g_z_0_x_0_yyzzz_xyy, g_z_0_x_0_yyzzz_xyz, g_z_0_x_0_yyzzz_xz, g_z_0_x_0_yyzzz_xzz, g_z_0_x_0_yyzzz_yy, g_z_0_x_0_yyzzz_yz, g_z_0_x_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyzzz_xx[k] = -g_z_0_x_0_yyzzz_xx[k] * ab_x + g_z_0_x_0_yyzzz_xxx[k];

                g_z_0_x_0_xyyzzz_xy[k] = -g_z_0_x_0_yyzzz_xy[k] * ab_x + g_z_0_x_0_yyzzz_xxy[k];

                g_z_0_x_0_xyyzzz_xz[k] = -g_z_0_x_0_yyzzz_xz[k] * ab_x + g_z_0_x_0_yyzzz_xxz[k];

                g_z_0_x_0_xyyzzz_yy[k] = -g_z_0_x_0_yyzzz_yy[k] * ab_x + g_z_0_x_0_yyzzz_xyy[k];

                g_z_0_x_0_xyyzzz_yz[k] = -g_z_0_x_0_yyzzz_yz[k] * ab_x + g_z_0_x_0_yyzzz_xyz[k];

                g_z_0_x_0_xyyzzz_zz[k] = -g_z_0_x_0_yyzzz_zz[k] * ab_x + g_z_0_x_0_yyzzz_xzz[k];
            }

            /// Set up 1122-1128 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyzzzz_xx = cbuffer.data(id_geom_1010_off + 1122 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzzz_xy = cbuffer.data(id_geom_1010_off + 1123 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzzz_xz = cbuffer.data(id_geom_1010_off + 1124 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzzz_yy = cbuffer.data(id_geom_1010_off + 1125 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzzz_yz = cbuffer.data(id_geom_1010_off + 1126 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzzz_zz = cbuffer.data(id_geom_1010_off + 1127 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyzzzz_xx, g_z_0_x_0_xyzzzz_xy, g_z_0_x_0_xyzzzz_xz, g_z_0_x_0_xyzzzz_yy, g_z_0_x_0_xyzzzz_yz, g_z_0_x_0_xyzzzz_zz, g_z_0_x_0_yzzzz_xx, g_z_0_x_0_yzzzz_xxx, g_z_0_x_0_yzzzz_xxy, g_z_0_x_0_yzzzz_xxz, g_z_0_x_0_yzzzz_xy, g_z_0_x_0_yzzzz_xyy, g_z_0_x_0_yzzzz_xyz, g_z_0_x_0_yzzzz_xz, g_z_0_x_0_yzzzz_xzz, g_z_0_x_0_yzzzz_yy, g_z_0_x_0_yzzzz_yz, g_z_0_x_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyzzzz_xx[k] = -g_z_0_x_0_yzzzz_xx[k] * ab_x + g_z_0_x_0_yzzzz_xxx[k];

                g_z_0_x_0_xyzzzz_xy[k] = -g_z_0_x_0_yzzzz_xy[k] * ab_x + g_z_0_x_0_yzzzz_xxy[k];

                g_z_0_x_0_xyzzzz_xz[k] = -g_z_0_x_0_yzzzz_xz[k] * ab_x + g_z_0_x_0_yzzzz_xxz[k];

                g_z_0_x_0_xyzzzz_yy[k] = -g_z_0_x_0_yzzzz_yy[k] * ab_x + g_z_0_x_0_yzzzz_xyy[k];

                g_z_0_x_0_xyzzzz_yz[k] = -g_z_0_x_0_yzzzz_yz[k] * ab_x + g_z_0_x_0_yzzzz_xyz[k];

                g_z_0_x_0_xyzzzz_zz[k] = -g_z_0_x_0_yzzzz_zz[k] * ab_x + g_z_0_x_0_yzzzz_xzz[k];
            }

            /// Set up 1128-1134 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xzzzzz_xx = cbuffer.data(id_geom_1010_off + 1128 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzzz_xy = cbuffer.data(id_geom_1010_off + 1129 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzzz_xz = cbuffer.data(id_geom_1010_off + 1130 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzzz_yy = cbuffer.data(id_geom_1010_off + 1131 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzzz_yz = cbuffer.data(id_geom_1010_off + 1132 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzzz_zz = cbuffer.data(id_geom_1010_off + 1133 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xzzzzz_xx, g_z_0_x_0_xzzzzz_xy, g_z_0_x_0_xzzzzz_xz, g_z_0_x_0_xzzzzz_yy, g_z_0_x_0_xzzzzz_yz, g_z_0_x_0_xzzzzz_zz, g_z_0_x_0_zzzzz_xx, g_z_0_x_0_zzzzz_xxx, g_z_0_x_0_zzzzz_xxy, g_z_0_x_0_zzzzz_xxz, g_z_0_x_0_zzzzz_xy, g_z_0_x_0_zzzzz_xyy, g_z_0_x_0_zzzzz_xyz, g_z_0_x_0_zzzzz_xz, g_z_0_x_0_zzzzz_xzz, g_z_0_x_0_zzzzz_yy, g_z_0_x_0_zzzzz_yz, g_z_0_x_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xzzzzz_xx[k] = -g_z_0_x_0_zzzzz_xx[k] * ab_x + g_z_0_x_0_zzzzz_xxx[k];

                g_z_0_x_0_xzzzzz_xy[k] = -g_z_0_x_0_zzzzz_xy[k] * ab_x + g_z_0_x_0_zzzzz_xxy[k];

                g_z_0_x_0_xzzzzz_xz[k] = -g_z_0_x_0_zzzzz_xz[k] * ab_x + g_z_0_x_0_zzzzz_xxz[k];

                g_z_0_x_0_xzzzzz_yy[k] = -g_z_0_x_0_zzzzz_yy[k] * ab_x + g_z_0_x_0_zzzzz_xyy[k];

                g_z_0_x_0_xzzzzz_yz[k] = -g_z_0_x_0_zzzzz_yz[k] * ab_x + g_z_0_x_0_zzzzz_xyz[k];

                g_z_0_x_0_xzzzzz_zz[k] = -g_z_0_x_0_zzzzz_zz[k] * ab_x + g_z_0_x_0_zzzzz_xzz[k];
            }

            /// Set up 1134-1140 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyyyy_xx = cbuffer.data(id_geom_1010_off + 1134 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyy_xy = cbuffer.data(id_geom_1010_off + 1135 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyy_xz = cbuffer.data(id_geom_1010_off + 1136 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyy_yy = cbuffer.data(id_geom_1010_off + 1137 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyy_yz = cbuffer.data(id_geom_1010_off + 1138 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyy_zz = cbuffer.data(id_geom_1010_off + 1139 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyyy_xx, g_z_0_x_0_yyyyy_xxy, g_z_0_x_0_yyyyy_xy, g_z_0_x_0_yyyyy_xyy, g_z_0_x_0_yyyyy_xyz, g_z_0_x_0_yyyyy_xz, g_z_0_x_0_yyyyy_yy, g_z_0_x_0_yyyyy_yyy, g_z_0_x_0_yyyyy_yyz, g_z_0_x_0_yyyyy_yz, g_z_0_x_0_yyyyy_yzz, g_z_0_x_0_yyyyy_zz, g_z_0_x_0_yyyyyy_xx, g_z_0_x_0_yyyyyy_xy, g_z_0_x_0_yyyyyy_xz, g_z_0_x_0_yyyyyy_yy, g_z_0_x_0_yyyyyy_yz, g_z_0_x_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyyyy_xx[k] = -g_z_0_x_0_yyyyy_xx[k] * ab_y + g_z_0_x_0_yyyyy_xxy[k];

                g_z_0_x_0_yyyyyy_xy[k] = -g_z_0_x_0_yyyyy_xy[k] * ab_y + g_z_0_x_0_yyyyy_xyy[k];

                g_z_0_x_0_yyyyyy_xz[k] = -g_z_0_x_0_yyyyy_xz[k] * ab_y + g_z_0_x_0_yyyyy_xyz[k];

                g_z_0_x_0_yyyyyy_yy[k] = -g_z_0_x_0_yyyyy_yy[k] * ab_y + g_z_0_x_0_yyyyy_yyy[k];

                g_z_0_x_0_yyyyyy_yz[k] = -g_z_0_x_0_yyyyy_yz[k] * ab_y + g_z_0_x_0_yyyyy_yyz[k];

                g_z_0_x_0_yyyyyy_zz[k] = -g_z_0_x_0_yyyyy_zz[k] * ab_y + g_z_0_x_0_yyyyy_yzz[k];
            }

            /// Set up 1140-1146 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyyyz_xx = cbuffer.data(id_geom_1010_off + 1140 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyz_xy = cbuffer.data(id_geom_1010_off + 1141 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyz_xz = cbuffer.data(id_geom_1010_off + 1142 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyz_yy = cbuffer.data(id_geom_1010_off + 1143 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyz_yz = cbuffer.data(id_geom_1010_off + 1144 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyz_zz = cbuffer.data(id_geom_1010_off + 1145 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyyyz_xx, g_z_0_x_0_yyyyyz_xy, g_z_0_x_0_yyyyyz_xz, g_z_0_x_0_yyyyyz_yy, g_z_0_x_0_yyyyyz_yz, g_z_0_x_0_yyyyyz_zz, g_z_0_x_0_yyyyz_xx, g_z_0_x_0_yyyyz_xxy, g_z_0_x_0_yyyyz_xy, g_z_0_x_0_yyyyz_xyy, g_z_0_x_0_yyyyz_xyz, g_z_0_x_0_yyyyz_xz, g_z_0_x_0_yyyyz_yy, g_z_0_x_0_yyyyz_yyy, g_z_0_x_0_yyyyz_yyz, g_z_0_x_0_yyyyz_yz, g_z_0_x_0_yyyyz_yzz, g_z_0_x_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyyyz_xx[k] = -g_z_0_x_0_yyyyz_xx[k] * ab_y + g_z_0_x_0_yyyyz_xxy[k];

                g_z_0_x_0_yyyyyz_xy[k] = -g_z_0_x_0_yyyyz_xy[k] * ab_y + g_z_0_x_0_yyyyz_xyy[k];

                g_z_0_x_0_yyyyyz_xz[k] = -g_z_0_x_0_yyyyz_xz[k] * ab_y + g_z_0_x_0_yyyyz_xyz[k];

                g_z_0_x_0_yyyyyz_yy[k] = -g_z_0_x_0_yyyyz_yy[k] * ab_y + g_z_0_x_0_yyyyz_yyy[k];

                g_z_0_x_0_yyyyyz_yz[k] = -g_z_0_x_0_yyyyz_yz[k] * ab_y + g_z_0_x_0_yyyyz_yyz[k];

                g_z_0_x_0_yyyyyz_zz[k] = -g_z_0_x_0_yyyyz_zz[k] * ab_y + g_z_0_x_0_yyyyz_yzz[k];
            }

            /// Set up 1146-1152 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyyzz_xx = cbuffer.data(id_geom_1010_off + 1146 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyzz_xy = cbuffer.data(id_geom_1010_off + 1147 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyzz_xz = cbuffer.data(id_geom_1010_off + 1148 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyzz_yy = cbuffer.data(id_geom_1010_off + 1149 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyzz_yz = cbuffer.data(id_geom_1010_off + 1150 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyzz_zz = cbuffer.data(id_geom_1010_off + 1151 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyyzz_xx, g_z_0_x_0_yyyyzz_xy, g_z_0_x_0_yyyyzz_xz, g_z_0_x_0_yyyyzz_yy, g_z_0_x_0_yyyyzz_yz, g_z_0_x_0_yyyyzz_zz, g_z_0_x_0_yyyzz_xx, g_z_0_x_0_yyyzz_xxy, g_z_0_x_0_yyyzz_xy, g_z_0_x_0_yyyzz_xyy, g_z_0_x_0_yyyzz_xyz, g_z_0_x_0_yyyzz_xz, g_z_0_x_0_yyyzz_yy, g_z_0_x_0_yyyzz_yyy, g_z_0_x_0_yyyzz_yyz, g_z_0_x_0_yyyzz_yz, g_z_0_x_0_yyyzz_yzz, g_z_0_x_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyyzz_xx[k] = -g_z_0_x_0_yyyzz_xx[k] * ab_y + g_z_0_x_0_yyyzz_xxy[k];

                g_z_0_x_0_yyyyzz_xy[k] = -g_z_0_x_0_yyyzz_xy[k] * ab_y + g_z_0_x_0_yyyzz_xyy[k];

                g_z_0_x_0_yyyyzz_xz[k] = -g_z_0_x_0_yyyzz_xz[k] * ab_y + g_z_0_x_0_yyyzz_xyz[k];

                g_z_0_x_0_yyyyzz_yy[k] = -g_z_0_x_0_yyyzz_yy[k] * ab_y + g_z_0_x_0_yyyzz_yyy[k];

                g_z_0_x_0_yyyyzz_yz[k] = -g_z_0_x_0_yyyzz_yz[k] * ab_y + g_z_0_x_0_yyyzz_yyz[k];

                g_z_0_x_0_yyyyzz_zz[k] = -g_z_0_x_0_yyyzz_zz[k] * ab_y + g_z_0_x_0_yyyzz_yzz[k];
            }

            /// Set up 1152-1158 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyzzz_xx = cbuffer.data(id_geom_1010_off + 1152 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzzz_xy = cbuffer.data(id_geom_1010_off + 1153 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzzz_xz = cbuffer.data(id_geom_1010_off + 1154 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzzz_yy = cbuffer.data(id_geom_1010_off + 1155 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzzz_yz = cbuffer.data(id_geom_1010_off + 1156 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzzz_zz = cbuffer.data(id_geom_1010_off + 1157 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyzzz_xx, g_z_0_x_0_yyyzzz_xy, g_z_0_x_0_yyyzzz_xz, g_z_0_x_0_yyyzzz_yy, g_z_0_x_0_yyyzzz_yz, g_z_0_x_0_yyyzzz_zz, g_z_0_x_0_yyzzz_xx, g_z_0_x_0_yyzzz_xxy, g_z_0_x_0_yyzzz_xy, g_z_0_x_0_yyzzz_xyy, g_z_0_x_0_yyzzz_xyz, g_z_0_x_0_yyzzz_xz, g_z_0_x_0_yyzzz_yy, g_z_0_x_0_yyzzz_yyy, g_z_0_x_0_yyzzz_yyz, g_z_0_x_0_yyzzz_yz, g_z_0_x_0_yyzzz_yzz, g_z_0_x_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyzzz_xx[k] = -g_z_0_x_0_yyzzz_xx[k] * ab_y + g_z_0_x_0_yyzzz_xxy[k];

                g_z_0_x_0_yyyzzz_xy[k] = -g_z_0_x_0_yyzzz_xy[k] * ab_y + g_z_0_x_0_yyzzz_xyy[k];

                g_z_0_x_0_yyyzzz_xz[k] = -g_z_0_x_0_yyzzz_xz[k] * ab_y + g_z_0_x_0_yyzzz_xyz[k];

                g_z_0_x_0_yyyzzz_yy[k] = -g_z_0_x_0_yyzzz_yy[k] * ab_y + g_z_0_x_0_yyzzz_yyy[k];

                g_z_0_x_0_yyyzzz_yz[k] = -g_z_0_x_0_yyzzz_yz[k] * ab_y + g_z_0_x_0_yyzzz_yyz[k];

                g_z_0_x_0_yyyzzz_zz[k] = -g_z_0_x_0_yyzzz_zz[k] * ab_y + g_z_0_x_0_yyzzz_yzz[k];
            }

            /// Set up 1158-1164 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyzzzz_xx = cbuffer.data(id_geom_1010_off + 1158 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzzz_xy = cbuffer.data(id_geom_1010_off + 1159 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzzz_xz = cbuffer.data(id_geom_1010_off + 1160 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzzz_yy = cbuffer.data(id_geom_1010_off + 1161 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzzz_yz = cbuffer.data(id_geom_1010_off + 1162 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzzz_zz = cbuffer.data(id_geom_1010_off + 1163 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyzzzz_xx, g_z_0_x_0_yyzzzz_xy, g_z_0_x_0_yyzzzz_xz, g_z_0_x_0_yyzzzz_yy, g_z_0_x_0_yyzzzz_yz, g_z_0_x_0_yyzzzz_zz, g_z_0_x_0_yzzzz_xx, g_z_0_x_0_yzzzz_xxy, g_z_0_x_0_yzzzz_xy, g_z_0_x_0_yzzzz_xyy, g_z_0_x_0_yzzzz_xyz, g_z_0_x_0_yzzzz_xz, g_z_0_x_0_yzzzz_yy, g_z_0_x_0_yzzzz_yyy, g_z_0_x_0_yzzzz_yyz, g_z_0_x_0_yzzzz_yz, g_z_0_x_0_yzzzz_yzz, g_z_0_x_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyzzzz_xx[k] = -g_z_0_x_0_yzzzz_xx[k] * ab_y + g_z_0_x_0_yzzzz_xxy[k];

                g_z_0_x_0_yyzzzz_xy[k] = -g_z_0_x_0_yzzzz_xy[k] * ab_y + g_z_0_x_0_yzzzz_xyy[k];

                g_z_0_x_0_yyzzzz_xz[k] = -g_z_0_x_0_yzzzz_xz[k] * ab_y + g_z_0_x_0_yzzzz_xyz[k];

                g_z_0_x_0_yyzzzz_yy[k] = -g_z_0_x_0_yzzzz_yy[k] * ab_y + g_z_0_x_0_yzzzz_yyy[k];

                g_z_0_x_0_yyzzzz_yz[k] = -g_z_0_x_0_yzzzz_yz[k] * ab_y + g_z_0_x_0_yzzzz_yyz[k];

                g_z_0_x_0_yyzzzz_zz[k] = -g_z_0_x_0_yzzzz_zz[k] * ab_y + g_z_0_x_0_yzzzz_yzz[k];
            }

            /// Set up 1164-1170 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yzzzzz_xx = cbuffer.data(id_geom_1010_off + 1164 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzzz_xy = cbuffer.data(id_geom_1010_off + 1165 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzzz_xz = cbuffer.data(id_geom_1010_off + 1166 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzzz_yy = cbuffer.data(id_geom_1010_off + 1167 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzzz_yz = cbuffer.data(id_geom_1010_off + 1168 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzzz_zz = cbuffer.data(id_geom_1010_off + 1169 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yzzzzz_xx, g_z_0_x_0_yzzzzz_xy, g_z_0_x_0_yzzzzz_xz, g_z_0_x_0_yzzzzz_yy, g_z_0_x_0_yzzzzz_yz, g_z_0_x_0_yzzzzz_zz, g_z_0_x_0_zzzzz_xx, g_z_0_x_0_zzzzz_xxy, g_z_0_x_0_zzzzz_xy, g_z_0_x_0_zzzzz_xyy, g_z_0_x_0_zzzzz_xyz, g_z_0_x_0_zzzzz_xz, g_z_0_x_0_zzzzz_yy, g_z_0_x_0_zzzzz_yyy, g_z_0_x_0_zzzzz_yyz, g_z_0_x_0_zzzzz_yz, g_z_0_x_0_zzzzz_yzz, g_z_0_x_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yzzzzz_xx[k] = -g_z_0_x_0_zzzzz_xx[k] * ab_y + g_z_0_x_0_zzzzz_xxy[k];

                g_z_0_x_0_yzzzzz_xy[k] = -g_z_0_x_0_zzzzz_xy[k] * ab_y + g_z_0_x_0_zzzzz_xyy[k];

                g_z_0_x_0_yzzzzz_xz[k] = -g_z_0_x_0_zzzzz_xz[k] * ab_y + g_z_0_x_0_zzzzz_xyz[k];

                g_z_0_x_0_yzzzzz_yy[k] = -g_z_0_x_0_zzzzz_yy[k] * ab_y + g_z_0_x_0_zzzzz_yyy[k];

                g_z_0_x_0_yzzzzz_yz[k] = -g_z_0_x_0_zzzzz_yz[k] * ab_y + g_z_0_x_0_zzzzz_yyz[k];

                g_z_0_x_0_yzzzzz_zz[k] = -g_z_0_x_0_zzzzz_zz[k] * ab_y + g_z_0_x_0_zzzzz_yzz[k];
            }

            /// Set up 1170-1176 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_zzzzzz_xx = cbuffer.data(id_geom_1010_off + 1170 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzzz_xy = cbuffer.data(id_geom_1010_off + 1171 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzzz_xz = cbuffer.data(id_geom_1010_off + 1172 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzzz_yy = cbuffer.data(id_geom_1010_off + 1173 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzzz_yz = cbuffer.data(id_geom_1010_off + 1174 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzzz_zz = cbuffer.data(id_geom_1010_off + 1175 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_zzzzz_xx, g_0_0_x_0_zzzzz_xy, g_0_0_x_0_zzzzz_xz, g_0_0_x_0_zzzzz_yy, g_0_0_x_0_zzzzz_yz, g_0_0_x_0_zzzzz_zz, g_z_0_x_0_zzzzz_xx, g_z_0_x_0_zzzzz_xxz, g_z_0_x_0_zzzzz_xy, g_z_0_x_0_zzzzz_xyz, g_z_0_x_0_zzzzz_xz, g_z_0_x_0_zzzzz_xzz, g_z_0_x_0_zzzzz_yy, g_z_0_x_0_zzzzz_yyz, g_z_0_x_0_zzzzz_yz, g_z_0_x_0_zzzzz_yzz, g_z_0_x_0_zzzzz_zz, g_z_0_x_0_zzzzz_zzz, g_z_0_x_0_zzzzzz_xx, g_z_0_x_0_zzzzzz_xy, g_z_0_x_0_zzzzzz_xz, g_z_0_x_0_zzzzzz_yy, g_z_0_x_0_zzzzzz_yz, g_z_0_x_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_zzzzzz_xx[k] = -g_0_0_x_0_zzzzz_xx[k] - g_z_0_x_0_zzzzz_xx[k] * ab_z + g_z_0_x_0_zzzzz_xxz[k];

                g_z_0_x_0_zzzzzz_xy[k] = -g_0_0_x_0_zzzzz_xy[k] - g_z_0_x_0_zzzzz_xy[k] * ab_z + g_z_0_x_0_zzzzz_xyz[k];

                g_z_0_x_0_zzzzzz_xz[k] = -g_0_0_x_0_zzzzz_xz[k] - g_z_0_x_0_zzzzz_xz[k] * ab_z + g_z_0_x_0_zzzzz_xzz[k];

                g_z_0_x_0_zzzzzz_yy[k] = -g_0_0_x_0_zzzzz_yy[k] - g_z_0_x_0_zzzzz_yy[k] * ab_z + g_z_0_x_0_zzzzz_yyz[k];

                g_z_0_x_0_zzzzzz_yz[k] = -g_0_0_x_0_zzzzz_yz[k] - g_z_0_x_0_zzzzz_yz[k] * ab_z + g_z_0_x_0_zzzzz_yzz[k];

                g_z_0_x_0_zzzzzz_zz[k] = -g_0_0_x_0_zzzzz_zz[k] - g_z_0_x_0_zzzzz_zz[k] * ab_z + g_z_0_x_0_zzzzz_zzz[k];
            }

            /// Set up 1176-1182 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxxx_xx = cbuffer.data(id_geom_1010_off + 1176 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxx_xy = cbuffer.data(id_geom_1010_off + 1177 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxx_xz = cbuffer.data(id_geom_1010_off + 1178 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxx_yy = cbuffer.data(id_geom_1010_off + 1179 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxx_yz = cbuffer.data(id_geom_1010_off + 1180 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxx_zz = cbuffer.data(id_geom_1010_off + 1181 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxx_xx, g_z_0_y_0_xxxxx_xxx, g_z_0_y_0_xxxxx_xxy, g_z_0_y_0_xxxxx_xxz, g_z_0_y_0_xxxxx_xy, g_z_0_y_0_xxxxx_xyy, g_z_0_y_0_xxxxx_xyz, g_z_0_y_0_xxxxx_xz, g_z_0_y_0_xxxxx_xzz, g_z_0_y_0_xxxxx_yy, g_z_0_y_0_xxxxx_yz, g_z_0_y_0_xxxxx_zz, g_z_0_y_0_xxxxxx_xx, g_z_0_y_0_xxxxxx_xy, g_z_0_y_0_xxxxxx_xz, g_z_0_y_0_xxxxxx_yy, g_z_0_y_0_xxxxxx_yz, g_z_0_y_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxxx_xx[k] = -g_z_0_y_0_xxxxx_xx[k] * ab_x + g_z_0_y_0_xxxxx_xxx[k];

                g_z_0_y_0_xxxxxx_xy[k] = -g_z_0_y_0_xxxxx_xy[k] * ab_x + g_z_0_y_0_xxxxx_xxy[k];

                g_z_0_y_0_xxxxxx_xz[k] = -g_z_0_y_0_xxxxx_xz[k] * ab_x + g_z_0_y_0_xxxxx_xxz[k];

                g_z_0_y_0_xxxxxx_yy[k] = -g_z_0_y_0_xxxxx_yy[k] * ab_x + g_z_0_y_0_xxxxx_xyy[k];

                g_z_0_y_0_xxxxxx_yz[k] = -g_z_0_y_0_xxxxx_yz[k] * ab_x + g_z_0_y_0_xxxxx_xyz[k];

                g_z_0_y_0_xxxxxx_zz[k] = -g_z_0_y_0_xxxxx_zz[k] * ab_x + g_z_0_y_0_xxxxx_xzz[k];
            }

            /// Set up 1182-1188 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxxy_xx = cbuffer.data(id_geom_1010_off + 1182 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxy_xy = cbuffer.data(id_geom_1010_off + 1183 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxy_xz = cbuffer.data(id_geom_1010_off + 1184 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxy_yy = cbuffer.data(id_geom_1010_off + 1185 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxy_yz = cbuffer.data(id_geom_1010_off + 1186 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxy_zz = cbuffer.data(id_geom_1010_off + 1187 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxxy_xx, g_z_0_y_0_xxxxxy_xy, g_z_0_y_0_xxxxxy_xz, g_z_0_y_0_xxxxxy_yy, g_z_0_y_0_xxxxxy_yz, g_z_0_y_0_xxxxxy_zz, g_z_0_y_0_xxxxy_xx, g_z_0_y_0_xxxxy_xxx, g_z_0_y_0_xxxxy_xxy, g_z_0_y_0_xxxxy_xxz, g_z_0_y_0_xxxxy_xy, g_z_0_y_0_xxxxy_xyy, g_z_0_y_0_xxxxy_xyz, g_z_0_y_0_xxxxy_xz, g_z_0_y_0_xxxxy_xzz, g_z_0_y_0_xxxxy_yy, g_z_0_y_0_xxxxy_yz, g_z_0_y_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxxy_xx[k] = -g_z_0_y_0_xxxxy_xx[k] * ab_x + g_z_0_y_0_xxxxy_xxx[k];

                g_z_0_y_0_xxxxxy_xy[k] = -g_z_0_y_0_xxxxy_xy[k] * ab_x + g_z_0_y_0_xxxxy_xxy[k];

                g_z_0_y_0_xxxxxy_xz[k] = -g_z_0_y_0_xxxxy_xz[k] * ab_x + g_z_0_y_0_xxxxy_xxz[k];

                g_z_0_y_0_xxxxxy_yy[k] = -g_z_0_y_0_xxxxy_yy[k] * ab_x + g_z_0_y_0_xxxxy_xyy[k];

                g_z_0_y_0_xxxxxy_yz[k] = -g_z_0_y_0_xxxxy_yz[k] * ab_x + g_z_0_y_0_xxxxy_xyz[k];

                g_z_0_y_0_xxxxxy_zz[k] = -g_z_0_y_0_xxxxy_zz[k] * ab_x + g_z_0_y_0_xxxxy_xzz[k];
            }

            /// Set up 1188-1194 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxxz_xx = cbuffer.data(id_geom_1010_off + 1188 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxz_xy = cbuffer.data(id_geom_1010_off + 1189 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxz_xz = cbuffer.data(id_geom_1010_off + 1190 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxz_yy = cbuffer.data(id_geom_1010_off + 1191 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxz_yz = cbuffer.data(id_geom_1010_off + 1192 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxz_zz = cbuffer.data(id_geom_1010_off + 1193 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxxz_xx, g_z_0_y_0_xxxxxz_xy, g_z_0_y_0_xxxxxz_xz, g_z_0_y_0_xxxxxz_yy, g_z_0_y_0_xxxxxz_yz, g_z_0_y_0_xxxxxz_zz, g_z_0_y_0_xxxxz_xx, g_z_0_y_0_xxxxz_xxx, g_z_0_y_0_xxxxz_xxy, g_z_0_y_0_xxxxz_xxz, g_z_0_y_0_xxxxz_xy, g_z_0_y_0_xxxxz_xyy, g_z_0_y_0_xxxxz_xyz, g_z_0_y_0_xxxxz_xz, g_z_0_y_0_xxxxz_xzz, g_z_0_y_0_xxxxz_yy, g_z_0_y_0_xxxxz_yz, g_z_0_y_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxxz_xx[k] = -g_z_0_y_0_xxxxz_xx[k] * ab_x + g_z_0_y_0_xxxxz_xxx[k];

                g_z_0_y_0_xxxxxz_xy[k] = -g_z_0_y_0_xxxxz_xy[k] * ab_x + g_z_0_y_0_xxxxz_xxy[k];

                g_z_0_y_0_xxxxxz_xz[k] = -g_z_0_y_0_xxxxz_xz[k] * ab_x + g_z_0_y_0_xxxxz_xxz[k];

                g_z_0_y_0_xxxxxz_yy[k] = -g_z_0_y_0_xxxxz_yy[k] * ab_x + g_z_0_y_0_xxxxz_xyy[k];

                g_z_0_y_0_xxxxxz_yz[k] = -g_z_0_y_0_xxxxz_yz[k] * ab_x + g_z_0_y_0_xxxxz_xyz[k];

                g_z_0_y_0_xxxxxz_zz[k] = -g_z_0_y_0_xxxxz_zz[k] * ab_x + g_z_0_y_0_xxxxz_xzz[k];
            }

            /// Set up 1194-1200 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxyy_xx = cbuffer.data(id_geom_1010_off + 1194 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyy_xy = cbuffer.data(id_geom_1010_off + 1195 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyy_xz = cbuffer.data(id_geom_1010_off + 1196 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyy_yy = cbuffer.data(id_geom_1010_off + 1197 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyy_yz = cbuffer.data(id_geom_1010_off + 1198 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyy_zz = cbuffer.data(id_geom_1010_off + 1199 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxyy_xx, g_z_0_y_0_xxxxyy_xy, g_z_0_y_0_xxxxyy_xz, g_z_0_y_0_xxxxyy_yy, g_z_0_y_0_xxxxyy_yz, g_z_0_y_0_xxxxyy_zz, g_z_0_y_0_xxxyy_xx, g_z_0_y_0_xxxyy_xxx, g_z_0_y_0_xxxyy_xxy, g_z_0_y_0_xxxyy_xxz, g_z_0_y_0_xxxyy_xy, g_z_0_y_0_xxxyy_xyy, g_z_0_y_0_xxxyy_xyz, g_z_0_y_0_xxxyy_xz, g_z_0_y_0_xxxyy_xzz, g_z_0_y_0_xxxyy_yy, g_z_0_y_0_xxxyy_yz, g_z_0_y_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxyy_xx[k] = -g_z_0_y_0_xxxyy_xx[k] * ab_x + g_z_0_y_0_xxxyy_xxx[k];

                g_z_0_y_0_xxxxyy_xy[k] = -g_z_0_y_0_xxxyy_xy[k] * ab_x + g_z_0_y_0_xxxyy_xxy[k];

                g_z_0_y_0_xxxxyy_xz[k] = -g_z_0_y_0_xxxyy_xz[k] * ab_x + g_z_0_y_0_xxxyy_xxz[k];

                g_z_0_y_0_xxxxyy_yy[k] = -g_z_0_y_0_xxxyy_yy[k] * ab_x + g_z_0_y_0_xxxyy_xyy[k];

                g_z_0_y_0_xxxxyy_yz[k] = -g_z_0_y_0_xxxyy_yz[k] * ab_x + g_z_0_y_0_xxxyy_xyz[k];

                g_z_0_y_0_xxxxyy_zz[k] = -g_z_0_y_0_xxxyy_zz[k] * ab_x + g_z_0_y_0_xxxyy_xzz[k];
            }

            /// Set up 1200-1206 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxyz_xx = cbuffer.data(id_geom_1010_off + 1200 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyz_xy = cbuffer.data(id_geom_1010_off + 1201 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyz_xz = cbuffer.data(id_geom_1010_off + 1202 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyz_yy = cbuffer.data(id_geom_1010_off + 1203 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyz_yz = cbuffer.data(id_geom_1010_off + 1204 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyz_zz = cbuffer.data(id_geom_1010_off + 1205 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxyz_xx, g_z_0_y_0_xxxxyz_xy, g_z_0_y_0_xxxxyz_xz, g_z_0_y_0_xxxxyz_yy, g_z_0_y_0_xxxxyz_yz, g_z_0_y_0_xxxxyz_zz, g_z_0_y_0_xxxyz_xx, g_z_0_y_0_xxxyz_xxx, g_z_0_y_0_xxxyz_xxy, g_z_0_y_0_xxxyz_xxz, g_z_0_y_0_xxxyz_xy, g_z_0_y_0_xxxyz_xyy, g_z_0_y_0_xxxyz_xyz, g_z_0_y_0_xxxyz_xz, g_z_0_y_0_xxxyz_xzz, g_z_0_y_0_xxxyz_yy, g_z_0_y_0_xxxyz_yz, g_z_0_y_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxyz_xx[k] = -g_z_0_y_0_xxxyz_xx[k] * ab_x + g_z_0_y_0_xxxyz_xxx[k];

                g_z_0_y_0_xxxxyz_xy[k] = -g_z_0_y_0_xxxyz_xy[k] * ab_x + g_z_0_y_0_xxxyz_xxy[k];

                g_z_0_y_0_xxxxyz_xz[k] = -g_z_0_y_0_xxxyz_xz[k] * ab_x + g_z_0_y_0_xxxyz_xxz[k];

                g_z_0_y_0_xxxxyz_yy[k] = -g_z_0_y_0_xxxyz_yy[k] * ab_x + g_z_0_y_0_xxxyz_xyy[k];

                g_z_0_y_0_xxxxyz_yz[k] = -g_z_0_y_0_xxxyz_yz[k] * ab_x + g_z_0_y_0_xxxyz_xyz[k];

                g_z_0_y_0_xxxxyz_zz[k] = -g_z_0_y_0_xxxyz_zz[k] * ab_x + g_z_0_y_0_xxxyz_xzz[k];
            }

            /// Set up 1206-1212 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxzz_xx = cbuffer.data(id_geom_1010_off + 1206 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxzz_xy = cbuffer.data(id_geom_1010_off + 1207 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxzz_xz = cbuffer.data(id_geom_1010_off + 1208 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxzz_yy = cbuffer.data(id_geom_1010_off + 1209 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxzz_yz = cbuffer.data(id_geom_1010_off + 1210 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxzz_zz = cbuffer.data(id_geom_1010_off + 1211 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxzz_xx, g_z_0_y_0_xxxxzz_xy, g_z_0_y_0_xxxxzz_xz, g_z_0_y_0_xxxxzz_yy, g_z_0_y_0_xxxxzz_yz, g_z_0_y_0_xxxxzz_zz, g_z_0_y_0_xxxzz_xx, g_z_0_y_0_xxxzz_xxx, g_z_0_y_0_xxxzz_xxy, g_z_0_y_0_xxxzz_xxz, g_z_0_y_0_xxxzz_xy, g_z_0_y_0_xxxzz_xyy, g_z_0_y_0_xxxzz_xyz, g_z_0_y_0_xxxzz_xz, g_z_0_y_0_xxxzz_xzz, g_z_0_y_0_xxxzz_yy, g_z_0_y_0_xxxzz_yz, g_z_0_y_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxzz_xx[k] = -g_z_0_y_0_xxxzz_xx[k] * ab_x + g_z_0_y_0_xxxzz_xxx[k];

                g_z_0_y_0_xxxxzz_xy[k] = -g_z_0_y_0_xxxzz_xy[k] * ab_x + g_z_0_y_0_xxxzz_xxy[k];

                g_z_0_y_0_xxxxzz_xz[k] = -g_z_0_y_0_xxxzz_xz[k] * ab_x + g_z_0_y_0_xxxzz_xxz[k];

                g_z_0_y_0_xxxxzz_yy[k] = -g_z_0_y_0_xxxzz_yy[k] * ab_x + g_z_0_y_0_xxxzz_xyy[k];

                g_z_0_y_0_xxxxzz_yz[k] = -g_z_0_y_0_xxxzz_yz[k] * ab_x + g_z_0_y_0_xxxzz_xyz[k];

                g_z_0_y_0_xxxxzz_zz[k] = -g_z_0_y_0_xxxzz_zz[k] * ab_x + g_z_0_y_0_xxxzz_xzz[k];
            }

            /// Set up 1212-1218 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxyyy_xx = cbuffer.data(id_geom_1010_off + 1212 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyy_xy = cbuffer.data(id_geom_1010_off + 1213 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyy_xz = cbuffer.data(id_geom_1010_off + 1214 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyy_yy = cbuffer.data(id_geom_1010_off + 1215 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyy_yz = cbuffer.data(id_geom_1010_off + 1216 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyy_zz = cbuffer.data(id_geom_1010_off + 1217 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxyyy_xx, g_z_0_y_0_xxxyyy_xy, g_z_0_y_0_xxxyyy_xz, g_z_0_y_0_xxxyyy_yy, g_z_0_y_0_xxxyyy_yz, g_z_0_y_0_xxxyyy_zz, g_z_0_y_0_xxyyy_xx, g_z_0_y_0_xxyyy_xxx, g_z_0_y_0_xxyyy_xxy, g_z_0_y_0_xxyyy_xxz, g_z_0_y_0_xxyyy_xy, g_z_0_y_0_xxyyy_xyy, g_z_0_y_0_xxyyy_xyz, g_z_0_y_0_xxyyy_xz, g_z_0_y_0_xxyyy_xzz, g_z_0_y_0_xxyyy_yy, g_z_0_y_0_xxyyy_yz, g_z_0_y_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxyyy_xx[k] = -g_z_0_y_0_xxyyy_xx[k] * ab_x + g_z_0_y_0_xxyyy_xxx[k];

                g_z_0_y_0_xxxyyy_xy[k] = -g_z_0_y_0_xxyyy_xy[k] * ab_x + g_z_0_y_0_xxyyy_xxy[k];

                g_z_0_y_0_xxxyyy_xz[k] = -g_z_0_y_0_xxyyy_xz[k] * ab_x + g_z_0_y_0_xxyyy_xxz[k];

                g_z_0_y_0_xxxyyy_yy[k] = -g_z_0_y_0_xxyyy_yy[k] * ab_x + g_z_0_y_0_xxyyy_xyy[k];

                g_z_0_y_0_xxxyyy_yz[k] = -g_z_0_y_0_xxyyy_yz[k] * ab_x + g_z_0_y_0_xxyyy_xyz[k];

                g_z_0_y_0_xxxyyy_zz[k] = -g_z_0_y_0_xxyyy_zz[k] * ab_x + g_z_0_y_0_xxyyy_xzz[k];
            }

            /// Set up 1218-1224 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxyyz_xx = cbuffer.data(id_geom_1010_off + 1218 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyz_xy = cbuffer.data(id_geom_1010_off + 1219 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyz_xz = cbuffer.data(id_geom_1010_off + 1220 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyz_yy = cbuffer.data(id_geom_1010_off + 1221 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyz_yz = cbuffer.data(id_geom_1010_off + 1222 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyz_zz = cbuffer.data(id_geom_1010_off + 1223 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxyyz_xx, g_z_0_y_0_xxxyyz_xy, g_z_0_y_0_xxxyyz_xz, g_z_0_y_0_xxxyyz_yy, g_z_0_y_0_xxxyyz_yz, g_z_0_y_0_xxxyyz_zz, g_z_0_y_0_xxyyz_xx, g_z_0_y_0_xxyyz_xxx, g_z_0_y_0_xxyyz_xxy, g_z_0_y_0_xxyyz_xxz, g_z_0_y_0_xxyyz_xy, g_z_0_y_0_xxyyz_xyy, g_z_0_y_0_xxyyz_xyz, g_z_0_y_0_xxyyz_xz, g_z_0_y_0_xxyyz_xzz, g_z_0_y_0_xxyyz_yy, g_z_0_y_0_xxyyz_yz, g_z_0_y_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxyyz_xx[k] = -g_z_0_y_0_xxyyz_xx[k] * ab_x + g_z_0_y_0_xxyyz_xxx[k];

                g_z_0_y_0_xxxyyz_xy[k] = -g_z_0_y_0_xxyyz_xy[k] * ab_x + g_z_0_y_0_xxyyz_xxy[k];

                g_z_0_y_0_xxxyyz_xz[k] = -g_z_0_y_0_xxyyz_xz[k] * ab_x + g_z_0_y_0_xxyyz_xxz[k];

                g_z_0_y_0_xxxyyz_yy[k] = -g_z_0_y_0_xxyyz_yy[k] * ab_x + g_z_0_y_0_xxyyz_xyy[k];

                g_z_0_y_0_xxxyyz_yz[k] = -g_z_0_y_0_xxyyz_yz[k] * ab_x + g_z_0_y_0_xxyyz_xyz[k];

                g_z_0_y_0_xxxyyz_zz[k] = -g_z_0_y_0_xxyyz_zz[k] * ab_x + g_z_0_y_0_xxyyz_xzz[k];
            }

            /// Set up 1224-1230 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxyzz_xx = cbuffer.data(id_geom_1010_off + 1224 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyzz_xy = cbuffer.data(id_geom_1010_off + 1225 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyzz_xz = cbuffer.data(id_geom_1010_off + 1226 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyzz_yy = cbuffer.data(id_geom_1010_off + 1227 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyzz_yz = cbuffer.data(id_geom_1010_off + 1228 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyzz_zz = cbuffer.data(id_geom_1010_off + 1229 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxyzz_xx, g_z_0_y_0_xxxyzz_xy, g_z_0_y_0_xxxyzz_xz, g_z_0_y_0_xxxyzz_yy, g_z_0_y_0_xxxyzz_yz, g_z_0_y_0_xxxyzz_zz, g_z_0_y_0_xxyzz_xx, g_z_0_y_0_xxyzz_xxx, g_z_0_y_0_xxyzz_xxy, g_z_0_y_0_xxyzz_xxz, g_z_0_y_0_xxyzz_xy, g_z_0_y_0_xxyzz_xyy, g_z_0_y_0_xxyzz_xyz, g_z_0_y_0_xxyzz_xz, g_z_0_y_0_xxyzz_xzz, g_z_0_y_0_xxyzz_yy, g_z_0_y_0_xxyzz_yz, g_z_0_y_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxyzz_xx[k] = -g_z_0_y_0_xxyzz_xx[k] * ab_x + g_z_0_y_0_xxyzz_xxx[k];

                g_z_0_y_0_xxxyzz_xy[k] = -g_z_0_y_0_xxyzz_xy[k] * ab_x + g_z_0_y_0_xxyzz_xxy[k];

                g_z_0_y_0_xxxyzz_xz[k] = -g_z_0_y_0_xxyzz_xz[k] * ab_x + g_z_0_y_0_xxyzz_xxz[k];

                g_z_0_y_0_xxxyzz_yy[k] = -g_z_0_y_0_xxyzz_yy[k] * ab_x + g_z_0_y_0_xxyzz_xyy[k];

                g_z_0_y_0_xxxyzz_yz[k] = -g_z_0_y_0_xxyzz_yz[k] * ab_x + g_z_0_y_0_xxyzz_xyz[k];

                g_z_0_y_0_xxxyzz_zz[k] = -g_z_0_y_0_xxyzz_zz[k] * ab_x + g_z_0_y_0_xxyzz_xzz[k];
            }

            /// Set up 1230-1236 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxzzz_xx = cbuffer.data(id_geom_1010_off + 1230 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzzz_xy = cbuffer.data(id_geom_1010_off + 1231 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzzz_xz = cbuffer.data(id_geom_1010_off + 1232 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzzz_yy = cbuffer.data(id_geom_1010_off + 1233 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzzz_yz = cbuffer.data(id_geom_1010_off + 1234 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzzz_zz = cbuffer.data(id_geom_1010_off + 1235 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxzzz_xx, g_z_0_y_0_xxxzzz_xy, g_z_0_y_0_xxxzzz_xz, g_z_0_y_0_xxxzzz_yy, g_z_0_y_0_xxxzzz_yz, g_z_0_y_0_xxxzzz_zz, g_z_0_y_0_xxzzz_xx, g_z_0_y_0_xxzzz_xxx, g_z_0_y_0_xxzzz_xxy, g_z_0_y_0_xxzzz_xxz, g_z_0_y_0_xxzzz_xy, g_z_0_y_0_xxzzz_xyy, g_z_0_y_0_xxzzz_xyz, g_z_0_y_0_xxzzz_xz, g_z_0_y_0_xxzzz_xzz, g_z_0_y_0_xxzzz_yy, g_z_0_y_0_xxzzz_yz, g_z_0_y_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxzzz_xx[k] = -g_z_0_y_0_xxzzz_xx[k] * ab_x + g_z_0_y_0_xxzzz_xxx[k];

                g_z_0_y_0_xxxzzz_xy[k] = -g_z_0_y_0_xxzzz_xy[k] * ab_x + g_z_0_y_0_xxzzz_xxy[k];

                g_z_0_y_0_xxxzzz_xz[k] = -g_z_0_y_0_xxzzz_xz[k] * ab_x + g_z_0_y_0_xxzzz_xxz[k];

                g_z_0_y_0_xxxzzz_yy[k] = -g_z_0_y_0_xxzzz_yy[k] * ab_x + g_z_0_y_0_xxzzz_xyy[k];

                g_z_0_y_0_xxxzzz_yz[k] = -g_z_0_y_0_xxzzz_yz[k] * ab_x + g_z_0_y_0_xxzzz_xyz[k];

                g_z_0_y_0_xxxzzz_zz[k] = -g_z_0_y_0_xxzzz_zz[k] * ab_x + g_z_0_y_0_xxzzz_xzz[k];
            }

            /// Set up 1236-1242 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyyyy_xx = cbuffer.data(id_geom_1010_off + 1236 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyy_xy = cbuffer.data(id_geom_1010_off + 1237 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyy_xz = cbuffer.data(id_geom_1010_off + 1238 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyy_yy = cbuffer.data(id_geom_1010_off + 1239 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyy_yz = cbuffer.data(id_geom_1010_off + 1240 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyy_zz = cbuffer.data(id_geom_1010_off + 1241 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyyyy_xx, g_z_0_y_0_xxyyyy_xy, g_z_0_y_0_xxyyyy_xz, g_z_0_y_0_xxyyyy_yy, g_z_0_y_0_xxyyyy_yz, g_z_0_y_0_xxyyyy_zz, g_z_0_y_0_xyyyy_xx, g_z_0_y_0_xyyyy_xxx, g_z_0_y_0_xyyyy_xxy, g_z_0_y_0_xyyyy_xxz, g_z_0_y_0_xyyyy_xy, g_z_0_y_0_xyyyy_xyy, g_z_0_y_0_xyyyy_xyz, g_z_0_y_0_xyyyy_xz, g_z_0_y_0_xyyyy_xzz, g_z_0_y_0_xyyyy_yy, g_z_0_y_0_xyyyy_yz, g_z_0_y_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyyyy_xx[k] = -g_z_0_y_0_xyyyy_xx[k] * ab_x + g_z_0_y_0_xyyyy_xxx[k];

                g_z_0_y_0_xxyyyy_xy[k] = -g_z_0_y_0_xyyyy_xy[k] * ab_x + g_z_0_y_0_xyyyy_xxy[k];

                g_z_0_y_0_xxyyyy_xz[k] = -g_z_0_y_0_xyyyy_xz[k] * ab_x + g_z_0_y_0_xyyyy_xxz[k];

                g_z_0_y_0_xxyyyy_yy[k] = -g_z_0_y_0_xyyyy_yy[k] * ab_x + g_z_0_y_0_xyyyy_xyy[k];

                g_z_0_y_0_xxyyyy_yz[k] = -g_z_0_y_0_xyyyy_yz[k] * ab_x + g_z_0_y_0_xyyyy_xyz[k];

                g_z_0_y_0_xxyyyy_zz[k] = -g_z_0_y_0_xyyyy_zz[k] * ab_x + g_z_0_y_0_xyyyy_xzz[k];
            }

            /// Set up 1242-1248 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyyyz_xx = cbuffer.data(id_geom_1010_off + 1242 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyz_xy = cbuffer.data(id_geom_1010_off + 1243 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyz_xz = cbuffer.data(id_geom_1010_off + 1244 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyz_yy = cbuffer.data(id_geom_1010_off + 1245 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyz_yz = cbuffer.data(id_geom_1010_off + 1246 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyz_zz = cbuffer.data(id_geom_1010_off + 1247 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyyyz_xx, g_z_0_y_0_xxyyyz_xy, g_z_0_y_0_xxyyyz_xz, g_z_0_y_0_xxyyyz_yy, g_z_0_y_0_xxyyyz_yz, g_z_0_y_0_xxyyyz_zz, g_z_0_y_0_xyyyz_xx, g_z_0_y_0_xyyyz_xxx, g_z_0_y_0_xyyyz_xxy, g_z_0_y_0_xyyyz_xxz, g_z_0_y_0_xyyyz_xy, g_z_0_y_0_xyyyz_xyy, g_z_0_y_0_xyyyz_xyz, g_z_0_y_0_xyyyz_xz, g_z_0_y_0_xyyyz_xzz, g_z_0_y_0_xyyyz_yy, g_z_0_y_0_xyyyz_yz, g_z_0_y_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyyyz_xx[k] = -g_z_0_y_0_xyyyz_xx[k] * ab_x + g_z_0_y_0_xyyyz_xxx[k];

                g_z_0_y_0_xxyyyz_xy[k] = -g_z_0_y_0_xyyyz_xy[k] * ab_x + g_z_0_y_0_xyyyz_xxy[k];

                g_z_0_y_0_xxyyyz_xz[k] = -g_z_0_y_0_xyyyz_xz[k] * ab_x + g_z_0_y_0_xyyyz_xxz[k];

                g_z_0_y_0_xxyyyz_yy[k] = -g_z_0_y_0_xyyyz_yy[k] * ab_x + g_z_0_y_0_xyyyz_xyy[k];

                g_z_0_y_0_xxyyyz_yz[k] = -g_z_0_y_0_xyyyz_yz[k] * ab_x + g_z_0_y_0_xyyyz_xyz[k];

                g_z_0_y_0_xxyyyz_zz[k] = -g_z_0_y_0_xyyyz_zz[k] * ab_x + g_z_0_y_0_xyyyz_xzz[k];
            }

            /// Set up 1248-1254 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyyzz_xx = cbuffer.data(id_geom_1010_off + 1248 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyzz_xy = cbuffer.data(id_geom_1010_off + 1249 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyzz_xz = cbuffer.data(id_geom_1010_off + 1250 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyzz_yy = cbuffer.data(id_geom_1010_off + 1251 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyzz_yz = cbuffer.data(id_geom_1010_off + 1252 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyzz_zz = cbuffer.data(id_geom_1010_off + 1253 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyyzz_xx, g_z_0_y_0_xxyyzz_xy, g_z_0_y_0_xxyyzz_xz, g_z_0_y_0_xxyyzz_yy, g_z_0_y_0_xxyyzz_yz, g_z_0_y_0_xxyyzz_zz, g_z_0_y_0_xyyzz_xx, g_z_0_y_0_xyyzz_xxx, g_z_0_y_0_xyyzz_xxy, g_z_0_y_0_xyyzz_xxz, g_z_0_y_0_xyyzz_xy, g_z_0_y_0_xyyzz_xyy, g_z_0_y_0_xyyzz_xyz, g_z_0_y_0_xyyzz_xz, g_z_0_y_0_xyyzz_xzz, g_z_0_y_0_xyyzz_yy, g_z_0_y_0_xyyzz_yz, g_z_0_y_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyyzz_xx[k] = -g_z_0_y_0_xyyzz_xx[k] * ab_x + g_z_0_y_0_xyyzz_xxx[k];

                g_z_0_y_0_xxyyzz_xy[k] = -g_z_0_y_0_xyyzz_xy[k] * ab_x + g_z_0_y_0_xyyzz_xxy[k];

                g_z_0_y_0_xxyyzz_xz[k] = -g_z_0_y_0_xyyzz_xz[k] * ab_x + g_z_0_y_0_xyyzz_xxz[k];

                g_z_0_y_0_xxyyzz_yy[k] = -g_z_0_y_0_xyyzz_yy[k] * ab_x + g_z_0_y_0_xyyzz_xyy[k];

                g_z_0_y_0_xxyyzz_yz[k] = -g_z_0_y_0_xyyzz_yz[k] * ab_x + g_z_0_y_0_xyyzz_xyz[k];

                g_z_0_y_0_xxyyzz_zz[k] = -g_z_0_y_0_xyyzz_zz[k] * ab_x + g_z_0_y_0_xyyzz_xzz[k];
            }

            /// Set up 1254-1260 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyzzz_xx = cbuffer.data(id_geom_1010_off + 1254 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzzz_xy = cbuffer.data(id_geom_1010_off + 1255 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzzz_xz = cbuffer.data(id_geom_1010_off + 1256 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzzz_yy = cbuffer.data(id_geom_1010_off + 1257 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzzz_yz = cbuffer.data(id_geom_1010_off + 1258 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzzz_zz = cbuffer.data(id_geom_1010_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyzzz_xx, g_z_0_y_0_xxyzzz_xy, g_z_0_y_0_xxyzzz_xz, g_z_0_y_0_xxyzzz_yy, g_z_0_y_0_xxyzzz_yz, g_z_0_y_0_xxyzzz_zz, g_z_0_y_0_xyzzz_xx, g_z_0_y_0_xyzzz_xxx, g_z_0_y_0_xyzzz_xxy, g_z_0_y_0_xyzzz_xxz, g_z_0_y_0_xyzzz_xy, g_z_0_y_0_xyzzz_xyy, g_z_0_y_0_xyzzz_xyz, g_z_0_y_0_xyzzz_xz, g_z_0_y_0_xyzzz_xzz, g_z_0_y_0_xyzzz_yy, g_z_0_y_0_xyzzz_yz, g_z_0_y_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyzzz_xx[k] = -g_z_0_y_0_xyzzz_xx[k] * ab_x + g_z_0_y_0_xyzzz_xxx[k];

                g_z_0_y_0_xxyzzz_xy[k] = -g_z_0_y_0_xyzzz_xy[k] * ab_x + g_z_0_y_0_xyzzz_xxy[k];

                g_z_0_y_0_xxyzzz_xz[k] = -g_z_0_y_0_xyzzz_xz[k] * ab_x + g_z_0_y_0_xyzzz_xxz[k];

                g_z_0_y_0_xxyzzz_yy[k] = -g_z_0_y_0_xyzzz_yy[k] * ab_x + g_z_0_y_0_xyzzz_xyy[k];

                g_z_0_y_0_xxyzzz_yz[k] = -g_z_0_y_0_xyzzz_yz[k] * ab_x + g_z_0_y_0_xyzzz_xyz[k];

                g_z_0_y_0_xxyzzz_zz[k] = -g_z_0_y_0_xyzzz_zz[k] * ab_x + g_z_0_y_0_xyzzz_xzz[k];
            }

            /// Set up 1260-1266 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxzzzz_xx = cbuffer.data(id_geom_1010_off + 1260 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzzz_xy = cbuffer.data(id_geom_1010_off + 1261 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzzz_xz = cbuffer.data(id_geom_1010_off + 1262 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzzz_yy = cbuffer.data(id_geom_1010_off + 1263 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzzz_yz = cbuffer.data(id_geom_1010_off + 1264 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzzz_zz = cbuffer.data(id_geom_1010_off + 1265 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxzzzz_xx, g_z_0_y_0_xxzzzz_xy, g_z_0_y_0_xxzzzz_xz, g_z_0_y_0_xxzzzz_yy, g_z_0_y_0_xxzzzz_yz, g_z_0_y_0_xxzzzz_zz, g_z_0_y_0_xzzzz_xx, g_z_0_y_0_xzzzz_xxx, g_z_0_y_0_xzzzz_xxy, g_z_0_y_0_xzzzz_xxz, g_z_0_y_0_xzzzz_xy, g_z_0_y_0_xzzzz_xyy, g_z_0_y_0_xzzzz_xyz, g_z_0_y_0_xzzzz_xz, g_z_0_y_0_xzzzz_xzz, g_z_0_y_0_xzzzz_yy, g_z_0_y_0_xzzzz_yz, g_z_0_y_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxzzzz_xx[k] = -g_z_0_y_0_xzzzz_xx[k] * ab_x + g_z_0_y_0_xzzzz_xxx[k];

                g_z_0_y_0_xxzzzz_xy[k] = -g_z_0_y_0_xzzzz_xy[k] * ab_x + g_z_0_y_0_xzzzz_xxy[k];

                g_z_0_y_0_xxzzzz_xz[k] = -g_z_0_y_0_xzzzz_xz[k] * ab_x + g_z_0_y_0_xzzzz_xxz[k];

                g_z_0_y_0_xxzzzz_yy[k] = -g_z_0_y_0_xzzzz_yy[k] * ab_x + g_z_0_y_0_xzzzz_xyy[k];

                g_z_0_y_0_xxzzzz_yz[k] = -g_z_0_y_0_xzzzz_yz[k] * ab_x + g_z_0_y_0_xzzzz_xyz[k];

                g_z_0_y_0_xxzzzz_zz[k] = -g_z_0_y_0_xzzzz_zz[k] * ab_x + g_z_0_y_0_xzzzz_xzz[k];
            }

            /// Set up 1266-1272 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyyyy_xx = cbuffer.data(id_geom_1010_off + 1266 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyy_xy = cbuffer.data(id_geom_1010_off + 1267 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyy_xz = cbuffer.data(id_geom_1010_off + 1268 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyy_yy = cbuffer.data(id_geom_1010_off + 1269 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyy_yz = cbuffer.data(id_geom_1010_off + 1270 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyy_zz = cbuffer.data(id_geom_1010_off + 1271 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyyyy_xx, g_z_0_y_0_xyyyyy_xy, g_z_0_y_0_xyyyyy_xz, g_z_0_y_0_xyyyyy_yy, g_z_0_y_0_xyyyyy_yz, g_z_0_y_0_xyyyyy_zz, g_z_0_y_0_yyyyy_xx, g_z_0_y_0_yyyyy_xxx, g_z_0_y_0_yyyyy_xxy, g_z_0_y_0_yyyyy_xxz, g_z_0_y_0_yyyyy_xy, g_z_0_y_0_yyyyy_xyy, g_z_0_y_0_yyyyy_xyz, g_z_0_y_0_yyyyy_xz, g_z_0_y_0_yyyyy_xzz, g_z_0_y_0_yyyyy_yy, g_z_0_y_0_yyyyy_yz, g_z_0_y_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyyyy_xx[k] = -g_z_0_y_0_yyyyy_xx[k] * ab_x + g_z_0_y_0_yyyyy_xxx[k];

                g_z_0_y_0_xyyyyy_xy[k] = -g_z_0_y_0_yyyyy_xy[k] * ab_x + g_z_0_y_0_yyyyy_xxy[k];

                g_z_0_y_0_xyyyyy_xz[k] = -g_z_0_y_0_yyyyy_xz[k] * ab_x + g_z_0_y_0_yyyyy_xxz[k];

                g_z_0_y_0_xyyyyy_yy[k] = -g_z_0_y_0_yyyyy_yy[k] * ab_x + g_z_0_y_0_yyyyy_xyy[k];

                g_z_0_y_0_xyyyyy_yz[k] = -g_z_0_y_0_yyyyy_yz[k] * ab_x + g_z_0_y_0_yyyyy_xyz[k];

                g_z_0_y_0_xyyyyy_zz[k] = -g_z_0_y_0_yyyyy_zz[k] * ab_x + g_z_0_y_0_yyyyy_xzz[k];
            }

            /// Set up 1272-1278 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyyyz_xx = cbuffer.data(id_geom_1010_off + 1272 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyz_xy = cbuffer.data(id_geom_1010_off + 1273 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyz_xz = cbuffer.data(id_geom_1010_off + 1274 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyz_yy = cbuffer.data(id_geom_1010_off + 1275 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyz_yz = cbuffer.data(id_geom_1010_off + 1276 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyz_zz = cbuffer.data(id_geom_1010_off + 1277 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyyyz_xx, g_z_0_y_0_xyyyyz_xy, g_z_0_y_0_xyyyyz_xz, g_z_0_y_0_xyyyyz_yy, g_z_0_y_0_xyyyyz_yz, g_z_0_y_0_xyyyyz_zz, g_z_0_y_0_yyyyz_xx, g_z_0_y_0_yyyyz_xxx, g_z_0_y_0_yyyyz_xxy, g_z_0_y_0_yyyyz_xxz, g_z_0_y_0_yyyyz_xy, g_z_0_y_0_yyyyz_xyy, g_z_0_y_0_yyyyz_xyz, g_z_0_y_0_yyyyz_xz, g_z_0_y_0_yyyyz_xzz, g_z_0_y_0_yyyyz_yy, g_z_0_y_0_yyyyz_yz, g_z_0_y_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyyyz_xx[k] = -g_z_0_y_0_yyyyz_xx[k] * ab_x + g_z_0_y_0_yyyyz_xxx[k];

                g_z_0_y_0_xyyyyz_xy[k] = -g_z_0_y_0_yyyyz_xy[k] * ab_x + g_z_0_y_0_yyyyz_xxy[k];

                g_z_0_y_0_xyyyyz_xz[k] = -g_z_0_y_0_yyyyz_xz[k] * ab_x + g_z_0_y_0_yyyyz_xxz[k];

                g_z_0_y_0_xyyyyz_yy[k] = -g_z_0_y_0_yyyyz_yy[k] * ab_x + g_z_0_y_0_yyyyz_xyy[k];

                g_z_0_y_0_xyyyyz_yz[k] = -g_z_0_y_0_yyyyz_yz[k] * ab_x + g_z_0_y_0_yyyyz_xyz[k];

                g_z_0_y_0_xyyyyz_zz[k] = -g_z_0_y_0_yyyyz_zz[k] * ab_x + g_z_0_y_0_yyyyz_xzz[k];
            }

            /// Set up 1278-1284 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyyzz_xx = cbuffer.data(id_geom_1010_off + 1278 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyzz_xy = cbuffer.data(id_geom_1010_off + 1279 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyzz_xz = cbuffer.data(id_geom_1010_off + 1280 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyzz_yy = cbuffer.data(id_geom_1010_off + 1281 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyzz_yz = cbuffer.data(id_geom_1010_off + 1282 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyzz_zz = cbuffer.data(id_geom_1010_off + 1283 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyyzz_xx, g_z_0_y_0_xyyyzz_xy, g_z_0_y_0_xyyyzz_xz, g_z_0_y_0_xyyyzz_yy, g_z_0_y_0_xyyyzz_yz, g_z_0_y_0_xyyyzz_zz, g_z_0_y_0_yyyzz_xx, g_z_0_y_0_yyyzz_xxx, g_z_0_y_0_yyyzz_xxy, g_z_0_y_0_yyyzz_xxz, g_z_0_y_0_yyyzz_xy, g_z_0_y_0_yyyzz_xyy, g_z_0_y_0_yyyzz_xyz, g_z_0_y_0_yyyzz_xz, g_z_0_y_0_yyyzz_xzz, g_z_0_y_0_yyyzz_yy, g_z_0_y_0_yyyzz_yz, g_z_0_y_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyyzz_xx[k] = -g_z_0_y_0_yyyzz_xx[k] * ab_x + g_z_0_y_0_yyyzz_xxx[k];

                g_z_0_y_0_xyyyzz_xy[k] = -g_z_0_y_0_yyyzz_xy[k] * ab_x + g_z_0_y_0_yyyzz_xxy[k];

                g_z_0_y_0_xyyyzz_xz[k] = -g_z_0_y_0_yyyzz_xz[k] * ab_x + g_z_0_y_0_yyyzz_xxz[k];

                g_z_0_y_0_xyyyzz_yy[k] = -g_z_0_y_0_yyyzz_yy[k] * ab_x + g_z_0_y_0_yyyzz_xyy[k];

                g_z_0_y_0_xyyyzz_yz[k] = -g_z_0_y_0_yyyzz_yz[k] * ab_x + g_z_0_y_0_yyyzz_xyz[k];

                g_z_0_y_0_xyyyzz_zz[k] = -g_z_0_y_0_yyyzz_zz[k] * ab_x + g_z_0_y_0_yyyzz_xzz[k];
            }

            /// Set up 1284-1290 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyzzz_xx = cbuffer.data(id_geom_1010_off + 1284 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzzz_xy = cbuffer.data(id_geom_1010_off + 1285 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzzz_xz = cbuffer.data(id_geom_1010_off + 1286 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzzz_yy = cbuffer.data(id_geom_1010_off + 1287 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzzz_yz = cbuffer.data(id_geom_1010_off + 1288 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzzz_zz = cbuffer.data(id_geom_1010_off + 1289 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyzzz_xx, g_z_0_y_0_xyyzzz_xy, g_z_0_y_0_xyyzzz_xz, g_z_0_y_0_xyyzzz_yy, g_z_0_y_0_xyyzzz_yz, g_z_0_y_0_xyyzzz_zz, g_z_0_y_0_yyzzz_xx, g_z_0_y_0_yyzzz_xxx, g_z_0_y_0_yyzzz_xxy, g_z_0_y_0_yyzzz_xxz, g_z_0_y_0_yyzzz_xy, g_z_0_y_0_yyzzz_xyy, g_z_0_y_0_yyzzz_xyz, g_z_0_y_0_yyzzz_xz, g_z_0_y_0_yyzzz_xzz, g_z_0_y_0_yyzzz_yy, g_z_0_y_0_yyzzz_yz, g_z_0_y_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyzzz_xx[k] = -g_z_0_y_0_yyzzz_xx[k] * ab_x + g_z_0_y_0_yyzzz_xxx[k];

                g_z_0_y_0_xyyzzz_xy[k] = -g_z_0_y_0_yyzzz_xy[k] * ab_x + g_z_0_y_0_yyzzz_xxy[k];

                g_z_0_y_0_xyyzzz_xz[k] = -g_z_0_y_0_yyzzz_xz[k] * ab_x + g_z_0_y_0_yyzzz_xxz[k];

                g_z_0_y_0_xyyzzz_yy[k] = -g_z_0_y_0_yyzzz_yy[k] * ab_x + g_z_0_y_0_yyzzz_xyy[k];

                g_z_0_y_0_xyyzzz_yz[k] = -g_z_0_y_0_yyzzz_yz[k] * ab_x + g_z_0_y_0_yyzzz_xyz[k];

                g_z_0_y_0_xyyzzz_zz[k] = -g_z_0_y_0_yyzzz_zz[k] * ab_x + g_z_0_y_0_yyzzz_xzz[k];
            }

            /// Set up 1290-1296 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyzzzz_xx = cbuffer.data(id_geom_1010_off + 1290 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzzz_xy = cbuffer.data(id_geom_1010_off + 1291 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzzz_xz = cbuffer.data(id_geom_1010_off + 1292 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzzz_yy = cbuffer.data(id_geom_1010_off + 1293 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzzz_yz = cbuffer.data(id_geom_1010_off + 1294 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzzz_zz = cbuffer.data(id_geom_1010_off + 1295 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyzzzz_xx, g_z_0_y_0_xyzzzz_xy, g_z_0_y_0_xyzzzz_xz, g_z_0_y_0_xyzzzz_yy, g_z_0_y_0_xyzzzz_yz, g_z_0_y_0_xyzzzz_zz, g_z_0_y_0_yzzzz_xx, g_z_0_y_0_yzzzz_xxx, g_z_0_y_0_yzzzz_xxy, g_z_0_y_0_yzzzz_xxz, g_z_0_y_0_yzzzz_xy, g_z_0_y_0_yzzzz_xyy, g_z_0_y_0_yzzzz_xyz, g_z_0_y_0_yzzzz_xz, g_z_0_y_0_yzzzz_xzz, g_z_0_y_0_yzzzz_yy, g_z_0_y_0_yzzzz_yz, g_z_0_y_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyzzzz_xx[k] = -g_z_0_y_0_yzzzz_xx[k] * ab_x + g_z_0_y_0_yzzzz_xxx[k];

                g_z_0_y_0_xyzzzz_xy[k] = -g_z_0_y_0_yzzzz_xy[k] * ab_x + g_z_0_y_0_yzzzz_xxy[k];

                g_z_0_y_0_xyzzzz_xz[k] = -g_z_0_y_0_yzzzz_xz[k] * ab_x + g_z_0_y_0_yzzzz_xxz[k];

                g_z_0_y_0_xyzzzz_yy[k] = -g_z_0_y_0_yzzzz_yy[k] * ab_x + g_z_0_y_0_yzzzz_xyy[k];

                g_z_0_y_0_xyzzzz_yz[k] = -g_z_0_y_0_yzzzz_yz[k] * ab_x + g_z_0_y_0_yzzzz_xyz[k];

                g_z_0_y_0_xyzzzz_zz[k] = -g_z_0_y_0_yzzzz_zz[k] * ab_x + g_z_0_y_0_yzzzz_xzz[k];
            }

            /// Set up 1296-1302 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xzzzzz_xx = cbuffer.data(id_geom_1010_off + 1296 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzzz_xy = cbuffer.data(id_geom_1010_off + 1297 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzzz_xz = cbuffer.data(id_geom_1010_off + 1298 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzzz_yy = cbuffer.data(id_geom_1010_off + 1299 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzzz_yz = cbuffer.data(id_geom_1010_off + 1300 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzzz_zz = cbuffer.data(id_geom_1010_off + 1301 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xzzzzz_xx, g_z_0_y_0_xzzzzz_xy, g_z_0_y_0_xzzzzz_xz, g_z_0_y_0_xzzzzz_yy, g_z_0_y_0_xzzzzz_yz, g_z_0_y_0_xzzzzz_zz, g_z_0_y_0_zzzzz_xx, g_z_0_y_0_zzzzz_xxx, g_z_0_y_0_zzzzz_xxy, g_z_0_y_0_zzzzz_xxz, g_z_0_y_0_zzzzz_xy, g_z_0_y_0_zzzzz_xyy, g_z_0_y_0_zzzzz_xyz, g_z_0_y_0_zzzzz_xz, g_z_0_y_0_zzzzz_xzz, g_z_0_y_0_zzzzz_yy, g_z_0_y_0_zzzzz_yz, g_z_0_y_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xzzzzz_xx[k] = -g_z_0_y_0_zzzzz_xx[k] * ab_x + g_z_0_y_0_zzzzz_xxx[k];

                g_z_0_y_0_xzzzzz_xy[k] = -g_z_0_y_0_zzzzz_xy[k] * ab_x + g_z_0_y_0_zzzzz_xxy[k];

                g_z_0_y_0_xzzzzz_xz[k] = -g_z_0_y_0_zzzzz_xz[k] * ab_x + g_z_0_y_0_zzzzz_xxz[k];

                g_z_0_y_0_xzzzzz_yy[k] = -g_z_0_y_0_zzzzz_yy[k] * ab_x + g_z_0_y_0_zzzzz_xyy[k];

                g_z_0_y_0_xzzzzz_yz[k] = -g_z_0_y_0_zzzzz_yz[k] * ab_x + g_z_0_y_0_zzzzz_xyz[k];

                g_z_0_y_0_xzzzzz_zz[k] = -g_z_0_y_0_zzzzz_zz[k] * ab_x + g_z_0_y_0_zzzzz_xzz[k];
            }

            /// Set up 1302-1308 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyyyy_xx = cbuffer.data(id_geom_1010_off + 1302 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyy_xy = cbuffer.data(id_geom_1010_off + 1303 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyy_xz = cbuffer.data(id_geom_1010_off + 1304 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyy_yy = cbuffer.data(id_geom_1010_off + 1305 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyy_yz = cbuffer.data(id_geom_1010_off + 1306 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyy_zz = cbuffer.data(id_geom_1010_off + 1307 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyyy_xx, g_z_0_y_0_yyyyy_xxy, g_z_0_y_0_yyyyy_xy, g_z_0_y_0_yyyyy_xyy, g_z_0_y_0_yyyyy_xyz, g_z_0_y_0_yyyyy_xz, g_z_0_y_0_yyyyy_yy, g_z_0_y_0_yyyyy_yyy, g_z_0_y_0_yyyyy_yyz, g_z_0_y_0_yyyyy_yz, g_z_0_y_0_yyyyy_yzz, g_z_0_y_0_yyyyy_zz, g_z_0_y_0_yyyyyy_xx, g_z_0_y_0_yyyyyy_xy, g_z_0_y_0_yyyyyy_xz, g_z_0_y_0_yyyyyy_yy, g_z_0_y_0_yyyyyy_yz, g_z_0_y_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyyyy_xx[k] = -g_z_0_y_0_yyyyy_xx[k] * ab_y + g_z_0_y_0_yyyyy_xxy[k];

                g_z_0_y_0_yyyyyy_xy[k] = -g_z_0_y_0_yyyyy_xy[k] * ab_y + g_z_0_y_0_yyyyy_xyy[k];

                g_z_0_y_0_yyyyyy_xz[k] = -g_z_0_y_0_yyyyy_xz[k] * ab_y + g_z_0_y_0_yyyyy_xyz[k];

                g_z_0_y_0_yyyyyy_yy[k] = -g_z_0_y_0_yyyyy_yy[k] * ab_y + g_z_0_y_0_yyyyy_yyy[k];

                g_z_0_y_0_yyyyyy_yz[k] = -g_z_0_y_0_yyyyy_yz[k] * ab_y + g_z_0_y_0_yyyyy_yyz[k];

                g_z_0_y_0_yyyyyy_zz[k] = -g_z_0_y_0_yyyyy_zz[k] * ab_y + g_z_0_y_0_yyyyy_yzz[k];
            }

            /// Set up 1308-1314 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyyyz_xx = cbuffer.data(id_geom_1010_off + 1308 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyz_xy = cbuffer.data(id_geom_1010_off + 1309 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyz_xz = cbuffer.data(id_geom_1010_off + 1310 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyz_yy = cbuffer.data(id_geom_1010_off + 1311 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyz_yz = cbuffer.data(id_geom_1010_off + 1312 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyz_zz = cbuffer.data(id_geom_1010_off + 1313 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyyyz_xx, g_z_0_y_0_yyyyyz_xy, g_z_0_y_0_yyyyyz_xz, g_z_0_y_0_yyyyyz_yy, g_z_0_y_0_yyyyyz_yz, g_z_0_y_0_yyyyyz_zz, g_z_0_y_0_yyyyz_xx, g_z_0_y_0_yyyyz_xxy, g_z_0_y_0_yyyyz_xy, g_z_0_y_0_yyyyz_xyy, g_z_0_y_0_yyyyz_xyz, g_z_0_y_0_yyyyz_xz, g_z_0_y_0_yyyyz_yy, g_z_0_y_0_yyyyz_yyy, g_z_0_y_0_yyyyz_yyz, g_z_0_y_0_yyyyz_yz, g_z_0_y_0_yyyyz_yzz, g_z_0_y_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyyyz_xx[k] = -g_z_0_y_0_yyyyz_xx[k] * ab_y + g_z_0_y_0_yyyyz_xxy[k];

                g_z_0_y_0_yyyyyz_xy[k] = -g_z_0_y_0_yyyyz_xy[k] * ab_y + g_z_0_y_0_yyyyz_xyy[k];

                g_z_0_y_0_yyyyyz_xz[k] = -g_z_0_y_0_yyyyz_xz[k] * ab_y + g_z_0_y_0_yyyyz_xyz[k];

                g_z_0_y_0_yyyyyz_yy[k] = -g_z_0_y_0_yyyyz_yy[k] * ab_y + g_z_0_y_0_yyyyz_yyy[k];

                g_z_0_y_0_yyyyyz_yz[k] = -g_z_0_y_0_yyyyz_yz[k] * ab_y + g_z_0_y_0_yyyyz_yyz[k];

                g_z_0_y_0_yyyyyz_zz[k] = -g_z_0_y_0_yyyyz_zz[k] * ab_y + g_z_0_y_0_yyyyz_yzz[k];
            }

            /// Set up 1314-1320 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyyzz_xx = cbuffer.data(id_geom_1010_off + 1314 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyzz_xy = cbuffer.data(id_geom_1010_off + 1315 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyzz_xz = cbuffer.data(id_geom_1010_off + 1316 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyzz_yy = cbuffer.data(id_geom_1010_off + 1317 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyzz_yz = cbuffer.data(id_geom_1010_off + 1318 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyzz_zz = cbuffer.data(id_geom_1010_off + 1319 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyyzz_xx, g_z_0_y_0_yyyyzz_xy, g_z_0_y_0_yyyyzz_xz, g_z_0_y_0_yyyyzz_yy, g_z_0_y_0_yyyyzz_yz, g_z_0_y_0_yyyyzz_zz, g_z_0_y_0_yyyzz_xx, g_z_0_y_0_yyyzz_xxy, g_z_0_y_0_yyyzz_xy, g_z_0_y_0_yyyzz_xyy, g_z_0_y_0_yyyzz_xyz, g_z_0_y_0_yyyzz_xz, g_z_0_y_0_yyyzz_yy, g_z_0_y_0_yyyzz_yyy, g_z_0_y_0_yyyzz_yyz, g_z_0_y_0_yyyzz_yz, g_z_0_y_0_yyyzz_yzz, g_z_0_y_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyyzz_xx[k] = -g_z_0_y_0_yyyzz_xx[k] * ab_y + g_z_0_y_0_yyyzz_xxy[k];

                g_z_0_y_0_yyyyzz_xy[k] = -g_z_0_y_0_yyyzz_xy[k] * ab_y + g_z_0_y_0_yyyzz_xyy[k];

                g_z_0_y_0_yyyyzz_xz[k] = -g_z_0_y_0_yyyzz_xz[k] * ab_y + g_z_0_y_0_yyyzz_xyz[k];

                g_z_0_y_0_yyyyzz_yy[k] = -g_z_0_y_0_yyyzz_yy[k] * ab_y + g_z_0_y_0_yyyzz_yyy[k];

                g_z_0_y_0_yyyyzz_yz[k] = -g_z_0_y_0_yyyzz_yz[k] * ab_y + g_z_0_y_0_yyyzz_yyz[k];

                g_z_0_y_0_yyyyzz_zz[k] = -g_z_0_y_0_yyyzz_zz[k] * ab_y + g_z_0_y_0_yyyzz_yzz[k];
            }

            /// Set up 1320-1326 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyzzz_xx = cbuffer.data(id_geom_1010_off + 1320 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzzz_xy = cbuffer.data(id_geom_1010_off + 1321 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzzz_xz = cbuffer.data(id_geom_1010_off + 1322 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzzz_yy = cbuffer.data(id_geom_1010_off + 1323 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzzz_yz = cbuffer.data(id_geom_1010_off + 1324 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzzz_zz = cbuffer.data(id_geom_1010_off + 1325 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyzzz_xx, g_z_0_y_0_yyyzzz_xy, g_z_0_y_0_yyyzzz_xz, g_z_0_y_0_yyyzzz_yy, g_z_0_y_0_yyyzzz_yz, g_z_0_y_0_yyyzzz_zz, g_z_0_y_0_yyzzz_xx, g_z_0_y_0_yyzzz_xxy, g_z_0_y_0_yyzzz_xy, g_z_0_y_0_yyzzz_xyy, g_z_0_y_0_yyzzz_xyz, g_z_0_y_0_yyzzz_xz, g_z_0_y_0_yyzzz_yy, g_z_0_y_0_yyzzz_yyy, g_z_0_y_0_yyzzz_yyz, g_z_0_y_0_yyzzz_yz, g_z_0_y_0_yyzzz_yzz, g_z_0_y_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyzzz_xx[k] = -g_z_0_y_0_yyzzz_xx[k] * ab_y + g_z_0_y_0_yyzzz_xxy[k];

                g_z_0_y_0_yyyzzz_xy[k] = -g_z_0_y_0_yyzzz_xy[k] * ab_y + g_z_0_y_0_yyzzz_xyy[k];

                g_z_0_y_0_yyyzzz_xz[k] = -g_z_0_y_0_yyzzz_xz[k] * ab_y + g_z_0_y_0_yyzzz_xyz[k];

                g_z_0_y_0_yyyzzz_yy[k] = -g_z_0_y_0_yyzzz_yy[k] * ab_y + g_z_0_y_0_yyzzz_yyy[k];

                g_z_0_y_0_yyyzzz_yz[k] = -g_z_0_y_0_yyzzz_yz[k] * ab_y + g_z_0_y_0_yyzzz_yyz[k];

                g_z_0_y_0_yyyzzz_zz[k] = -g_z_0_y_0_yyzzz_zz[k] * ab_y + g_z_0_y_0_yyzzz_yzz[k];
            }

            /// Set up 1326-1332 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyzzzz_xx = cbuffer.data(id_geom_1010_off + 1326 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzzz_xy = cbuffer.data(id_geom_1010_off + 1327 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzzz_xz = cbuffer.data(id_geom_1010_off + 1328 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzzz_yy = cbuffer.data(id_geom_1010_off + 1329 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzzz_yz = cbuffer.data(id_geom_1010_off + 1330 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzzz_zz = cbuffer.data(id_geom_1010_off + 1331 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyzzzz_xx, g_z_0_y_0_yyzzzz_xy, g_z_0_y_0_yyzzzz_xz, g_z_0_y_0_yyzzzz_yy, g_z_0_y_0_yyzzzz_yz, g_z_0_y_0_yyzzzz_zz, g_z_0_y_0_yzzzz_xx, g_z_0_y_0_yzzzz_xxy, g_z_0_y_0_yzzzz_xy, g_z_0_y_0_yzzzz_xyy, g_z_0_y_0_yzzzz_xyz, g_z_0_y_0_yzzzz_xz, g_z_0_y_0_yzzzz_yy, g_z_0_y_0_yzzzz_yyy, g_z_0_y_0_yzzzz_yyz, g_z_0_y_0_yzzzz_yz, g_z_0_y_0_yzzzz_yzz, g_z_0_y_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyzzzz_xx[k] = -g_z_0_y_0_yzzzz_xx[k] * ab_y + g_z_0_y_0_yzzzz_xxy[k];

                g_z_0_y_0_yyzzzz_xy[k] = -g_z_0_y_0_yzzzz_xy[k] * ab_y + g_z_0_y_0_yzzzz_xyy[k];

                g_z_0_y_0_yyzzzz_xz[k] = -g_z_0_y_0_yzzzz_xz[k] * ab_y + g_z_0_y_0_yzzzz_xyz[k];

                g_z_0_y_0_yyzzzz_yy[k] = -g_z_0_y_0_yzzzz_yy[k] * ab_y + g_z_0_y_0_yzzzz_yyy[k];

                g_z_0_y_0_yyzzzz_yz[k] = -g_z_0_y_0_yzzzz_yz[k] * ab_y + g_z_0_y_0_yzzzz_yyz[k];

                g_z_0_y_0_yyzzzz_zz[k] = -g_z_0_y_0_yzzzz_zz[k] * ab_y + g_z_0_y_0_yzzzz_yzz[k];
            }

            /// Set up 1332-1338 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yzzzzz_xx = cbuffer.data(id_geom_1010_off + 1332 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzzz_xy = cbuffer.data(id_geom_1010_off + 1333 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzzz_xz = cbuffer.data(id_geom_1010_off + 1334 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzzz_yy = cbuffer.data(id_geom_1010_off + 1335 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzzz_yz = cbuffer.data(id_geom_1010_off + 1336 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzzz_zz = cbuffer.data(id_geom_1010_off + 1337 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yzzzzz_xx, g_z_0_y_0_yzzzzz_xy, g_z_0_y_0_yzzzzz_xz, g_z_0_y_0_yzzzzz_yy, g_z_0_y_0_yzzzzz_yz, g_z_0_y_0_yzzzzz_zz, g_z_0_y_0_zzzzz_xx, g_z_0_y_0_zzzzz_xxy, g_z_0_y_0_zzzzz_xy, g_z_0_y_0_zzzzz_xyy, g_z_0_y_0_zzzzz_xyz, g_z_0_y_0_zzzzz_xz, g_z_0_y_0_zzzzz_yy, g_z_0_y_0_zzzzz_yyy, g_z_0_y_0_zzzzz_yyz, g_z_0_y_0_zzzzz_yz, g_z_0_y_0_zzzzz_yzz, g_z_0_y_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yzzzzz_xx[k] = -g_z_0_y_0_zzzzz_xx[k] * ab_y + g_z_0_y_0_zzzzz_xxy[k];

                g_z_0_y_0_yzzzzz_xy[k] = -g_z_0_y_0_zzzzz_xy[k] * ab_y + g_z_0_y_0_zzzzz_xyy[k];

                g_z_0_y_0_yzzzzz_xz[k] = -g_z_0_y_0_zzzzz_xz[k] * ab_y + g_z_0_y_0_zzzzz_xyz[k];

                g_z_0_y_0_yzzzzz_yy[k] = -g_z_0_y_0_zzzzz_yy[k] * ab_y + g_z_0_y_0_zzzzz_yyy[k];

                g_z_0_y_0_yzzzzz_yz[k] = -g_z_0_y_0_zzzzz_yz[k] * ab_y + g_z_0_y_0_zzzzz_yyz[k];

                g_z_0_y_0_yzzzzz_zz[k] = -g_z_0_y_0_zzzzz_zz[k] * ab_y + g_z_0_y_0_zzzzz_yzz[k];
            }

            /// Set up 1338-1344 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_zzzzzz_xx = cbuffer.data(id_geom_1010_off + 1338 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzzz_xy = cbuffer.data(id_geom_1010_off + 1339 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzzz_xz = cbuffer.data(id_geom_1010_off + 1340 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzzz_yy = cbuffer.data(id_geom_1010_off + 1341 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzzz_yz = cbuffer.data(id_geom_1010_off + 1342 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzzz_zz = cbuffer.data(id_geom_1010_off + 1343 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_zzzzz_xx, g_0_0_y_0_zzzzz_xy, g_0_0_y_0_zzzzz_xz, g_0_0_y_0_zzzzz_yy, g_0_0_y_0_zzzzz_yz, g_0_0_y_0_zzzzz_zz, g_z_0_y_0_zzzzz_xx, g_z_0_y_0_zzzzz_xxz, g_z_0_y_0_zzzzz_xy, g_z_0_y_0_zzzzz_xyz, g_z_0_y_0_zzzzz_xz, g_z_0_y_0_zzzzz_xzz, g_z_0_y_0_zzzzz_yy, g_z_0_y_0_zzzzz_yyz, g_z_0_y_0_zzzzz_yz, g_z_0_y_0_zzzzz_yzz, g_z_0_y_0_zzzzz_zz, g_z_0_y_0_zzzzz_zzz, g_z_0_y_0_zzzzzz_xx, g_z_0_y_0_zzzzzz_xy, g_z_0_y_0_zzzzzz_xz, g_z_0_y_0_zzzzzz_yy, g_z_0_y_0_zzzzzz_yz, g_z_0_y_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_zzzzzz_xx[k] = -g_0_0_y_0_zzzzz_xx[k] - g_z_0_y_0_zzzzz_xx[k] * ab_z + g_z_0_y_0_zzzzz_xxz[k];

                g_z_0_y_0_zzzzzz_xy[k] = -g_0_0_y_0_zzzzz_xy[k] - g_z_0_y_0_zzzzz_xy[k] * ab_z + g_z_0_y_0_zzzzz_xyz[k];

                g_z_0_y_0_zzzzzz_xz[k] = -g_0_0_y_0_zzzzz_xz[k] - g_z_0_y_0_zzzzz_xz[k] * ab_z + g_z_0_y_0_zzzzz_xzz[k];

                g_z_0_y_0_zzzzzz_yy[k] = -g_0_0_y_0_zzzzz_yy[k] - g_z_0_y_0_zzzzz_yy[k] * ab_z + g_z_0_y_0_zzzzz_yyz[k];

                g_z_0_y_0_zzzzzz_yz[k] = -g_0_0_y_0_zzzzz_yz[k] - g_z_0_y_0_zzzzz_yz[k] * ab_z + g_z_0_y_0_zzzzz_yzz[k];

                g_z_0_y_0_zzzzzz_zz[k] = -g_0_0_y_0_zzzzz_zz[k] - g_z_0_y_0_zzzzz_zz[k] * ab_z + g_z_0_y_0_zzzzz_zzz[k];
            }

            /// Set up 1344-1350 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxxx_xx = cbuffer.data(id_geom_1010_off + 1344 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxx_xy = cbuffer.data(id_geom_1010_off + 1345 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxx_xz = cbuffer.data(id_geom_1010_off + 1346 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxx_yy = cbuffer.data(id_geom_1010_off + 1347 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxx_yz = cbuffer.data(id_geom_1010_off + 1348 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxx_zz = cbuffer.data(id_geom_1010_off + 1349 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxx_xx, g_z_0_z_0_xxxxx_xxx, g_z_0_z_0_xxxxx_xxy, g_z_0_z_0_xxxxx_xxz, g_z_0_z_0_xxxxx_xy, g_z_0_z_0_xxxxx_xyy, g_z_0_z_0_xxxxx_xyz, g_z_0_z_0_xxxxx_xz, g_z_0_z_0_xxxxx_xzz, g_z_0_z_0_xxxxx_yy, g_z_0_z_0_xxxxx_yz, g_z_0_z_0_xxxxx_zz, g_z_0_z_0_xxxxxx_xx, g_z_0_z_0_xxxxxx_xy, g_z_0_z_0_xxxxxx_xz, g_z_0_z_0_xxxxxx_yy, g_z_0_z_0_xxxxxx_yz, g_z_0_z_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxxx_xx[k] = -g_z_0_z_0_xxxxx_xx[k] * ab_x + g_z_0_z_0_xxxxx_xxx[k];

                g_z_0_z_0_xxxxxx_xy[k] = -g_z_0_z_0_xxxxx_xy[k] * ab_x + g_z_0_z_0_xxxxx_xxy[k];

                g_z_0_z_0_xxxxxx_xz[k] = -g_z_0_z_0_xxxxx_xz[k] * ab_x + g_z_0_z_0_xxxxx_xxz[k];

                g_z_0_z_0_xxxxxx_yy[k] = -g_z_0_z_0_xxxxx_yy[k] * ab_x + g_z_0_z_0_xxxxx_xyy[k];

                g_z_0_z_0_xxxxxx_yz[k] = -g_z_0_z_0_xxxxx_yz[k] * ab_x + g_z_0_z_0_xxxxx_xyz[k];

                g_z_0_z_0_xxxxxx_zz[k] = -g_z_0_z_0_xxxxx_zz[k] * ab_x + g_z_0_z_0_xxxxx_xzz[k];
            }

            /// Set up 1350-1356 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxxy_xx = cbuffer.data(id_geom_1010_off + 1350 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxy_xy = cbuffer.data(id_geom_1010_off + 1351 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxy_xz = cbuffer.data(id_geom_1010_off + 1352 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxy_yy = cbuffer.data(id_geom_1010_off + 1353 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxy_yz = cbuffer.data(id_geom_1010_off + 1354 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxy_zz = cbuffer.data(id_geom_1010_off + 1355 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxxy_xx, g_z_0_z_0_xxxxxy_xy, g_z_0_z_0_xxxxxy_xz, g_z_0_z_0_xxxxxy_yy, g_z_0_z_0_xxxxxy_yz, g_z_0_z_0_xxxxxy_zz, g_z_0_z_0_xxxxy_xx, g_z_0_z_0_xxxxy_xxx, g_z_0_z_0_xxxxy_xxy, g_z_0_z_0_xxxxy_xxz, g_z_0_z_0_xxxxy_xy, g_z_0_z_0_xxxxy_xyy, g_z_0_z_0_xxxxy_xyz, g_z_0_z_0_xxxxy_xz, g_z_0_z_0_xxxxy_xzz, g_z_0_z_0_xxxxy_yy, g_z_0_z_0_xxxxy_yz, g_z_0_z_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxxy_xx[k] = -g_z_0_z_0_xxxxy_xx[k] * ab_x + g_z_0_z_0_xxxxy_xxx[k];

                g_z_0_z_0_xxxxxy_xy[k] = -g_z_0_z_0_xxxxy_xy[k] * ab_x + g_z_0_z_0_xxxxy_xxy[k];

                g_z_0_z_0_xxxxxy_xz[k] = -g_z_0_z_0_xxxxy_xz[k] * ab_x + g_z_0_z_0_xxxxy_xxz[k];

                g_z_0_z_0_xxxxxy_yy[k] = -g_z_0_z_0_xxxxy_yy[k] * ab_x + g_z_0_z_0_xxxxy_xyy[k];

                g_z_0_z_0_xxxxxy_yz[k] = -g_z_0_z_0_xxxxy_yz[k] * ab_x + g_z_0_z_0_xxxxy_xyz[k];

                g_z_0_z_0_xxxxxy_zz[k] = -g_z_0_z_0_xxxxy_zz[k] * ab_x + g_z_0_z_0_xxxxy_xzz[k];
            }

            /// Set up 1356-1362 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxxz_xx = cbuffer.data(id_geom_1010_off + 1356 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxz_xy = cbuffer.data(id_geom_1010_off + 1357 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxz_xz = cbuffer.data(id_geom_1010_off + 1358 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxz_yy = cbuffer.data(id_geom_1010_off + 1359 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxz_yz = cbuffer.data(id_geom_1010_off + 1360 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxz_zz = cbuffer.data(id_geom_1010_off + 1361 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxxz_xx, g_z_0_z_0_xxxxxz_xy, g_z_0_z_0_xxxxxz_xz, g_z_0_z_0_xxxxxz_yy, g_z_0_z_0_xxxxxz_yz, g_z_0_z_0_xxxxxz_zz, g_z_0_z_0_xxxxz_xx, g_z_0_z_0_xxxxz_xxx, g_z_0_z_0_xxxxz_xxy, g_z_0_z_0_xxxxz_xxz, g_z_0_z_0_xxxxz_xy, g_z_0_z_0_xxxxz_xyy, g_z_0_z_0_xxxxz_xyz, g_z_0_z_0_xxxxz_xz, g_z_0_z_0_xxxxz_xzz, g_z_0_z_0_xxxxz_yy, g_z_0_z_0_xxxxz_yz, g_z_0_z_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxxz_xx[k] = -g_z_0_z_0_xxxxz_xx[k] * ab_x + g_z_0_z_0_xxxxz_xxx[k];

                g_z_0_z_0_xxxxxz_xy[k] = -g_z_0_z_0_xxxxz_xy[k] * ab_x + g_z_0_z_0_xxxxz_xxy[k];

                g_z_0_z_0_xxxxxz_xz[k] = -g_z_0_z_0_xxxxz_xz[k] * ab_x + g_z_0_z_0_xxxxz_xxz[k];

                g_z_0_z_0_xxxxxz_yy[k] = -g_z_0_z_0_xxxxz_yy[k] * ab_x + g_z_0_z_0_xxxxz_xyy[k];

                g_z_0_z_0_xxxxxz_yz[k] = -g_z_0_z_0_xxxxz_yz[k] * ab_x + g_z_0_z_0_xxxxz_xyz[k];

                g_z_0_z_0_xxxxxz_zz[k] = -g_z_0_z_0_xxxxz_zz[k] * ab_x + g_z_0_z_0_xxxxz_xzz[k];
            }

            /// Set up 1362-1368 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxyy_xx = cbuffer.data(id_geom_1010_off + 1362 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyy_xy = cbuffer.data(id_geom_1010_off + 1363 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyy_xz = cbuffer.data(id_geom_1010_off + 1364 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyy_yy = cbuffer.data(id_geom_1010_off + 1365 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyy_yz = cbuffer.data(id_geom_1010_off + 1366 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyy_zz = cbuffer.data(id_geom_1010_off + 1367 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxyy_xx, g_z_0_z_0_xxxxyy_xy, g_z_0_z_0_xxxxyy_xz, g_z_0_z_0_xxxxyy_yy, g_z_0_z_0_xxxxyy_yz, g_z_0_z_0_xxxxyy_zz, g_z_0_z_0_xxxyy_xx, g_z_0_z_0_xxxyy_xxx, g_z_0_z_0_xxxyy_xxy, g_z_0_z_0_xxxyy_xxz, g_z_0_z_0_xxxyy_xy, g_z_0_z_0_xxxyy_xyy, g_z_0_z_0_xxxyy_xyz, g_z_0_z_0_xxxyy_xz, g_z_0_z_0_xxxyy_xzz, g_z_0_z_0_xxxyy_yy, g_z_0_z_0_xxxyy_yz, g_z_0_z_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxyy_xx[k] = -g_z_0_z_0_xxxyy_xx[k] * ab_x + g_z_0_z_0_xxxyy_xxx[k];

                g_z_0_z_0_xxxxyy_xy[k] = -g_z_0_z_0_xxxyy_xy[k] * ab_x + g_z_0_z_0_xxxyy_xxy[k];

                g_z_0_z_0_xxxxyy_xz[k] = -g_z_0_z_0_xxxyy_xz[k] * ab_x + g_z_0_z_0_xxxyy_xxz[k];

                g_z_0_z_0_xxxxyy_yy[k] = -g_z_0_z_0_xxxyy_yy[k] * ab_x + g_z_0_z_0_xxxyy_xyy[k];

                g_z_0_z_0_xxxxyy_yz[k] = -g_z_0_z_0_xxxyy_yz[k] * ab_x + g_z_0_z_0_xxxyy_xyz[k];

                g_z_0_z_0_xxxxyy_zz[k] = -g_z_0_z_0_xxxyy_zz[k] * ab_x + g_z_0_z_0_xxxyy_xzz[k];
            }

            /// Set up 1368-1374 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxyz_xx = cbuffer.data(id_geom_1010_off + 1368 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyz_xy = cbuffer.data(id_geom_1010_off + 1369 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyz_xz = cbuffer.data(id_geom_1010_off + 1370 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyz_yy = cbuffer.data(id_geom_1010_off + 1371 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyz_yz = cbuffer.data(id_geom_1010_off + 1372 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyz_zz = cbuffer.data(id_geom_1010_off + 1373 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxyz_xx, g_z_0_z_0_xxxxyz_xy, g_z_0_z_0_xxxxyz_xz, g_z_0_z_0_xxxxyz_yy, g_z_0_z_0_xxxxyz_yz, g_z_0_z_0_xxxxyz_zz, g_z_0_z_0_xxxyz_xx, g_z_0_z_0_xxxyz_xxx, g_z_0_z_0_xxxyz_xxy, g_z_0_z_0_xxxyz_xxz, g_z_0_z_0_xxxyz_xy, g_z_0_z_0_xxxyz_xyy, g_z_0_z_0_xxxyz_xyz, g_z_0_z_0_xxxyz_xz, g_z_0_z_0_xxxyz_xzz, g_z_0_z_0_xxxyz_yy, g_z_0_z_0_xxxyz_yz, g_z_0_z_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxyz_xx[k] = -g_z_0_z_0_xxxyz_xx[k] * ab_x + g_z_0_z_0_xxxyz_xxx[k];

                g_z_0_z_0_xxxxyz_xy[k] = -g_z_0_z_0_xxxyz_xy[k] * ab_x + g_z_0_z_0_xxxyz_xxy[k];

                g_z_0_z_0_xxxxyz_xz[k] = -g_z_0_z_0_xxxyz_xz[k] * ab_x + g_z_0_z_0_xxxyz_xxz[k];

                g_z_0_z_0_xxxxyz_yy[k] = -g_z_0_z_0_xxxyz_yy[k] * ab_x + g_z_0_z_0_xxxyz_xyy[k];

                g_z_0_z_0_xxxxyz_yz[k] = -g_z_0_z_0_xxxyz_yz[k] * ab_x + g_z_0_z_0_xxxyz_xyz[k];

                g_z_0_z_0_xxxxyz_zz[k] = -g_z_0_z_0_xxxyz_zz[k] * ab_x + g_z_0_z_0_xxxyz_xzz[k];
            }

            /// Set up 1374-1380 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxzz_xx = cbuffer.data(id_geom_1010_off + 1374 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxzz_xy = cbuffer.data(id_geom_1010_off + 1375 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxzz_xz = cbuffer.data(id_geom_1010_off + 1376 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxzz_yy = cbuffer.data(id_geom_1010_off + 1377 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxzz_yz = cbuffer.data(id_geom_1010_off + 1378 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxzz_zz = cbuffer.data(id_geom_1010_off + 1379 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxzz_xx, g_z_0_z_0_xxxxzz_xy, g_z_0_z_0_xxxxzz_xz, g_z_0_z_0_xxxxzz_yy, g_z_0_z_0_xxxxzz_yz, g_z_0_z_0_xxxxzz_zz, g_z_0_z_0_xxxzz_xx, g_z_0_z_0_xxxzz_xxx, g_z_0_z_0_xxxzz_xxy, g_z_0_z_0_xxxzz_xxz, g_z_0_z_0_xxxzz_xy, g_z_0_z_0_xxxzz_xyy, g_z_0_z_0_xxxzz_xyz, g_z_0_z_0_xxxzz_xz, g_z_0_z_0_xxxzz_xzz, g_z_0_z_0_xxxzz_yy, g_z_0_z_0_xxxzz_yz, g_z_0_z_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxzz_xx[k] = -g_z_0_z_0_xxxzz_xx[k] * ab_x + g_z_0_z_0_xxxzz_xxx[k];

                g_z_0_z_0_xxxxzz_xy[k] = -g_z_0_z_0_xxxzz_xy[k] * ab_x + g_z_0_z_0_xxxzz_xxy[k];

                g_z_0_z_0_xxxxzz_xz[k] = -g_z_0_z_0_xxxzz_xz[k] * ab_x + g_z_0_z_0_xxxzz_xxz[k];

                g_z_0_z_0_xxxxzz_yy[k] = -g_z_0_z_0_xxxzz_yy[k] * ab_x + g_z_0_z_0_xxxzz_xyy[k];

                g_z_0_z_0_xxxxzz_yz[k] = -g_z_0_z_0_xxxzz_yz[k] * ab_x + g_z_0_z_0_xxxzz_xyz[k];

                g_z_0_z_0_xxxxzz_zz[k] = -g_z_0_z_0_xxxzz_zz[k] * ab_x + g_z_0_z_0_xxxzz_xzz[k];
            }

            /// Set up 1380-1386 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxyyy_xx = cbuffer.data(id_geom_1010_off + 1380 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyy_xy = cbuffer.data(id_geom_1010_off + 1381 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyy_xz = cbuffer.data(id_geom_1010_off + 1382 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyy_yy = cbuffer.data(id_geom_1010_off + 1383 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyy_yz = cbuffer.data(id_geom_1010_off + 1384 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyy_zz = cbuffer.data(id_geom_1010_off + 1385 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxyyy_xx, g_z_0_z_0_xxxyyy_xy, g_z_0_z_0_xxxyyy_xz, g_z_0_z_0_xxxyyy_yy, g_z_0_z_0_xxxyyy_yz, g_z_0_z_0_xxxyyy_zz, g_z_0_z_0_xxyyy_xx, g_z_0_z_0_xxyyy_xxx, g_z_0_z_0_xxyyy_xxy, g_z_0_z_0_xxyyy_xxz, g_z_0_z_0_xxyyy_xy, g_z_0_z_0_xxyyy_xyy, g_z_0_z_0_xxyyy_xyz, g_z_0_z_0_xxyyy_xz, g_z_0_z_0_xxyyy_xzz, g_z_0_z_0_xxyyy_yy, g_z_0_z_0_xxyyy_yz, g_z_0_z_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxyyy_xx[k] = -g_z_0_z_0_xxyyy_xx[k] * ab_x + g_z_0_z_0_xxyyy_xxx[k];

                g_z_0_z_0_xxxyyy_xy[k] = -g_z_0_z_0_xxyyy_xy[k] * ab_x + g_z_0_z_0_xxyyy_xxy[k];

                g_z_0_z_0_xxxyyy_xz[k] = -g_z_0_z_0_xxyyy_xz[k] * ab_x + g_z_0_z_0_xxyyy_xxz[k];

                g_z_0_z_0_xxxyyy_yy[k] = -g_z_0_z_0_xxyyy_yy[k] * ab_x + g_z_0_z_0_xxyyy_xyy[k];

                g_z_0_z_0_xxxyyy_yz[k] = -g_z_0_z_0_xxyyy_yz[k] * ab_x + g_z_0_z_0_xxyyy_xyz[k];

                g_z_0_z_0_xxxyyy_zz[k] = -g_z_0_z_0_xxyyy_zz[k] * ab_x + g_z_0_z_0_xxyyy_xzz[k];
            }

            /// Set up 1386-1392 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxyyz_xx = cbuffer.data(id_geom_1010_off + 1386 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyz_xy = cbuffer.data(id_geom_1010_off + 1387 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyz_xz = cbuffer.data(id_geom_1010_off + 1388 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyz_yy = cbuffer.data(id_geom_1010_off + 1389 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyz_yz = cbuffer.data(id_geom_1010_off + 1390 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyz_zz = cbuffer.data(id_geom_1010_off + 1391 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxyyz_xx, g_z_0_z_0_xxxyyz_xy, g_z_0_z_0_xxxyyz_xz, g_z_0_z_0_xxxyyz_yy, g_z_0_z_0_xxxyyz_yz, g_z_0_z_0_xxxyyz_zz, g_z_0_z_0_xxyyz_xx, g_z_0_z_0_xxyyz_xxx, g_z_0_z_0_xxyyz_xxy, g_z_0_z_0_xxyyz_xxz, g_z_0_z_0_xxyyz_xy, g_z_0_z_0_xxyyz_xyy, g_z_0_z_0_xxyyz_xyz, g_z_0_z_0_xxyyz_xz, g_z_0_z_0_xxyyz_xzz, g_z_0_z_0_xxyyz_yy, g_z_0_z_0_xxyyz_yz, g_z_0_z_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxyyz_xx[k] = -g_z_0_z_0_xxyyz_xx[k] * ab_x + g_z_0_z_0_xxyyz_xxx[k];

                g_z_0_z_0_xxxyyz_xy[k] = -g_z_0_z_0_xxyyz_xy[k] * ab_x + g_z_0_z_0_xxyyz_xxy[k];

                g_z_0_z_0_xxxyyz_xz[k] = -g_z_0_z_0_xxyyz_xz[k] * ab_x + g_z_0_z_0_xxyyz_xxz[k];

                g_z_0_z_0_xxxyyz_yy[k] = -g_z_0_z_0_xxyyz_yy[k] * ab_x + g_z_0_z_0_xxyyz_xyy[k];

                g_z_0_z_0_xxxyyz_yz[k] = -g_z_0_z_0_xxyyz_yz[k] * ab_x + g_z_0_z_0_xxyyz_xyz[k];

                g_z_0_z_0_xxxyyz_zz[k] = -g_z_0_z_0_xxyyz_zz[k] * ab_x + g_z_0_z_0_xxyyz_xzz[k];
            }

            /// Set up 1392-1398 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxyzz_xx = cbuffer.data(id_geom_1010_off + 1392 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyzz_xy = cbuffer.data(id_geom_1010_off + 1393 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyzz_xz = cbuffer.data(id_geom_1010_off + 1394 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyzz_yy = cbuffer.data(id_geom_1010_off + 1395 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyzz_yz = cbuffer.data(id_geom_1010_off + 1396 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyzz_zz = cbuffer.data(id_geom_1010_off + 1397 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxyzz_xx, g_z_0_z_0_xxxyzz_xy, g_z_0_z_0_xxxyzz_xz, g_z_0_z_0_xxxyzz_yy, g_z_0_z_0_xxxyzz_yz, g_z_0_z_0_xxxyzz_zz, g_z_0_z_0_xxyzz_xx, g_z_0_z_0_xxyzz_xxx, g_z_0_z_0_xxyzz_xxy, g_z_0_z_0_xxyzz_xxz, g_z_0_z_0_xxyzz_xy, g_z_0_z_0_xxyzz_xyy, g_z_0_z_0_xxyzz_xyz, g_z_0_z_0_xxyzz_xz, g_z_0_z_0_xxyzz_xzz, g_z_0_z_0_xxyzz_yy, g_z_0_z_0_xxyzz_yz, g_z_0_z_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxyzz_xx[k] = -g_z_0_z_0_xxyzz_xx[k] * ab_x + g_z_0_z_0_xxyzz_xxx[k];

                g_z_0_z_0_xxxyzz_xy[k] = -g_z_0_z_0_xxyzz_xy[k] * ab_x + g_z_0_z_0_xxyzz_xxy[k];

                g_z_0_z_0_xxxyzz_xz[k] = -g_z_0_z_0_xxyzz_xz[k] * ab_x + g_z_0_z_0_xxyzz_xxz[k];

                g_z_0_z_0_xxxyzz_yy[k] = -g_z_0_z_0_xxyzz_yy[k] * ab_x + g_z_0_z_0_xxyzz_xyy[k];

                g_z_0_z_0_xxxyzz_yz[k] = -g_z_0_z_0_xxyzz_yz[k] * ab_x + g_z_0_z_0_xxyzz_xyz[k];

                g_z_0_z_0_xxxyzz_zz[k] = -g_z_0_z_0_xxyzz_zz[k] * ab_x + g_z_0_z_0_xxyzz_xzz[k];
            }

            /// Set up 1398-1404 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxzzz_xx = cbuffer.data(id_geom_1010_off + 1398 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzzz_xy = cbuffer.data(id_geom_1010_off + 1399 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzzz_xz = cbuffer.data(id_geom_1010_off + 1400 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzzz_yy = cbuffer.data(id_geom_1010_off + 1401 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzzz_yz = cbuffer.data(id_geom_1010_off + 1402 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzzz_zz = cbuffer.data(id_geom_1010_off + 1403 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxzzz_xx, g_z_0_z_0_xxxzzz_xy, g_z_0_z_0_xxxzzz_xz, g_z_0_z_0_xxxzzz_yy, g_z_0_z_0_xxxzzz_yz, g_z_0_z_0_xxxzzz_zz, g_z_0_z_0_xxzzz_xx, g_z_0_z_0_xxzzz_xxx, g_z_0_z_0_xxzzz_xxy, g_z_0_z_0_xxzzz_xxz, g_z_0_z_0_xxzzz_xy, g_z_0_z_0_xxzzz_xyy, g_z_0_z_0_xxzzz_xyz, g_z_0_z_0_xxzzz_xz, g_z_0_z_0_xxzzz_xzz, g_z_0_z_0_xxzzz_yy, g_z_0_z_0_xxzzz_yz, g_z_0_z_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxzzz_xx[k] = -g_z_0_z_0_xxzzz_xx[k] * ab_x + g_z_0_z_0_xxzzz_xxx[k];

                g_z_0_z_0_xxxzzz_xy[k] = -g_z_0_z_0_xxzzz_xy[k] * ab_x + g_z_0_z_0_xxzzz_xxy[k];

                g_z_0_z_0_xxxzzz_xz[k] = -g_z_0_z_0_xxzzz_xz[k] * ab_x + g_z_0_z_0_xxzzz_xxz[k];

                g_z_0_z_0_xxxzzz_yy[k] = -g_z_0_z_0_xxzzz_yy[k] * ab_x + g_z_0_z_0_xxzzz_xyy[k];

                g_z_0_z_0_xxxzzz_yz[k] = -g_z_0_z_0_xxzzz_yz[k] * ab_x + g_z_0_z_0_xxzzz_xyz[k];

                g_z_0_z_0_xxxzzz_zz[k] = -g_z_0_z_0_xxzzz_zz[k] * ab_x + g_z_0_z_0_xxzzz_xzz[k];
            }

            /// Set up 1404-1410 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyyyy_xx = cbuffer.data(id_geom_1010_off + 1404 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyy_xy = cbuffer.data(id_geom_1010_off + 1405 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyy_xz = cbuffer.data(id_geom_1010_off + 1406 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyy_yy = cbuffer.data(id_geom_1010_off + 1407 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyy_yz = cbuffer.data(id_geom_1010_off + 1408 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyy_zz = cbuffer.data(id_geom_1010_off + 1409 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyyyy_xx, g_z_0_z_0_xxyyyy_xy, g_z_0_z_0_xxyyyy_xz, g_z_0_z_0_xxyyyy_yy, g_z_0_z_0_xxyyyy_yz, g_z_0_z_0_xxyyyy_zz, g_z_0_z_0_xyyyy_xx, g_z_0_z_0_xyyyy_xxx, g_z_0_z_0_xyyyy_xxy, g_z_0_z_0_xyyyy_xxz, g_z_0_z_0_xyyyy_xy, g_z_0_z_0_xyyyy_xyy, g_z_0_z_0_xyyyy_xyz, g_z_0_z_0_xyyyy_xz, g_z_0_z_0_xyyyy_xzz, g_z_0_z_0_xyyyy_yy, g_z_0_z_0_xyyyy_yz, g_z_0_z_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyyyy_xx[k] = -g_z_0_z_0_xyyyy_xx[k] * ab_x + g_z_0_z_0_xyyyy_xxx[k];

                g_z_0_z_0_xxyyyy_xy[k] = -g_z_0_z_0_xyyyy_xy[k] * ab_x + g_z_0_z_0_xyyyy_xxy[k];

                g_z_0_z_0_xxyyyy_xz[k] = -g_z_0_z_0_xyyyy_xz[k] * ab_x + g_z_0_z_0_xyyyy_xxz[k];

                g_z_0_z_0_xxyyyy_yy[k] = -g_z_0_z_0_xyyyy_yy[k] * ab_x + g_z_0_z_0_xyyyy_xyy[k];

                g_z_0_z_0_xxyyyy_yz[k] = -g_z_0_z_0_xyyyy_yz[k] * ab_x + g_z_0_z_0_xyyyy_xyz[k];

                g_z_0_z_0_xxyyyy_zz[k] = -g_z_0_z_0_xyyyy_zz[k] * ab_x + g_z_0_z_0_xyyyy_xzz[k];
            }

            /// Set up 1410-1416 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyyyz_xx = cbuffer.data(id_geom_1010_off + 1410 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyz_xy = cbuffer.data(id_geom_1010_off + 1411 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyz_xz = cbuffer.data(id_geom_1010_off + 1412 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyz_yy = cbuffer.data(id_geom_1010_off + 1413 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyz_yz = cbuffer.data(id_geom_1010_off + 1414 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyz_zz = cbuffer.data(id_geom_1010_off + 1415 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyyyz_xx, g_z_0_z_0_xxyyyz_xy, g_z_0_z_0_xxyyyz_xz, g_z_0_z_0_xxyyyz_yy, g_z_0_z_0_xxyyyz_yz, g_z_0_z_0_xxyyyz_zz, g_z_0_z_0_xyyyz_xx, g_z_0_z_0_xyyyz_xxx, g_z_0_z_0_xyyyz_xxy, g_z_0_z_0_xyyyz_xxz, g_z_0_z_0_xyyyz_xy, g_z_0_z_0_xyyyz_xyy, g_z_0_z_0_xyyyz_xyz, g_z_0_z_0_xyyyz_xz, g_z_0_z_0_xyyyz_xzz, g_z_0_z_0_xyyyz_yy, g_z_0_z_0_xyyyz_yz, g_z_0_z_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyyyz_xx[k] = -g_z_0_z_0_xyyyz_xx[k] * ab_x + g_z_0_z_0_xyyyz_xxx[k];

                g_z_0_z_0_xxyyyz_xy[k] = -g_z_0_z_0_xyyyz_xy[k] * ab_x + g_z_0_z_0_xyyyz_xxy[k];

                g_z_0_z_0_xxyyyz_xz[k] = -g_z_0_z_0_xyyyz_xz[k] * ab_x + g_z_0_z_0_xyyyz_xxz[k];

                g_z_0_z_0_xxyyyz_yy[k] = -g_z_0_z_0_xyyyz_yy[k] * ab_x + g_z_0_z_0_xyyyz_xyy[k];

                g_z_0_z_0_xxyyyz_yz[k] = -g_z_0_z_0_xyyyz_yz[k] * ab_x + g_z_0_z_0_xyyyz_xyz[k];

                g_z_0_z_0_xxyyyz_zz[k] = -g_z_0_z_0_xyyyz_zz[k] * ab_x + g_z_0_z_0_xyyyz_xzz[k];
            }

            /// Set up 1416-1422 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyyzz_xx = cbuffer.data(id_geom_1010_off + 1416 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyzz_xy = cbuffer.data(id_geom_1010_off + 1417 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyzz_xz = cbuffer.data(id_geom_1010_off + 1418 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyzz_yy = cbuffer.data(id_geom_1010_off + 1419 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyzz_yz = cbuffer.data(id_geom_1010_off + 1420 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyzz_zz = cbuffer.data(id_geom_1010_off + 1421 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyyzz_xx, g_z_0_z_0_xxyyzz_xy, g_z_0_z_0_xxyyzz_xz, g_z_0_z_0_xxyyzz_yy, g_z_0_z_0_xxyyzz_yz, g_z_0_z_0_xxyyzz_zz, g_z_0_z_0_xyyzz_xx, g_z_0_z_0_xyyzz_xxx, g_z_0_z_0_xyyzz_xxy, g_z_0_z_0_xyyzz_xxz, g_z_0_z_0_xyyzz_xy, g_z_0_z_0_xyyzz_xyy, g_z_0_z_0_xyyzz_xyz, g_z_0_z_0_xyyzz_xz, g_z_0_z_0_xyyzz_xzz, g_z_0_z_0_xyyzz_yy, g_z_0_z_0_xyyzz_yz, g_z_0_z_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyyzz_xx[k] = -g_z_0_z_0_xyyzz_xx[k] * ab_x + g_z_0_z_0_xyyzz_xxx[k];

                g_z_0_z_0_xxyyzz_xy[k] = -g_z_0_z_0_xyyzz_xy[k] * ab_x + g_z_0_z_0_xyyzz_xxy[k];

                g_z_0_z_0_xxyyzz_xz[k] = -g_z_0_z_0_xyyzz_xz[k] * ab_x + g_z_0_z_0_xyyzz_xxz[k];

                g_z_0_z_0_xxyyzz_yy[k] = -g_z_0_z_0_xyyzz_yy[k] * ab_x + g_z_0_z_0_xyyzz_xyy[k];

                g_z_0_z_0_xxyyzz_yz[k] = -g_z_0_z_0_xyyzz_yz[k] * ab_x + g_z_0_z_0_xyyzz_xyz[k];

                g_z_0_z_0_xxyyzz_zz[k] = -g_z_0_z_0_xyyzz_zz[k] * ab_x + g_z_0_z_0_xyyzz_xzz[k];
            }

            /// Set up 1422-1428 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyzzz_xx = cbuffer.data(id_geom_1010_off + 1422 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzzz_xy = cbuffer.data(id_geom_1010_off + 1423 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzzz_xz = cbuffer.data(id_geom_1010_off + 1424 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzzz_yy = cbuffer.data(id_geom_1010_off + 1425 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzzz_yz = cbuffer.data(id_geom_1010_off + 1426 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzzz_zz = cbuffer.data(id_geom_1010_off + 1427 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyzzz_xx, g_z_0_z_0_xxyzzz_xy, g_z_0_z_0_xxyzzz_xz, g_z_0_z_0_xxyzzz_yy, g_z_0_z_0_xxyzzz_yz, g_z_0_z_0_xxyzzz_zz, g_z_0_z_0_xyzzz_xx, g_z_0_z_0_xyzzz_xxx, g_z_0_z_0_xyzzz_xxy, g_z_0_z_0_xyzzz_xxz, g_z_0_z_0_xyzzz_xy, g_z_0_z_0_xyzzz_xyy, g_z_0_z_0_xyzzz_xyz, g_z_0_z_0_xyzzz_xz, g_z_0_z_0_xyzzz_xzz, g_z_0_z_0_xyzzz_yy, g_z_0_z_0_xyzzz_yz, g_z_0_z_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyzzz_xx[k] = -g_z_0_z_0_xyzzz_xx[k] * ab_x + g_z_0_z_0_xyzzz_xxx[k];

                g_z_0_z_0_xxyzzz_xy[k] = -g_z_0_z_0_xyzzz_xy[k] * ab_x + g_z_0_z_0_xyzzz_xxy[k];

                g_z_0_z_0_xxyzzz_xz[k] = -g_z_0_z_0_xyzzz_xz[k] * ab_x + g_z_0_z_0_xyzzz_xxz[k];

                g_z_0_z_0_xxyzzz_yy[k] = -g_z_0_z_0_xyzzz_yy[k] * ab_x + g_z_0_z_0_xyzzz_xyy[k];

                g_z_0_z_0_xxyzzz_yz[k] = -g_z_0_z_0_xyzzz_yz[k] * ab_x + g_z_0_z_0_xyzzz_xyz[k];

                g_z_0_z_0_xxyzzz_zz[k] = -g_z_0_z_0_xyzzz_zz[k] * ab_x + g_z_0_z_0_xyzzz_xzz[k];
            }

            /// Set up 1428-1434 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxzzzz_xx = cbuffer.data(id_geom_1010_off + 1428 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzzz_xy = cbuffer.data(id_geom_1010_off + 1429 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzzz_xz = cbuffer.data(id_geom_1010_off + 1430 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzzz_yy = cbuffer.data(id_geom_1010_off + 1431 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzzz_yz = cbuffer.data(id_geom_1010_off + 1432 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzzz_zz = cbuffer.data(id_geom_1010_off + 1433 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxzzzz_xx, g_z_0_z_0_xxzzzz_xy, g_z_0_z_0_xxzzzz_xz, g_z_0_z_0_xxzzzz_yy, g_z_0_z_0_xxzzzz_yz, g_z_0_z_0_xxzzzz_zz, g_z_0_z_0_xzzzz_xx, g_z_0_z_0_xzzzz_xxx, g_z_0_z_0_xzzzz_xxy, g_z_0_z_0_xzzzz_xxz, g_z_0_z_0_xzzzz_xy, g_z_0_z_0_xzzzz_xyy, g_z_0_z_0_xzzzz_xyz, g_z_0_z_0_xzzzz_xz, g_z_0_z_0_xzzzz_xzz, g_z_0_z_0_xzzzz_yy, g_z_0_z_0_xzzzz_yz, g_z_0_z_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxzzzz_xx[k] = -g_z_0_z_0_xzzzz_xx[k] * ab_x + g_z_0_z_0_xzzzz_xxx[k];

                g_z_0_z_0_xxzzzz_xy[k] = -g_z_0_z_0_xzzzz_xy[k] * ab_x + g_z_0_z_0_xzzzz_xxy[k];

                g_z_0_z_0_xxzzzz_xz[k] = -g_z_0_z_0_xzzzz_xz[k] * ab_x + g_z_0_z_0_xzzzz_xxz[k];

                g_z_0_z_0_xxzzzz_yy[k] = -g_z_0_z_0_xzzzz_yy[k] * ab_x + g_z_0_z_0_xzzzz_xyy[k];

                g_z_0_z_0_xxzzzz_yz[k] = -g_z_0_z_0_xzzzz_yz[k] * ab_x + g_z_0_z_0_xzzzz_xyz[k];

                g_z_0_z_0_xxzzzz_zz[k] = -g_z_0_z_0_xzzzz_zz[k] * ab_x + g_z_0_z_0_xzzzz_xzz[k];
            }

            /// Set up 1434-1440 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyyyy_xx = cbuffer.data(id_geom_1010_off + 1434 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyy_xy = cbuffer.data(id_geom_1010_off + 1435 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyy_xz = cbuffer.data(id_geom_1010_off + 1436 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyy_yy = cbuffer.data(id_geom_1010_off + 1437 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyy_yz = cbuffer.data(id_geom_1010_off + 1438 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyy_zz = cbuffer.data(id_geom_1010_off + 1439 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyyyy_xx, g_z_0_z_0_xyyyyy_xy, g_z_0_z_0_xyyyyy_xz, g_z_0_z_0_xyyyyy_yy, g_z_0_z_0_xyyyyy_yz, g_z_0_z_0_xyyyyy_zz, g_z_0_z_0_yyyyy_xx, g_z_0_z_0_yyyyy_xxx, g_z_0_z_0_yyyyy_xxy, g_z_0_z_0_yyyyy_xxz, g_z_0_z_0_yyyyy_xy, g_z_0_z_0_yyyyy_xyy, g_z_0_z_0_yyyyy_xyz, g_z_0_z_0_yyyyy_xz, g_z_0_z_0_yyyyy_xzz, g_z_0_z_0_yyyyy_yy, g_z_0_z_0_yyyyy_yz, g_z_0_z_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyyyy_xx[k] = -g_z_0_z_0_yyyyy_xx[k] * ab_x + g_z_0_z_0_yyyyy_xxx[k];

                g_z_0_z_0_xyyyyy_xy[k] = -g_z_0_z_0_yyyyy_xy[k] * ab_x + g_z_0_z_0_yyyyy_xxy[k];

                g_z_0_z_0_xyyyyy_xz[k] = -g_z_0_z_0_yyyyy_xz[k] * ab_x + g_z_0_z_0_yyyyy_xxz[k];

                g_z_0_z_0_xyyyyy_yy[k] = -g_z_0_z_0_yyyyy_yy[k] * ab_x + g_z_0_z_0_yyyyy_xyy[k];

                g_z_0_z_0_xyyyyy_yz[k] = -g_z_0_z_0_yyyyy_yz[k] * ab_x + g_z_0_z_0_yyyyy_xyz[k];

                g_z_0_z_0_xyyyyy_zz[k] = -g_z_0_z_0_yyyyy_zz[k] * ab_x + g_z_0_z_0_yyyyy_xzz[k];
            }

            /// Set up 1440-1446 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyyyz_xx = cbuffer.data(id_geom_1010_off + 1440 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyz_xy = cbuffer.data(id_geom_1010_off + 1441 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyz_xz = cbuffer.data(id_geom_1010_off + 1442 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyz_yy = cbuffer.data(id_geom_1010_off + 1443 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyz_yz = cbuffer.data(id_geom_1010_off + 1444 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyz_zz = cbuffer.data(id_geom_1010_off + 1445 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyyyz_xx, g_z_0_z_0_xyyyyz_xy, g_z_0_z_0_xyyyyz_xz, g_z_0_z_0_xyyyyz_yy, g_z_0_z_0_xyyyyz_yz, g_z_0_z_0_xyyyyz_zz, g_z_0_z_0_yyyyz_xx, g_z_0_z_0_yyyyz_xxx, g_z_0_z_0_yyyyz_xxy, g_z_0_z_0_yyyyz_xxz, g_z_0_z_0_yyyyz_xy, g_z_0_z_0_yyyyz_xyy, g_z_0_z_0_yyyyz_xyz, g_z_0_z_0_yyyyz_xz, g_z_0_z_0_yyyyz_xzz, g_z_0_z_0_yyyyz_yy, g_z_0_z_0_yyyyz_yz, g_z_0_z_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyyyz_xx[k] = -g_z_0_z_0_yyyyz_xx[k] * ab_x + g_z_0_z_0_yyyyz_xxx[k];

                g_z_0_z_0_xyyyyz_xy[k] = -g_z_0_z_0_yyyyz_xy[k] * ab_x + g_z_0_z_0_yyyyz_xxy[k];

                g_z_0_z_0_xyyyyz_xz[k] = -g_z_0_z_0_yyyyz_xz[k] * ab_x + g_z_0_z_0_yyyyz_xxz[k];

                g_z_0_z_0_xyyyyz_yy[k] = -g_z_0_z_0_yyyyz_yy[k] * ab_x + g_z_0_z_0_yyyyz_xyy[k];

                g_z_0_z_0_xyyyyz_yz[k] = -g_z_0_z_0_yyyyz_yz[k] * ab_x + g_z_0_z_0_yyyyz_xyz[k];

                g_z_0_z_0_xyyyyz_zz[k] = -g_z_0_z_0_yyyyz_zz[k] * ab_x + g_z_0_z_0_yyyyz_xzz[k];
            }

            /// Set up 1446-1452 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyyzz_xx = cbuffer.data(id_geom_1010_off + 1446 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyzz_xy = cbuffer.data(id_geom_1010_off + 1447 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyzz_xz = cbuffer.data(id_geom_1010_off + 1448 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyzz_yy = cbuffer.data(id_geom_1010_off + 1449 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyzz_yz = cbuffer.data(id_geom_1010_off + 1450 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyzz_zz = cbuffer.data(id_geom_1010_off + 1451 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyyzz_xx, g_z_0_z_0_xyyyzz_xy, g_z_0_z_0_xyyyzz_xz, g_z_0_z_0_xyyyzz_yy, g_z_0_z_0_xyyyzz_yz, g_z_0_z_0_xyyyzz_zz, g_z_0_z_0_yyyzz_xx, g_z_0_z_0_yyyzz_xxx, g_z_0_z_0_yyyzz_xxy, g_z_0_z_0_yyyzz_xxz, g_z_0_z_0_yyyzz_xy, g_z_0_z_0_yyyzz_xyy, g_z_0_z_0_yyyzz_xyz, g_z_0_z_0_yyyzz_xz, g_z_0_z_0_yyyzz_xzz, g_z_0_z_0_yyyzz_yy, g_z_0_z_0_yyyzz_yz, g_z_0_z_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyyzz_xx[k] = -g_z_0_z_0_yyyzz_xx[k] * ab_x + g_z_0_z_0_yyyzz_xxx[k];

                g_z_0_z_0_xyyyzz_xy[k] = -g_z_0_z_0_yyyzz_xy[k] * ab_x + g_z_0_z_0_yyyzz_xxy[k];

                g_z_0_z_0_xyyyzz_xz[k] = -g_z_0_z_0_yyyzz_xz[k] * ab_x + g_z_0_z_0_yyyzz_xxz[k];

                g_z_0_z_0_xyyyzz_yy[k] = -g_z_0_z_0_yyyzz_yy[k] * ab_x + g_z_0_z_0_yyyzz_xyy[k];

                g_z_0_z_0_xyyyzz_yz[k] = -g_z_0_z_0_yyyzz_yz[k] * ab_x + g_z_0_z_0_yyyzz_xyz[k];

                g_z_0_z_0_xyyyzz_zz[k] = -g_z_0_z_0_yyyzz_zz[k] * ab_x + g_z_0_z_0_yyyzz_xzz[k];
            }

            /// Set up 1452-1458 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyzzz_xx = cbuffer.data(id_geom_1010_off + 1452 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzzz_xy = cbuffer.data(id_geom_1010_off + 1453 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzzz_xz = cbuffer.data(id_geom_1010_off + 1454 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzzz_yy = cbuffer.data(id_geom_1010_off + 1455 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzzz_yz = cbuffer.data(id_geom_1010_off + 1456 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzzz_zz = cbuffer.data(id_geom_1010_off + 1457 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyzzz_xx, g_z_0_z_0_xyyzzz_xy, g_z_0_z_0_xyyzzz_xz, g_z_0_z_0_xyyzzz_yy, g_z_0_z_0_xyyzzz_yz, g_z_0_z_0_xyyzzz_zz, g_z_0_z_0_yyzzz_xx, g_z_0_z_0_yyzzz_xxx, g_z_0_z_0_yyzzz_xxy, g_z_0_z_0_yyzzz_xxz, g_z_0_z_0_yyzzz_xy, g_z_0_z_0_yyzzz_xyy, g_z_0_z_0_yyzzz_xyz, g_z_0_z_0_yyzzz_xz, g_z_0_z_0_yyzzz_xzz, g_z_0_z_0_yyzzz_yy, g_z_0_z_0_yyzzz_yz, g_z_0_z_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyzzz_xx[k] = -g_z_0_z_0_yyzzz_xx[k] * ab_x + g_z_0_z_0_yyzzz_xxx[k];

                g_z_0_z_0_xyyzzz_xy[k] = -g_z_0_z_0_yyzzz_xy[k] * ab_x + g_z_0_z_0_yyzzz_xxy[k];

                g_z_0_z_0_xyyzzz_xz[k] = -g_z_0_z_0_yyzzz_xz[k] * ab_x + g_z_0_z_0_yyzzz_xxz[k];

                g_z_0_z_0_xyyzzz_yy[k] = -g_z_0_z_0_yyzzz_yy[k] * ab_x + g_z_0_z_0_yyzzz_xyy[k];

                g_z_0_z_0_xyyzzz_yz[k] = -g_z_0_z_0_yyzzz_yz[k] * ab_x + g_z_0_z_0_yyzzz_xyz[k];

                g_z_0_z_0_xyyzzz_zz[k] = -g_z_0_z_0_yyzzz_zz[k] * ab_x + g_z_0_z_0_yyzzz_xzz[k];
            }

            /// Set up 1458-1464 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyzzzz_xx = cbuffer.data(id_geom_1010_off + 1458 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzzz_xy = cbuffer.data(id_geom_1010_off + 1459 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzzz_xz = cbuffer.data(id_geom_1010_off + 1460 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzzz_yy = cbuffer.data(id_geom_1010_off + 1461 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzzz_yz = cbuffer.data(id_geom_1010_off + 1462 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzzz_zz = cbuffer.data(id_geom_1010_off + 1463 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyzzzz_xx, g_z_0_z_0_xyzzzz_xy, g_z_0_z_0_xyzzzz_xz, g_z_0_z_0_xyzzzz_yy, g_z_0_z_0_xyzzzz_yz, g_z_0_z_0_xyzzzz_zz, g_z_0_z_0_yzzzz_xx, g_z_0_z_0_yzzzz_xxx, g_z_0_z_0_yzzzz_xxy, g_z_0_z_0_yzzzz_xxz, g_z_0_z_0_yzzzz_xy, g_z_0_z_0_yzzzz_xyy, g_z_0_z_0_yzzzz_xyz, g_z_0_z_0_yzzzz_xz, g_z_0_z_0_yzzzz_xzz, g_z_0_z_0_yzzzz_yy, g_z_0_z_0_yzzzz_yz, g_z_0_z_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyzzzz_xx[k] = -g_z_0_z_0_yzzzz_xx[k] * ab_x + g_z_0_z_0_yzzzz_xxx[k];

                g_z_0_z_0_xyzzzz_xy[k] = -g_z_0_z_0_yzzzz_xy[k] * ab_x + g_z_0_z_0_yzzzz_xxy[k];

                g_z_0_z_0_xyzzzz_xz[k] = -g_z_0_z_0_yzzzz_xz[k] * ab_x + g_z_0_z_0_yzzzz_xxz[k];

                g_z_0_z_0_xyzzzz_yy[k] = -g_z_0_z_0_yzzzz_yy[k] * ab_x + g_z_0_z_0_yzzzz_xyy[k];

                g_z_0_z_0_xyzzzz_yz[k] = -g_z_0_z_0_yzzzz_yz[k] * ab_x + g_z_0_z_0_yzzzz_xyz[k];

                g_z_0_z_0_xyzzzz_zz[k] = -g_z_0_z_0_yzzzz_zz[k] * ab_x + g_z_0_z_0_yzzzz_xzz[k];
            }

            /// Set up 1464-1470 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xzzzzz_xx = cbuffer.data(id_geom_1010_off + 1464 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzzz_xy = cbuffer.data(id_geom_1010_off + 1465 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzzz_xz = cbuffer.data(id_geom_1010_off + 1466 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzzz_yy = cbuffer.data(id_geom_1010_off + 1467 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzzz_yz = cbuffer.data(id_geom_1010_off + 1468 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzzz_zz = cbuffer.data(id_geom_1010_off + 1469 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xzzzzz_xx, g_z_0_z_0_xzzzzz_xy, g_z_0_z_0_xzzzzz_xz, g_z_0_z_0_xzzzzz_yy, g_z_0_z_0_xzzzzz_yz, g_z_0_z_0_xzzzzz_zz, g_z_0_z_0_zzzzz_xx, g_z_0_z_0_zzzzz_xxx, g_z_0_z_0_zzzzz_xxy, g_z_0_z_0_zzzzz_xxz, g_z_0_z_0_zzzzz_xy, g_z_0_z_0_zzzzz_xyy, g_z_0_z_0_zzzzz_xyz, g_z_0_z_0_zzzzz_xz, g_z_0_z_0_zzzzz_xzz, g_z_0_z_0_zzzzz_yy, g_z_0_z_0_zzzzz_yz, g_z_0_z_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xzzzzz_xx[k] = -g_z_0_z_0_zzzzz_xx[k] * ab_x + g_z_0_z_0_zzzzz_xxx[k];

                g_z_0_z_0_xzzzzz_xy[k] = -g_z_0_z_0_zzzzz_xy[k] * ab_x + g_z_0_z_0_zzzzz_xxy[k];

                g_z_0_z_0_xzzzzz_xz[k] = -g_z_0_z_0_zzzzz_xz[k] * ab_x + g_z_0_z_0_zzzzz_xxz[k];

                g_z_0_z_0_xzzzzz_yy[k] = -g_z_0_z_0_zzzzz_yy[k] * ab_x + g_z_0_z_0_zzzzz_xyy[k];

                g_z_0_z_0_xzzzzz_yz[k] = -g_z_0_z_0_zzzzz_yz[k] * ab_x + g_z_0_z_0_zzzzz_xyz[k];

                g_z_0_z_0_xzzzzz_zz[k] = -g_z_0_z_0_zzzzz_zz[k] * ab_x + g_z_0_z_0_zzzzz_xzz[k];
            }

            /// Set up 1470-1476 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyyyy_xx = cbuffer.data(id_geom_1010_off + 1470 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyy_xy = cbuffer.data(id_geom_1010_off + 1471 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyy_xz = cbuffer.data(id_geom_1010_off + 1472 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyy_yy = cbuffer.data(id_geom_1010_off + 1473 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyy_yz = cbuffer.data(id_geom_1010_off + 1474 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyy_zz = cbuffer.data(id_geom_1010_off + 1475 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyyy_xx, g_z_0_z_0_yyyyy_xxy, g_z_0_z_0_yyyyy_xy, g_z_0_z_0_yyyyy_xyy, g_z_0_z_0_yyyyy_xyz, g_z_0_z_0_yyyyy_xz, g_z_0_z_0_yyyyy_yy, g_z_0_z_0_yyyyy_yyy, g_z_0_z_0_yyyyy_yyz, g_z_0_z_0_yyyyy_yz, g_z_0_z_0_yyyyy_yzz, g_z_0_z_0_yyyyy_zz, g_z_0_z_0_yyyyyy_xx, g_z_0_z_0_yyyyyy_xy, g_z_0_z_0_yyyyyy_xz, g_z_0_z_0_yyyyyy_yy, g_z_0_z_0_yyyyyy_yz, g_z_0_z_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyyyy_xx[k] = -g_z_0_z_0_yyyyy_xx[k] * ab_y + g_z_0_z_0_yyyyy_xxy[k];

                g_z_0_z_0_yyyyyy_xy[k] = -g_z_0_z_0_yyyyy_xy[k] * ab_y + g_z_0_z_0_yyyyy_xyy[k];

                g_z_0_z_0_yyyyyy_xz[k] = -g_z_0_z_0_yyyyy_xz[k] * ab_y + g_z_0_z_0_yyyyy_xyz[k];

                g_z_0_z_0_yyyyyy_yy[k] = -g_z_0_z_0_yyyyy_yy[k] * ab_y + g_z_0_z_0_yyyyy_yyy[k];

                g_z_0_z_0_yyyyyy_yz[k] = -g_z_0_z_0_yyyyy_yz[k] * ab_y + g_z_0_z_0_yyyyy_yyz[k];

                g_z_0_z_0_yyyyyy_zz[k] = -g_z_0_z_0_yyyyy_zz[k] * ab_y + g_z_0_z_0_yyyyy_yzz[k];
            }

            /// Set up 1476-1482 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyyyz_xx = cbuffer.data(id_geom_1010_off + 1476 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyz_xy = cbuffer.data(id_geom_1010_off + 1477 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyz_xz = cbuffer.data(id_geom_1010_off + 1478 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyz_yy = cbuffer.data(id_geom_1010_off + 1479 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyz_yz = cbuffer.data(id_geom_1010_off + 1480 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyz_zz = cbuffer.data(id_geom_1010_off + 1481 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyyyz_xx, g_z_0_z_0_yyyyyz_xy, g_z_0_z_0_yyyyyz_xz, g_z_0_z_0_yyyyyz_yy, g_z_0_z_0_yyyyyz_yz, g_z_0_z_0_yyyyyz_zz, g_z_0_z_0_yyyyz_xx, g_z_0_z_0_yyyyz_xxy, g_z_0_z_0_yyyyz_xy, g_z_0_z_0_yyyyz_xyy, g_z_0_z_0_yyyyz_xyz, g_z_0_z_0_yyyyz_xz, g_z_0_z_0_yyyyz_yy, g_z_0_z_0_yyyyz_yyy, g_z_0_z_0_yyyyz_yyz, g_z_0_z_0_yyyyz_yz, g_z_0_z_0_yyyyz_yzz, g_z_0_z_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyyyz_xx[k] = -g_z_0_z_0_yyyyz_xx[k] * ab_y + g_z_0_z_0_yyyyz_xxy[k];

                g_z_0_z_0_yyyyyz_xy[k] = -g_z_0_z_0_yyyyz_xy[k] * ab_y + g_z_0_z_0_yyyyz_xyy[k];

                g_z_0_z_0_yyyyyz_xz[k] = -g_z_0_z_0_yyyyz_xz[k] * ab_y + g_z_0_z_0_yyyyz_xyz[k];

                g_z_0_z_0_yyyyyz_yy[k] = -g_z_0_z_0_yyyyz_yy[k] * ab_y + g_z_0_z_0_yyyyz_yyy[k];

                g_z_0_z_0_yyyyyz_yz[k] = -g_z_0_z_0_yyyyz_yz[k] * ab_y + g_z_0_z_0_yyyyz_yyz[k];

                g_z_0_z_0_yyyyyz_zz[k] = -g_z_0_z_0_yyyyz_zz[k] * ab_y + g_z_0_z_0_yyyyz_yzz[k];
            }

            /// Set up 1482-1488 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyyzz_xx = cbuffer.data(id_geom_1010_off + 1482 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyzz_xy = cbuffer.data(id_geom_1010_off + 1483 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyzz_xz = cbuffer.data(id_geom_1010_off + 1484 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyzz_yy = cbuffer.data(id_geom_1010_off + 1485 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyzz_yz = cbuffer.data(id_geom_1010_off + 1486 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyzz_zz = cbuffer.data(id_geom_1010_off + 1487 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyyzz_xx, g_z_0_z_0_yyyyzz_xy, g_z_0_z_0_yyyyzz_xz, g_z_0_z_0_yyyyzz_yy, g_z_0_z_0_yyyyzz_yz, g_z_0_z_0_yyyyzz_zz, g_z_0_z_0_yyyzz_xx, g_z_0_z_0_yyyzz_xxy, g_z_0_z_0_yyyzz_xy, g_z_0_z_0_yyyzz_xyy, g_z_0_z_0_yyyzz_xyz, g_z_0_z_0_yyyzz_xz, g_z_0_z_0_yyyzz_yy, g_z_0_z_0_yyyzz_yyy, g_z_0_z_0_yyyzz_yyz, g_z_0_z_0_yyyzz_yz, g_z_0_z_0_yyyzz_yzz, g_z_0_z_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyyzz_xx[k] = -g_z_0_z_0_yyyzz_xx[k] * ab_y + g_z_0_z_0_yyyzz_xxy[k];

                g_z_0_z_0_yyyyzz_xy[k] = -g_z_0_z_0_yyyzz_xy[k] * ab_y + g_z_0_z_0_yyyzz_xyy[k];

                g_z_0_z_0_yyyyzz_xz[k] = -g_z_0_z_0_yyyzz_xz[k] * ab_y + g_z_0_z_0_yyyzz_xyz[k];

                g_z_0_z_0_yyyyzz_yy[k] = -g_z_0_z_0_yyyzz_yy[k] * ab_y + g_z_0_z_0_yyyzz_yyy[k];

                g_z_0_z_0_yyyyzz_yz[k] = -g_z_0_z_0_yyyzz_yz[k] * ab_y + g_z_0_z_0_yyyzz_yyz[k];

                g_z_0_z_0_yyyyzz_zz[k] = -g_z_0_z_0_yyyzz_zz[k] * ab_y + g_z_0_z_0_yyyzz_yzz[k];
            }

            /// Set up 1488-1494 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyzzz_xx = cbuffer.data(id_geom_1010_off + 1488 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzzz_xy = cbuffer.data(id_geom_1010_off + 1489 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzzz_xz = cbuffer.data(id_geom_1010_off + 1490 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzzz_yy = cbuffer.data(id_geom_1010_off + 1491 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzzz_yz = cbuffer.data(id_geom_1010_off + 1492 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzzz_zz = cbuffer.data(id_geom_1010_off + 1493 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyzzz_xx, g_z_0_z_0_yyyzzz_xy, g_z_0_z_0_yyyzzz_xz, g_z_0_z_0_yyyzzz_yy, g_z_0_z_0_yyyzzz_yz, g_z_0_z_0_yyyzzz_zz, g_z_0_z_0_yyzzz_xx, g_z_0_z_0_yyzzz_xxy, g_z_0_z_0_yyzzz_xy, g_z_0_z_0_yyzzz_xyy, g_z_0_z_0_yyzzz_xyz, g_z_0_z_0_yyzzz_xz, g_z_0_z_0_yyzzz_yy, g_z_0_z_0_yyzzz_yyy, g_z_0_z_0_yyzzz_yyz, g_z_0_z_0_yyzzz_yz, g_z_0_z_0_yyzzz_yzz, g_z_0_z_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyzzz_xx[k] = -g_z_0_z_0_yyzzz_xx[k] * ab_y + g_z_0_z_0_yyzzz_xxy[k];

                g_z_0_z_0_yyyzzz_xy[k] = -g_z_0_z_0_yyzzz_xy[k] * ab_y + g_z_0_z_0_yyzzz_xyy[k];

                g_z_0_z_0_yyyzzz_xz[k] = -g_z_0_z_0_yyzzz_xz[k] * ab_y + g_z_0_z_0_yyzzz_xyz[k];

                g_z_0_z_0_yyyzzz_yy[k] = -g_z_0_z_0_yyzzz_yy[k] * ab_y + g_z_0_z_0_yyzzz_yyy[k];

                g_z_0_z_0_yyyzzz_yz[k] = -g_z_0_z_0_yyzzz_yz[k] * ab_y + g_z_0_z_0_yyzzz_yyz[k];

                g_z_0_z_0_yyyzzz_zz[k] = -g_z_0_z_0_yyzzz_zz[k] * ab_y + g_z_0_z_0_yyzzz_yzz[k];
            }

            /// Set up 1494-1500 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyzzzz_xx = cbuffer.data(id_geom_1010_off + 1494 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzzz_xy = cbuffer.data(id_geom_1010_off + 1495 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzzz_xz = cbuffer.data(id_geom_1010_off + 1496 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzzz_yy = cbuffer.data(id_geom_1010_off + 1497 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzzz_yz = cbuffer.data(id_geom_1010_off + 1498 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzzz_zz = cbuffer.data(id_geom_1010_off + 1499 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyzzzz_xx, g_z_0_z_0_yyzzzz_xy, g_z_0_z_0_yyzzzz_xz, g_z_0_z_0_yyzzzz_yy, g_z_0_z_0_yyzzzz_yz, g_z_0_z_0_yyzzzz_zz, g_z_0_z_0_yzzzz_xx, g_z_0_z_0_yzzzz_xxy, g_z_0_z_0_yzzzz_xy, g_z_0_z_0_yzzzz_xyy, g_z_0_z_0_yzzzz_xyz, g_z_0_z_0_yzzzz_xz, g_z_0_z_0_yzzzz_yy, g_z_0_z_0_yzzzz_yyy, g_z_0_z_0_yzzzz_yyz, g_z_0_z_0_yzzzz_yz, g_z_0_z_0_yzzzz_yzz, g_z_0_z_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyzzzz_xx[k] = -g_z_0_z_0_yzzzz_xx[k] * ab_y + g_z_0_z_0_yzzzz_xxy[k];

                g_z_0_z_0_yyzzzz_xy[k] = -g_z_0_z_0_yzzzz_xy[k] * ab_y + g_z_0_z_0_yzzzz_xyy[k];

                g_z_0_z_0_yyzzzz_xz[k] = -g_z_0_z_0_yzzzz_xz[k] * ab_y + g_z_0_z_0_yzzzz_xyz[k];

                g_z_0_z_0_yyzzzz_yy[k] = -g_z_0_z_0_yzzzz_yy[k] * ab_y + g_z_0_z_0_yzzzz_yyy[k];

                g_z_0_z_0_yyzzzz_yz[k] = -g_z_0_z_0_yzzzz_yz[k] * ab_y + g_z_0_z_0_yzzzz_yyz[k];

                g_z_0_z_0_yyzzzz_zz[k] = -g_z_0_z_0_yzzzz_zz[k] * ab_y + g_z_0_z_0_yzzzz_yzz[k];
            }

            /// Set up 1500-1506 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yzzzzz_xx = cbuffer.data(id_geom_1010_off + 1500 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzzz_xy = cbuffer.data(id_geom_1010_off + 1501 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzzz_xz = cbuffer.data(id_geom_1010_off + 1502 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzzz_yy = cbuffer.data(id_geom_1010_off + 1503 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzzz_yz = cbuffer.data(id_geom_1010_off + 1504 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzzz_zz = cbuffer.data(id_geom_1010_off + 1505 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yzzzzz_xx, g_z_0_z_0_yzzzzz_xy, g_z_0_z_0_yzzzzz_xz, g_z_0_z_0_yzzzzz_yy, g_z_0_z_0_yzzzzz_yz, g_z_0_z_0_yzzzzz_zz, g_z_0_z_0_zzzzz_xx, g_z_0_z_0_zzzzz_xxy, g_z_0_z_0_zzzzz_xy, g_z_0_z_0_zzzzz_xyy, g_z_0_z_0_zzzzz_xyz, g_z_0_z_0_zzzzz_xz, g_z_0_z_0_zzzzz_yy, g_z_0_z_0_zzzzz_yyy, g_z_0_z_0_zzzzz_yyz, g_z_0_z_0_zzzzz_yz, g_z_0_z_0_zzzzz_yzz, g_z_0_z_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yzzzzz_xx[k] = -g_z_0_z_0_zzzzz_xx[k] * ab_y + g_z_0_z_0_zzzzz_xxy[k];

                g_z_0_z_0_yzzzzz_xy[k] = -g_z_0_z_0_zzzzz_xy[k] * ab_y + g_z_0_z_0_zzzzz_xyy[k];

                g_z_0_z_0_yzzzzz_xz[k] = -g_z_0_z_0_zzzzz_xz[k] * ab_y + g_z_0_z_0_zzzzz_xyz[k];

                g_z_0_z_0_yzzzzz_yy[k] = -g_z_0_z_0_zzzzz_yy[k] * ab_y + g_z_0_z_0_zzzzz_yyy[k];

                g_z_0_z_0_yzzzzz_yz[k] = -g_z_0_z_0_zzzzz_yz[k] * ab_y + g_z_0_z_0_zzzzz_yyz[k];

                g_z_0_z_0_yzzzzz_zz[k] = -g_z_0_z_0_zzzzz_zz[k] * ab_y + g_z_0_z_0_zzzzz_yzz[k];
            }

            /// Set up 1506-1512 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_zzzzzz_xx = cbuffer.data(id_geom_1010_off + 1506 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzzz_xy = cbuffer.data(id_geom_1010_off + 1507 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzzz_xz = cbuffer.data(id_geom_1010_off + 1508 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzzz_yy = cbuffer.data(id_geom_1010_off + 1509 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzzz_yz = cbuffer.data(id_geom_1010_off + 1510 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzzz_zz = cbuffer.data(id_geom_1010_off + 1511 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_zzzzz_xx, g_0_0_z_0_zzzzz_xy, g_0_0_z_0_zzzzz_xz, g_0_0_z_0_zzzzz_yy, g_0_0_z_0_zzzzz_yz, g_0_0_z_0_zzzzz_zz, g_z_0_z_0_zzzzz_xx, g_z_0_z_0_zzzzz_xxz, g_z_0_z_0_zzzzz_xy, g_z_0_z_0_zzzzz_xyz, g_z_0_z_0_zzzzz_xz, g_z_0_z_0_zzzzz_xzz, g_z_0_z_0_zzzzz_yy, g_z_0_z_0_zzzzz_yyz, g_z_0_z_0_zzzzz_yz, g_z_0_z_0_zzzzz_yzz, g_z_0_z_0_zzzzz_zz, g_z_0_z_0_zzzzz_zzz, g_z_0_z_0_zzzzzz_xx, g_z_0_z_0_zzzzzz_xy, g_z_0_z_0_zzzzzz_xz, g_z_0_z_0_zzzzzz_yy, g_z_0_z_0_zzzzzz_yz, g_z_0_z_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_zzzzzz_xx[k] = -g_0_0_z_0_zzzzz_xx[k] - g_z_0_z_0_zzzzz_xx[k] * ab_z + g_z_0_z_0_zzzzz_xxz[k];

                g_z_0_z_0_zzzzzz_xy[k] = -g_0_0_z_0_zzzzz_xy[k] - g_z_0_z_0_zzzzz_xy[k] * ab_z + g_z_0_z_0_zzzzz_xyz[k];

                g_z_0_z_0_zzzzzz_xz[k] = -g_0_0_z_0_zzzzz_xz[k] - g_z_0_z_0_zzzzz_xz[k] * ab_z + g_z_0_z_0_zzzzz_xzz[k];

                g_z_0_z_0_zzzzzz_yy[k] = -g_0_0_z_0_zzzzz_yy[k] - g_z_0_z_0_zzzzz_yy[k] * ab_z + g_z_0_z_0_zzzzz_yyz[k];

                g_z_0_z_0_zzzzzz_yz[k] = -g_0_0_z_0_zzzzz_yz[k] - g_z_0_z_0_zzzzz_yz[k] * ab_z + g_z_0_z_0_zzzzz_yzz[k];

                g_z_0_z_0_zzzzzz_zz[k] = -g_0_0_z_0_zzzzz_zz[k] - g_z_0_z_0_zzzzz_zz[k] * ab_z + g_z_0_z_0_zzzzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

