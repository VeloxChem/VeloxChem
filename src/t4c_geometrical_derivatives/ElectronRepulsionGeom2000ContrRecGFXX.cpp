#include "ElectronRepulsionGeom2000ContrRecGFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_gfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_gfxx,
                                            const size_t idx_geom_10_ffxx,
                                            const size_t idx_geom_20_ffxx,
                                            const size_t idx_geom_20_fgxx,
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
            /// Set up components of auxilary buffer : FFSS

            const auto ff_geom_10_off = idx_geom_10_ffxx + i * dcomps + j;

            auto g_x_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 125 * ccomps * dcomps);

            auto g_y_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 126 * ccomps * dcomps);

            auto g_y_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 127 * ccomps * dcomps);

            auto g_y_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 128 * ccomps * dcomps);

            auto g_y_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 129 * ccomps * dcomps);

            auto g_y_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 130 * ccomps * dcomps);

            auto g_y_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 131 * ccomps * dcomps);

            auto g_y_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 132 * ccomps * dcomps);

            auto g_y_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 133 * ccomps * dcomps);

            auto g_y_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 134 * ccomps * dcomps);

            auto g_y_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 135 * ccomps * dcomps);

            auto g_y_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 136 * ccomps * dcomps);

            auto g_y_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 137 * ccomps * dcomps);

            auto g_y_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 138 * ccomps * dcomps);

            auto g_y_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 139 * ccomps * dcomps);

            auto g_y_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 140 * ccomps * dcomps);

            auto g_y_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 141 * ccomps * dcomps);

            auto g_y_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 142 * ccomps * dcomps);

            auto g_y_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 143 * ccomps * dcomps);

            auto g_y_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 144 * ccomps * dcomps);

            auto g_y_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 145 * ccomps * dcomps);

            auto g_y_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 146 * ccomps * dcomps);

            auto g_y_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 147 * ccomps * dcomps);

            auto g_y_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 148 * ccomps * dcomps);

            auto g_y_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 149 * ccomps * dcomps);

            auto g_y_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 164 * ccomps * dcomps);

            auto g_y_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 179 * ccomps * dcomps);

            auto g_y_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 180 * ccomps * dcomps);

            auto g_y_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 181 * ccomps * dcomps);

            auto g_y_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 182 * ccomps * dcomps);

            auto g_y_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 183 * ccomps * dcomps);

            auto g_y_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 184 * ccomps * dcomps);

            auto g_y_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 185 * ccomps * dcomps);

            auto g_y_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 186 * ccomps * dcomps);

            auto g_y_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 187 * ccomps * dcomps);

            auto g_y_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 188 * ccomps * dcomps);

            auto g_y_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 189 * ccomps * dcomps);

            auto g_y_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 190 * ccomps * dcomps);

            auto g_y_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 191 * ccomps * dcomps);

            auto g_y_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 192 * ccomps * dcomps);

            auto g_y_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 193 * ccomps * dcomps);

            auto g_y_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 194 * ccomps * dcomps);

            auto g_y_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 195 * ccomps * dcomps);

            auto g_y_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 196 * ccomps * dcomps);

            auto g_y_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 197 * ccomps * dcomps);

            auto g_y_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 198 * ccomps * dcomps);

            auto g_y_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 199 * ccomps * dcomps);

            auto g_z_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 200 * ccomps * dcomps);

            auto g_z_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 201 * ccomps * dcomps);

            auto g_z_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 202 * ccomps * dcomps);

            auto g_z_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 203 * ccomps * dcomps);

            auto g_z_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 204 * ccomps * dcomps);

            auto g_z_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 205 * ccomps * dcomps);

            auto g_z_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 206 * ccomps * dcomps);

            auto g_z_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 207 * ccomps * dcomps);

            auto g_z_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 208 * ccomps * dcomps);

            auto g_z_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 209 * ccomps * dcomps);

            auto g_z_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 210 * ccomps * dcomps);

            auto g_z_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 211 * ccomps * dcomps);

            auto g_z_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 212 * ccomps * dcomps);

            auto g_z_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 213 * ccomps * dcomps);

            auto g_z_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 214 * ccomps * dcomps);

            auto g_z_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 215 * ccomps * dcomps);

            auto g_z_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 216 * ccomps * dcomps);

            auto g_z_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 217 * ccomps * dcomps);

            auto g_z_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 218 * ccomps * dcomps);

            auto g_z_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 219 * ccomps * dcomps);

            auto g_z_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 220 * ccomps * dcomps);

            auto g_z_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 221 * ccomps * dcomps);

            auto g_z_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 222 * ccomps * dcomps);

            auto g_z_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 223 * ccomps * dcomps);

            auto g_z_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 224 * ccomps * dcomps);

            auto g_z_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 225 * ccomps * dcomps);

            auto g_z_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 226 * ccomps * dcomps);

            auto g_z_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 227 * ccomps * dcomps);

            auto g_z_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 228 * ccomps * dcomps);

            auto g_z_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 229 * ccomps * dcomps);

            auto g_z_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 230 * ccomps * dcomps);

            auto g_z_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 231 * ccomps * dcomps);

            auto g_z_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 232 * ccomps * dcomps);

            auto g_z_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 233 * ccomps * dcomps);

            auto g_z_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 234 * ccomps * dcomps);

            auto g_z_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 235 * ccomps * dcomps);

            auto g_z_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 236 * ccomps * dcomps);

            auto g_z_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 237 * ccomps * dcomps);

            auto g_z_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 238 * ccomps * dcomps);

            auto g_z_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 239 * ccomps * dcomps);

            auto g_z_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 240 * ccomps * dcomps);

            auto g_z_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 241 * ccomps * dcomps);

            auto g_z_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 242 * ccomps * dcomps);

            auto g_z_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 243 * ccomps * dcomps);

            auto g_z_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 244 * ccomps * dcomps);

            auto g_z_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 245 * ccomps * dcomps);

            auto g_z_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 246 * ccomps * dcomps);

            auto g_z_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 247 * ccomps * dcomps);

            auto g_z_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 248 * ccomps * dcomps);

            auto g_z_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 249 * ccomps * dcomps);

            auto g_z_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 250 * ccomps * dcomps);

            auto g_z_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 251 * ccomps * dcomps);

            auto g_z_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 252 * ccomps * dcomps);

            auto g_z_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 253 * ccomps * dcomps);

            auto g_z_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 254 * ccomps * dcomps);

            auto g_z_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 255 * ccomps * dcomps);

            auto g_z_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 256 * ccomps * dcomps);

            auto g_z_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 257 * ccomps * dcomps);

            auto g_z_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 258 * ccomps * dcomps);

            auto g_z_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 259 * ccomps * dcomps);

            auto g_z_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 260 * ccomps * dcomps);

            auto g_z_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 261 * ccomps * dcomps);

            auto g_z_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 262 * ccomps * dcomps);

            auto g_z_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 263 * ccomps * dcomps);

            auto g_z_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 264 * ccomps * dcomps);

            auto g_z_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 265 * ccomps * dcomps);

            auto g_z_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 266 * ccomps * dcomps);

            auto g_z_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 267 * ccomps * dcomps);

            auto g_z_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 268 * ccomps * dcomps);

            auto g_z_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 269 * ccomps * dcomps);

            auto g_z_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 270 * ccomps * dcomps);

            auto g_z_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 271 * ccomps * dcomps);

            auto g_z_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 272 * ccomps * dcomps);

            auto g_z_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 273 * ccomps * dcomps);

            auto g_z_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 274 * ccomps * dcomps);

            auto g_z_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 275 * ccomps * dcomps);

            auto g_z_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 276 * ccomps * dcomps);

            auto g_z_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 277 * ccomps * dcomps);

            auto g_z_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 278 * ccomps * dcomps);

            auto g_z_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 279 * ccomps * dcomps);

            auto g_z_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 280 * ccomps * dcomps);

            auto g_z_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 281 * ccomps * dcomps);

            auto g_z_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 282 * ccomps * dcomps);

            auto g_z_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 283 * ccomps * dcomps);

            auto g_z_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 284 * ccomps * dcomps);

            auto g_z_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 285 * ccomps * dcomps);

            auto g_z_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 286 * ccomps * dcomps);

            auto g_z_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 287 * ccomps * dcomps);

            auto g_z_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 288 * ccomps * dcomps);

            auto g_z_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 289 * ccomps * dcomps);

            auto g_z_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 290 * ccomps * dcomps);

            auto g_z_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 291 * ccomps * dcomps);

            auto g_z_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 292 * ccomps * dcomps);

            auto g_z_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 293 * ccomps * dcomps);

            auto g_z_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 294 * ccomps * dcomps);

            auto g_z_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 295 * ccomps * dcomps);

            auto g_z_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 296 * ccomps * dcomps);

            auto g_z_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 297 * ccomps * dcomps);

            auto g_z_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 298 * ccomps * dcomps);

            auto g_z_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 299 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FFSS

            const auto ff_geom_20_off = idx_geom_20_ffxx + i * dcomps + j;

            auto g_xx_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 199 * ccomps * dcomps);

            auto g_xz_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 200 * ccomps * dcomps);

            auto g_xz_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 201 * ccomps * dcomps);

            auto g_xz_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 202 * ccomps * dcomps);

            auto g_xz_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 203 * ccomps * dcomps);

            auto g_xz_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 204 * ccomps * dcomps);

            auto g_xz_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 205 * ccomps * dcomps);

            auto g_xz_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 206 * ccomps * dcomps);

            auto g_xz_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 207 * ccomps * dcomps);

            auto g_xz_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 208 * ccomps * dcomps);

            auto g_xz_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 209 * ccomps * dcomps);

            auto g_xz_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 210 * ccomps * dcomps);

            auto g_xz_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 211 * ccomps * dcomps);

            auto g_xz_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 212 * ccomps * dcomps);

            auto g_xz_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 213 * ccomps * dcomps);

            auto g_xz_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 214 * ccomps * dcomps);

            auto g_xz_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 215 * ccomps * dcomps);

            auto g_xz_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 218 * ccomps * dcomps);

            auto g_xz_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 221 * ccomps * dcomps);

            auto g_xz_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 223 * ccomps * dcomps);

            auto g_xz_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 224 * ccomps * dcomps);

            auto g_xz_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 227 * ccomps * dcomps);

            auto g_xz_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 230 * ccomps * dcomps);

            auto g_xz_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 233 * ccomps * dcomps);

            auto g_xz_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 236 * ccomps * dcomps);

            auto g_xz_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 239 * ccomps * dcomps);

            auto g_xz_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 242 * ccomps * dcomps);

            auto g_xz_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 245 * ccomps * dcomps);

            auto g_xz_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 248 * ccomps * dcomps);

            auto g_xz_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 269 * ccomps * dcomps);

            auto g_xz_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 270 * ccomps * dcomps);

            auto g_xz_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 271 * ccomps * dcomps);

            auto g_xz_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 272 * ccomps * dcomps);

            auto g_xz_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 273 * ccomps * dcomps);

            auto g_xz_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 274 * ccomps * dcomps);

            auto g_xz_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 275 * ccomps * dcomps);

            auto g_xz_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 276 * ccomps * dcomps);

            auto g_xz_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 277 * ccomps * dcomps);

            auto g_xz_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 278 * ccomps * dcomps);

            auto g_xz_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 279 * ccomps * dcomps);

            auto g_xz_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 280 * ccomps * dcomps);

            auto g_xz_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 281 * ccomps * dcomps);

            auto g_xz_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 282 * ccomps * dcomps);

            auto g_xz_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 283 * ccomps * dcomps);

            auto g_xz_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 284 * ccomps * dcomps);

            auto g_xz_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 285 * ccomps * dcomps);

            auto g_xz_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 286 * ccomps * dcomps);

            auto g_xz_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 287 * ccomps * dcomps);

            auto g_xz_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 288 * ccomps * dcomps);

            auto g_xz_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 289 * ccomps * dcomps);

            auto g_xz_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 290 * ccomps * dcomps);

            auto g_xz_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 291 * ccomps * dcomps);

            auto g_xz_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 292 * ccomps * dcomps);

            auto g_xz_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 293 * ccomps * dcomps);

            auto g_xz_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 294 * ccomps * dcomps);

            auto g_xz_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 295 * ccomps * dcomps);

            auto g_xz_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 296 * ccomps * dcomps);

            auto g_xz_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 297 * ccomps * dcomps);

            auto g_xz_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 298 * ccomps * dcomps);

            auto g_xz_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 299 * ccomps * dcomps);

            auto g_yy_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 300 * ccomps * dcomps);

            auto g_yy_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 301 * ccomps * dcomps);

            auto g_yy_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 302 * ccomps * dcomps);

            auto g_yy_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 303 * ccomps * dcomps);

            auto g_yy_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 304 * ccomps * dcomps);

            auto g_yy_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 305 * ccomps * dcomps);

            auto g_yy_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 306 * ccomps * dcomps);

            auto g_yy_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 307 * ccomps * dcomps);

            auto g_yy_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 308 * ccomps * dcomps);

            auto g_yy_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 309 * ccomps * dcomps);

            auto g_yy_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 310 * ccomps * dcomps);

            auto g_yy_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 311 * ccomps * dcomps);

            auto g_yy_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 312 * ccomps * dcomps);

            auto g_yy_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 313 * ccomps * dcomps);

            auto g_yy_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 314 * ccomps * dcomps);

            auto g_yy_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 315 * ccomps * dcomps);

            auto g_yy_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 316 * ccomps * dcomps);

            auto g_yy_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 317 * ccomps * dcomps);

            auto g_yy_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 318 * ccomps * dcomps);

            auto g_yy_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 319 * ccomps * dcomps);

            auto g_yy_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 320 * ccomps * dcomps);

            auto g_yy_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 321 * ccomps * dcomps);

            auto g_yy_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 322 * ccomps * dcomps);

            auto g_yy_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 323 * ccomps * dcomps);

            auto g_yy_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 326 * ccomps * dcomps);

            auto g_yy_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 329 * ccomps * dcomps);

            auto g_yy_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 332 * ccomps * dcomps);

            auto g_yy_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 335 * ccomps * dcomps);

            auto g_yy_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 336 * ccomps * dcomps);

            auto g_yy_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 337 * ccomps * dcomps);

            auto g_yy_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 338 * ccomps * dcomps);

            auto g_yy_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 339 * ccomps * dcomps);

            auto g_yy_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 340 * ccomps * dcomps);

            auto g_yy_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 341 * ccomps * dcomps);

            auto g_yy_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 342 * ccomps * dcomps);

            auto g_yy_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 343 * ccomps * dcomps);

            auto g_yy_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 344 * ccomps * dcomps);

            auto g_yy_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 345 * ccomps * dcomps);

            auto g_yy_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 346 * ccomps * dcomps);

            auto g_yy_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 347 * ccomps * dcomps);

            auto g_yy_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 348 * ccomps * dcomps);

            auto g_yy_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 349 * ccomps * dcomps);

            auto g_yy_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 350 * ccomps * dcomps);

            auto g_yy_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 351 * ccomps * dcomps);

            auto g_yy_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 352 * ccomps * dcomps);

            auto g_yy_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 353 * ccomps * dcomps);

            auto g_yy_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 354 * ccomps * dcomps);

            auto g_yy_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 355 * ccomps * dcomps);

            auto g_yy_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 356 * ccomps * dcomps);

            auto g_yy_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 357 * ccomps * dcomps);

            auto g_yy_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 358 * ccomps * dcomps);

            auto g_yy_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 359 * ccomps * dcomps);

            auto g_yy_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 360 * ccomps * dcomps);

            auto g_yy_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 361 * ccomps * dcomps);

            auto g_yy_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 362 * ccomps * dcomps);

            auto g_yy_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 363 * ccomps * dcomps);

            auto g_yy_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 364 * ccomps * dcomps);

            auto g_yy_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 365 * ccomps * dcomps);

            auto g_yy_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 366 * ccomps * dcomps);

            auto g_yy_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 367 * ccomps * dcomps);

            auto g_yy_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 368 * ccomps * dcomps);

            auto g_yy_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 369 * ccomps * dcomps);

            auto g_yy_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 370 * ccomps * dcomps);

            auto g_yy_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 371 * ccomps * dcomps);

            auto g_yy_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 372 * ccomps * dcomps);

            auto g_yy_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 373 * ccomps * dcomps);

            auto g_yy_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 374 * ccomps * dcomps);

            auto g_yy_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 375 * ccomps * dcomps);

            auto g_yy_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 376 * ccomps * dcomps);

            auto g_yy_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 377 * ccomps * dcomps);

            auto g_yy_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 378 * ccomps * dcomps);

            auto g_yy_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 379 * ccomps * dcomps);

            auto g_yy_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 380 * ccomps * dcomps);

            auto g_yy_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 381 * ccomps * dcomps);

            auto g_yy_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 382 * ccomps * dcomps);

            auto g_yy_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 383 * ccomps * dcomps);

            auto g_yy_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 384 * ccomps * dcomps);

            auto g_yy_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 385 * ccomps * dcomps);

            auto g_yy_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 386 * ccomps * dcomps);

            auto g_yy_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 387 * ccomps * dcomps);

            auto g_yy_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 388 * ccomps * dcomps);

            auto g_yy_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 389 * ccomps * dcomps);

            auto g_yy_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 390 * ccomps * dcomps);

            auto g_yy_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 391 * ccomps * dcomps);

            auto g_yy_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 392 * ccomps * dcomps);

            auto g_yy_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 393 * ccomps * dcomps);

            auto g_yy_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 394 * ccomps * dcomps);

            auto g_yy_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 395 * ccomps * dcomps);

            auto g_yy_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 396 * ccomps * dcomps);

            auto g_yy_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 397 * ccomps * dcomps);

            auto g_yy_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 398 * ccomps * dcomps);

            auto g_yy_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 399 * ccomps * dcomps);

            auto g_yz_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 400 * ccomps * dcomps);

            auto g_yz_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 401 * ccomps * dcomps);

            auto g_yz_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 402 * ccomps * dcomps);

            auto g_yz_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 403 * ccomps * dcomps);

            auto g_yz_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 404 * ccomps * dcomps);

            auto g_yz_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 405 * ccomps * dcomps);

            auto g_yz_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 406 * ccomps * dcomps);

            auto g_yz_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 407 * ccomps * dcomps);

            auto g_yz_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 408 * ccomps * dcomps);

            auto g_yz_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 409 * ccomps * dcomps);

            auto g_yz_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 410 * ccomps * dcomps);

            auto g_yz_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 411 * ccomps * dcomps);

            auto g_yz_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 412 * ccomps * dcomps);

            auto g_yz_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 413 * ccomps * dcomps);

            auto g_yz_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 414 * ccomps * dcomps);

            auto g_yz_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 415 * ccomps * dcomps);

            auto g_yz_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 416 * ccomps * dcomps);

            auto g_yz_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 417 * ccomps * dcomps);

            auto g_yz_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 418 * ccomps * dcomps);

            auto g_yz_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 419 * ccomps * dcomps);

            auto g_yz_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 420 * ccomps * dcomps);

            auto g_yz_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 421 * ccomps * dcomps);

            auto g_yz_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 422 * ccomps * dcomps);

            auto g_yz_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 423 * ccomps * dcomps);

            auto g_yz_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 424 * ccomps * dcomps);

            auto g_yz_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 425 * ccomps * dcomps);

            auto g_yz_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 426 * ccomps * dcomps);

            auto g_yz_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 427 * ccomps * dcomps);

            auto g_yz_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 428 * ccomps * dcomps);

            auto g_yz_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 429 * ccomps * dcomps);

            auto g_yz_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 430 * ccomps * dcomps);

            auto g_yz_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 431 * ccomps * dcomps);

            auto g_yz_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 432 * ccomps * dcomps);

            auto g_yz_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 433 * ccomps * dcomps);

            auto g_yz_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 434 * ccomps * dcomps);

            auto g_yz_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 435 * ccomps * dcomps);

            auto g_yz_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 436 * ccomps * dcomps);

            auto g_yz_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 437 * ccomps * dcomps);

            auto g_yz_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 438 * ccomps * dcomps);

            auto g_yz_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 439 * ccomps * dcomps);

            auto g_yz_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 440 * ccomps * dcomps);

            auto g_yz_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 441 * ccomps * dcomps);

            auto g_yz_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 442 * ccomps * dcomps);

            auto g_yz_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 443 * ccomps * dcomps);

            auto g_yz_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 444 * ccomps * dcomps);

            auto g_yz_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 445 * ccomps * dcomps);

            auto g_yz_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 446 * ccomps * dcomps);

            auto g_yz_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 447 * ccomps * dcomps);

            auto g_yz_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 448 * ccomps * dcomps);

            auto g_yz_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 449 * ccomps * dcomps);

            auto g_yz_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 450 * ccomps * dcomps);

            auto g_yz_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 451 * ccomps * dcomps);

            auto g_yz_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 452 * ccomps * dcomps);

            auto g_yz_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 453 * ccomps * dcomps);

            auto g_yz_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 454 * ccomps * dcomps);

            auto g_yz_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 455 * ccomps * dcomps);

            auto g_yz_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 456 * ccomps * dcomps);

            auto g_yz_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 457 * ccomps * dcomps);

            auto g_yz_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 458 * ccomps * dcomps);

            auto g_yz_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 459 * ccomps * dcomps);

            auto g_yz_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 460 * ccomps * dcomps);

            auto g_yz_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 461 * ccomps * dcomps);

            auto g_yz_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 462 * ccomps * dcomps);

            auto g_yz_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 463 * ccomps * dcomps);

            auto g_yz_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 464 * ccomps * dcomps);

            auto g_yz_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 465 * ccomps * dcomps);

            auto g_yz_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 466 * ccomps * dcomps);

            auto g_yz_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 467 * ccomps * dcomps);

            auto g_yz_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 468 * ccomps * dcomps);

            auto g_yz_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 469 * ccomps * dcomps);

            auto g_yz_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 470 * ccomps * dcomps);

            auto g_yz_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 471 * ccomps * dcomps);

            auto g_yz_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 472 * ccomps * dcomps);

            auto g_yz_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 473 * ccomps * dcomps);

            auto g_yz_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 474 * ccomps * dcomps);

            auto g_yz_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 475 * ccomps * dcomps);

            auto g_yz_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 476 * ccomps * dcomps);

            auto g_yz_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 477 * ccomps * dcomps);

            auto g_yz_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 478 * ccomps * dcomps);

            auto g_yz_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 479 * ccomps * dcomps);

            auto g_yz_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 480 * ccomps * dcomps);

            auto g_yz_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 481 * ccomps * dcomps);

            auto g_yz_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 482 * ccomps * dcomps);

            auto g_yz_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 483 * ccomps * dcomps);

            auto g_yz_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 484 * ccomps * dcomps);

            auto g_yz_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 485 * ccomps * dcomps);

            auto g_yz_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 486 * ccomps * dcomps);

            auto g_yz_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 487 * ccomps * dcomps);

            auto g_yz_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 488 * ccomps * dcomps);

            auto g_yz_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 489 * ccomps * dcomps);

            auto g_yz_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 490 * ccomps * dcomps);

            auto g_yz_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 491 * ccomps * dcomps);

            auto g_yz_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 492 * ccomps * dcomps);

            auto g_yz_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 493 * ccomps * dcomps);

            auto g_yz_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 494 * ccomps * dcomps);

            auto g_yz_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 495 * ccomps * dcomps);

            auto g_yz_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 496 * ccomps * dcomps);

            auto g_yz_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 497 * ccomps * dcomps);

            auto g_yz_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 498 * ccomps * dcomps);

            auto g_yz_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 499 * ccomps * dcomps);

            auto g_zz_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 500 * ccomps * dcomps);

            auto g_zz_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 501 * ccomps * dcomps);

            auto g_zz_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 502 * ccomps * dcomps);

            auto g_zz_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 503 * ccomps * dcomps);

            auto g_zz_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 504 * ccomps * dcomps);

            auto g_zz_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 505 * ccomps * dcomps);

            auto g_zz_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 506 * ccomps * dcomps);

            auto g_zz_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 507 * ccomps * dcomps);

            auto g_zz_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 508 * ccomps * dcomps);

            auto g_zz_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 509 * ccomps * dcomps);

            auto g_zz_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 510 * ccomps * dcomps);

            auto g_zz_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 511 * ccomps * dcomps);

            auto g_zz_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 512 * ccomps * dcomps);

            auto g_zz_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 513 * ccomps * dcomps);

            auto g_zz_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 514 * ccomps * dcomps);

            auto g_zz_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 515 * ccomps * dcomps);

            auto g_zz_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 516 * ccomps * dcomps);

            auto g_zz_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 517 * ccomps * dcomps);

            auto g_zz_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 518 * ccomps * dcomps);

            auto g_zz_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 519 * ccomps * dcomps);

            auto g_zz_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 520 * ccomps * dcomps);

            auto g_zz_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 521 * ccomps * dcomps);

            auto g_zz_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 522 * ccomps * dcomps);

            auto g_zz_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 523 * ccomps * dcomps);

            auto g_zz_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 524 * ccomps * dcomps);

            auto g_zz_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 525 * ccomps * dcomps);

            auto g_zz_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 526 * ccomps * dcomps);

            auto g_zz_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 527 * ccomps * dcomps);

            auto g_zz_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 528 * ccomps * dcomps);

            auto g_zz_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 529 * ccomps * dcomps);

            auto g_zz_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 530 * ccomps * dcomps);

            auto g_zz_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 531 * ccomps * dcomps);

            auto g_zz_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 532 * ccomps * dcomps);

            auto g_zz_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 533 * ccomps * dcomps);

            auto g_zz_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 534 * ccomps * dcomps);

            auto g_zz_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 535 * ccomps * dcomps);

            auto g_zz_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 536 * ccomps * dcomps);

            auto g_zz_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 537 * ccomps * dcomps);

            auto g_zz_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 538 * ccomps * dcomps);

            auto g_zz_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 539 * ccomps * dcomps);

            auto g_zz_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 540 * ccomps * dcomps);

            auto g_zz_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 541 * ccomps * dcomps);

            auto g_zz_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 542 * ccomps * dcomps);

            auto g_zz_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 543 * ccomps * dcomps);

            auto g_zz_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 544 * ccomps * dcomps);

            auto g_zz_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 545 * ccomps * dcomps);

            auto g_zz_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 546 * ccomps * dcomps);

            auto g_zz_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 547 * ccomps * dcomps);

            auto g_zz_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 548 * ccomps * dcomps);

            auto g_zz_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 549 * ccomps * dcomps);

            auto g_zz_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 550 * ccomps * dcomps);

            auto g_zz_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 551 * ccomps * dcomps);

            auto g_zz_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 552 * ccomps * dcomps);

            auto g_zz_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 553 * ccomps * dcomps);

            auto g_zz_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 554 * ccomps * dcomps);

            auto g_zz_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 555 * ccomps * dcomps);

            auto g_zz_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 556 * ccomps * dcomps);

            auto g_zz_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 557 * ccomps * dcomps);

            auto g_zz_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 558 * ccomps * dcomps);

            auto g_zz_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 559 * ccomps * dcomps);

            auto g_zz_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 560 * ccomps * dcomps);

            auto g_zz_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 561 * ccomps * dcomps);

            auto g_zz_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 562 * ccomps * dcomps);

            auto g_zz_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 563 * ccomps * dcomps);

            auto g_zz_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 564 * ccomps * dcomps);

            auto g_zz_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 565 * ccomps * dcomps);

            auto g_zz_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 566 * ccomps * dcomps);

            auto g_zz_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 567 * ccomps * dcomps);

            auto g_zz_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 568 * ccomps * dcomps);

            auto g_zz_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 569 * ccomps * dcomps);

            auto g_zz_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 570 * ccomps * dcomps);

            auto g_zz_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 571 * ccomps * dcomps);

            auto g_zz_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 572 * ccomps * dcomps);

            auto g_zz_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 573 * ccomps * dcomps);

            auto g_zz_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 574 * ccomps * dcomps);

            auto g_zz_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 575 * ccomps * dcomps);

            auto g_zz_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 576 * ccomps * dcomps);

            auto g_zz_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 577 * ccomps * dcomps);

            auto g_zz_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 578 * ccomps * dcomps);

            auto g_zz_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 579 * ccomps * dcomps);

            auto g_zz_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 580 * ccomps * dcomps);

            auto g_zz_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 581 * ccomps * dcomps);

            auto g_zz_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 582 * ccomps * dcomps);

            auto g_zz_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 583 * ccomps * dcomps);

            auto g_zz_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 584 * ccomps * dcomps);

            auto g_zz_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 585 * ccomps * dcomps);

            auto g_zz_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 586 * ccomps * dcomps);

            auto g_zz_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 587 * ccomps * dcomps);

            auto g_zz_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 588 * ccomps * dcomps);

            auto g_zz_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 589 * ccomps * dcomps);

            auto g_zz_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 590 * ccomps * dcomps);

            auto g_zz_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 591 * ccomps * dcomps);

            auto g_zz_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 592 * ccomps * dcomps);

            auto g_zz_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 593 * ccomps * dcomps);

            auto g_zz_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 594 * ccomps * dcomps);

            auto g_zz_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 595 * ccomps * dcomps);

            auto g_zz_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 596 * ccomps * dcomps);

            auto g_zz_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 597 * ccomps * dcomps);

            auto g_zz_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 598 * ccomps * dcomps);

            auto g_zz_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 599 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FGSS

            const auto fg_geom_20_off = idx_geom_20_fgxx + i * dcomps + j;

            auto g_xx_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 125 * ccomps * dcomps);

            auto g_xx_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 129 * ccomps * dcomps);

            auto g_xx_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 131 * ccomps * dcomps);

            auto g_xx_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 134 * ccomps * dcomps);

            auto g_xx_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 137 * ccomps * dcomps);

            auto g_xx_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 139 * ccomps * dcomps);

            auto g_xx_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 143 * ccomps * dcomps);

            auto g_xx_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 209 * ccomps * dcomps);

            auto g_xy_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 239 * ccomps * dcomps);

            auto g_xy_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 251 * ccomps * dcomps);

            auto g_xy_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 252 * ccomps * dcomps);

            auto g_xy_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 253 * ccomps * dcomps);

            auto g_xy_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 254 * ccomps * dcomps);

            auto g_xy_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 255 * ccomps * dcomps);

            auto g_xy_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 256 * ccomps * dcomps);

            auto g_xy_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 257 * ccomps * dcomps);

            auto g_xy_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 258 * ccomps * dcomps);

            auto g_xy_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 259 * ccomps * dcomps);

            auto g_xy_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 260 * ccomps * dcomps);

            auto g_xy_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 261 * ccomps * dcomps);

            auto g_xy_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 262 * ccomps * dcomps);

            auto g_xy_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 263 * ccomps * dcomps);

            auto g_xy_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 264 * ccomps * dcomps);

            auto g_xy_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 265 * ccomps * dcomps);

            auto g_xy_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 266 * ccomps * dcomps);

            auto g_xy_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 267 * ccomps * dcomps);

            auto g_xy_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 268 * ccomps * dcomps);

            auto g_xy_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 269 * ccomps * dcomps);

            auto g_xy_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 270 * ccomps * dcomps);

            auto g_xy_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 271 * ccomps * dcomps);

            auto g_xy_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 272 * ccomps * dcomps);

            auto g_xy_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 273 * ccomps * dcomps);

            auto g_xy_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 274 * ccomps * dcomps);

            auto g_xy_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 275 * ccomps * dcomps);

            auto g_xy_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 276 * ccomps * dcomps);

            auto g_xy_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 277 * ccomps * dcomps);

            auto g_xy_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 278 * ccomps * dcomps);

            auto g_xy_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 279 * ccomps * dcomps);

            auto g_xy_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 280 * ccomps * dcomps);

            auto g_xy_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 281 * ccomps * dcomps);

            auto g_xy_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 282 * ccomps * dcomps);

            auto g_xy_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 283 * ccomps * dcomps);

            auto g_xy_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 284 * ccomps * dcomps);

            auto g_xy_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 285 * ccomps * dcomps);

            auto g_xy_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 286 * ccomps * dcomps);

            auto g_xy_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 287 * ccomps * dcomps);

            auto g_xy_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 288 * ccomps * dcomps);

            auto g_xy_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 289 * ccomps * dcomps);

            auto g_xy_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 290 * ccomps * dcomps);

            auto g_xy_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 291 * ccomps * dcomps);

            auto g_xy_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 292 * ccomps * dcomps);

            auto g_xy_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 293 * ccomps * dcomps);

            auto g_xy_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 294 * ccomps * dcomps);

            auto g_xy_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 295 * ccomps * dcomps);

            auto g_xy_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 296 * ccomps * dcomps);

            auto g_xy_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 297 * ccomps * dcomps);

            auto g_xy_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 298 * ccomps * dcomps);

            auto g_xy_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 299 * ccomps * dcomps);

            auto g_xz_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 305 * ccomps * dcomps);

            auto g_xz_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 309 * ccomps * dcomps);

            auto g_xz_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 311 * ccomps * dcomps);

            auto g_xz_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 314 * ccomps * dcomps);

            auto g_xz_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 317 * ccomps * dcomps);

            auto g_xz_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 319 * ccomps * dcomps);

            auto g_xz_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 323 * ccomps * dcomps);

            auto g_xz_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 324 * ccomps * dcomps);

            auto g_xz_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 325 * ccomps * dcomps);

            auto g_xz_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 326 * ccomps * dcomps);

            auto g_xz_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 327 * ccomps * dcomps);

            auto g_xz_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 328 * ccomps * dcomps);

            auto g_xz_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 329 * ccomps * dcomps);

            auto g_xz_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 330 * ccomps * dcomps);

            auto g_xz_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 331 * ccomps * dcomps);

            auto g_xz_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 332 * ccomps * dcomps);

            auto g_xz_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 333 * ccomps * dcomps);

            auto g_xz_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 334 * ccomps * dcomps);

            auto g_xz_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 335 * ccomps * dcomps);

            auto g_xz_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 341 * ccomps * dcomps);

            auto g_xz_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 347 * ccomps * dcomps);

            auto g_xz_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 353 * ccomps * dcomps);

            auto g_xz_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 359 * ccomps * dcomps);

            auto g_xz_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 365 * ccomps * dcomps);

            auto g_xz_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 371 * ccomps * dcomps);

            auto g_xz_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 377 * ccomps * dcomps);

            auto g_xz_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 378 * ccomps * dcomps);

            auto g_xz_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 379 * ccomps * dcomps);

            auto g_xz_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 380 * ccomps * dcomps);

            auto g_xz_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 381 * ccomps * dcomps);

            auto g_xz_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 382 * ccomps * dcomps);

            auto g_xz_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 383 * ccomps * dcomps);

            auto g_xz_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 384 * ccomps * dcomps);

            auto g_xz_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 385 * ccomps * dcomps);

            auto g_xz_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 386 * ccomps * dcomps);

            auto g_xz_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 387 * ccomps * dcomps);

            auto g_xz_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 388 * ccomps * dcomps);

            auto g_xz_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 389 * ccomps * dcomps);

            auto g_xz_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 390 * ccomps * dcomps);

            auto g_xz_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 391 * ccomps * dcomps);

            auto g_xz_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 392 * ccomps * dcomps);

            auto g_xz_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 393 * ccomps * dcomps);

            auto g_xz_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 394 * ccomps * dcomps);

            auto g_xz_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 395 * ccomps * dcomps);

            auto g_xz_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 396 * ccomps * dcomps);

            auto g_xz_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 397 * ccomps * dcomps);

            auto g_xz_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 398 * ccomps * dcomps);

            auto g_xz_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 399 * ccomps * dcomps);

            auto g_xz_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 400 * ccomps * dcomps);

            auto g_xz_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 401 * ccomps * dcomps);

            auto g_xz_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 402 * ccomps * dcomps);

            auto g_xz_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 403 * ccomps * dcomps);

            auto g_xz_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 404 * ccomps * dcomps);

            auto g_xz_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 405 * ccomps * dcomps);

            auto g_xz_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 406 * ccomps * dcomps);

            auto g_xz_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 407 * ccomps * dcomps);

            auto g_xz_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 408 * ccomps * dcomps);

            auto g_xz_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 409 * ccomps * dcomps);

            auto g_xz_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 410 * ccomps * dcomps);

            auto g_xz_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 411 * ccomps * dcomps);

            auto g_xz_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 412 * ccomps * dcomps);

            auto g_xz_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 413 * ccomps * dcomps);

            auto g_xz_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 414 * ccomps * dcomps);

            auto g_xz_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 415 * ccomps * dcomps);

            auto g_xz_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 416 * ccomps * dcomps);

            auto g_xz_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 417 * ccomps * dcomps);

            auto g_xz_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 418 * ccomps * dcomps);

            auto g_xz_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 419 * ccomps * dcomps);

            auto g_xz_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 420 * ccomps * dcomps);

            auto g_xz_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 421 * ccomps * dcomps);

            auto g_xz_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 422 * ccomps * dcomps);

            auto g_xz_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 423 * ccomps * dcomps);

            auto g_xz_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 424 * ccomps * dcomps);

            auto g_xz_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 425 * ccomps * dcomps);

            auto g_xz_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 426 * ccomps * dcomps);

            auto g_xz_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 427 * ccomps * dcomps);

            auto g_xz_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 428 * ccomps * dcomps);

            auto g_xz_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 429 * ccomps * dcomps);

            auto g_xz_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 430 * ccomps * dcomps);

            auto g_xz_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 431 * ccomps * dcomps);

            auto g_xz_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 432 * ccomps * dcomps);

            auto g_xz_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 433 * ccomps * dcomps);

            auto g_xz_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 434 * ccomps * dcomps);

            auto g_xz_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 435 * ccomps * dcomps);

            auto g_xz_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 436 * ccomps * dcomps);

            auto g_xz_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 437 * ccomps * dcomps);

            auto g_xz_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 438 * ccomps * dcomps);

            auto g_xz_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 439 * ccomps * dcomps);

            auto g_xz_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 440 * ccomps * dcomps);

            auto g_xz_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 441 * ccomps * dcomps);

            auto g_xz_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 442 * ccomps * dcomps);

            auto g_xz_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 443 * ccomps * dcomps);

            auto g_xz_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 444 * ccomps * dcomps);

            auto g_xz_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 445 * ccomps * dcomps);

            auto g_xz_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 446 * ccomps * dcomps);

            auto g_xz_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 447 * ccomps * dcomps);

            auto g_xz_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 448 * ccomps * dcomps);

            auto g_xz_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 449 * ccomps * dcomps);

            auto g_yy_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 450 * ccomps * dcomps);

            auto g_yy_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 451 * ccomps * dcomps);

            auto g_yy_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 452 * ccomps * dcomps);

            auto g_yy_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 453 * ccomps * dcomps);

            auto g_yy_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 454 * ccomps * dcomps);

            auto g_yy_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 455 * ccomps * dcomps);

            auto g_yy_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 456 * ccomps * dcomps);

            auto g_yy_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 457 * ccomps * dcomps);

            auto g_yy_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 458 * ccomps * dcomps);

            auto g_yy_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 459 * ccomps * dcomps);

            auto g_yy_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 460 * ccomps * dcomps);

            auto g_yy_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 461 * ccomps * dcomps);

            auto g_yy_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 462 * ccomps * dcomps);

            auto g_yy_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 463 * ccomps * dcomps);

            auto g_yy_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 464 * ccomps * dcomps);

            auto g_yy_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 465 * ccomps * dcomps);

            auto g_yy_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 466 * ccomps * dcomps);

            auto g_yy_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 467 * ccomps * dcomps);

            auto g_yy_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 468 * ccomps * dcomps);

            auto g_yy_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 469 * ccomps * dcomps);

            auto g_yy_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 470 * ccomps * dcomps);

            auto g_yy_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 471 * ccomps * dcomps);

            auto g_yy_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 472 * ccomps * dcomps);

            auto g_yy_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 473 * ccomps * dcomps);

            auto g_yy_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 474 * ccomps * dcomps);

            auto g_yy_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 475 * ccomps * dcomps);

            auto g_yy_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 476 * ccomps * dcomps);

            auto g_yy_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 477 * ccomps * dcomps);

            auto g_yy_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 478 * ccomps * dcomps);

            auto g_yy_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 479 * ccomps * dcomps);

            auto g_yy_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 480 * ccomps * dcomps);

            auto g_yy_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 481 * ccomps * dcomps);

            auto g_yy_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 482 * ccomps * dcomps);

            auto g_yy_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 483 * ccomps * dcomps);

            auto g_yy_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 484 * ccomps * dcomps);

            auto g_yy_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 485 * ccomps * dcomps);

            auto g_yy_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 486 * ccomps * dcomps);

            auto g_yy_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 487 * ccomps * dcomps);

            auto g_yy_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 488 * ccomps * dcomps);

            auto g_yy_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 489 * ccomps * dcomps);

            auto g_yy_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 490 * ccomps * dcomps);

            auto g_yy_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 491 * ccomps * dcomps);

            auto g_yy_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 492 * ccomps * dcomps);

            auto g_yy_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 493 * ccomps * dcomps);

            auto g_yy_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 494 * ccomps * dcomps);

            auto g_yy_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 495 * ccomps * dcomps);

            auto g_yy_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 496 * ccomps * dcomps);

            auto g_yy_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 497 * ccomps * dcomps);

            auto g_yy_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 498 * ccomps * dcomps);

            auto g_yy_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 499 * ccomps * dcomps);

            auto g_yy_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 500 * ccomps * dcomps);

            auto g_yy_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 501 * ccomps * dcomps);

            auto g_yy_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 502 * ccomps * dcomps);

            auto g_yy_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 503 * ccomps * dcomps);

            auto g_yy_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 504 * ccomps * dcomps);

            auto g_yy_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 505 * ccomps * dcomps);

            auto g_yy_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 506 * ccomps * dcomps);

            auto g_yy_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 507 * ccomps * dcomps);

            auto g_yy_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 508 * ccomps * dcomps);

            auto g_yy_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 509 * ccomps * dcomps);

            auto g_yy_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 510 * ccomps * dcomps);

            auto g_yy_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 511 * ccomps * dcomps);

            auto g_yy_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 512 * ccomps * dcomps);

            auto g_yy_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 513 * ccomps * dcomps);

            auto g_yy_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 514 * ccomps * dcomps);

            auto g_yy_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 515 * ccomps * dcomps);

            auto g_yy_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 516 * ccomps * dcomps);

            auto g_yy_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 517 * ccomps * dcomps);

            auto g_yy_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 518 * ccomps * dcomps);

            auto g_yy_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 519 * ccomps * dcomps);

            auto g_yy_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 520 * ccomps * dcomps);

            auto g_yy_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 521 * ccomps * dcomps);

            auto g_yy_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 522 * ccomps * dcomps);

            auto g_yy_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 523 * ccomps * dcomps);

            auto g_yy_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 524 * ccomps * dcomps);

            auto g_yy_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 525 * ccomps * dcomps);

            auto g_yy_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 526 * ccomps * dcomps);

            auto g_yy_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 527 * ccomps * dcomps);

            auto g_yy_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 528 * ccomps * dcomps);

            auto g_yy_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 529 * ccomps * dcomps);

            auto g_yy_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 530 * ccomps * dcomps);

            auto g_yy_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 531 * ccomps * dcomps);

            auto g_yy_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 532 * ccomps * dcomps);

            auto g_yy_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 533 * ccomps * dcomps);

            auto g_yy_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 534 * ccomps * dcomps);

            auto g_yy_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 535 * ccomps * dcomps);

            auto g_yy_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 536 * ccomps * dcomps);

            auto g_yy_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 537 * ccomps * dcomps);

            auto g_yy_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 538 * ccomps * dcomps);

            auto g_yy_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 539 * ccomps * dcomps);

            auto g_yy_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 540 * ccomps * dcomps);

            auto g_yy_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 541 * ccomps * dcomps);

            auto g_yy_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 542 * ccomps * dcomps);

            auto g_yy_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 543 * ccomps * dcomps);

            auto g_yy_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 544 * ccomps * dcomps);

            auto g_yy_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 545 * ccomps * dcomps);

            auto g_yy_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 546 * ccomps * dcomps);

            auto g_yy_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 547 * ccomps * dcomps);

            auto g_yy_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 548 * ccomps * dcomps);

            auto g_yy_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 549 * ccomps * dcomps);

            auto g_yy_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 550 * ccomps * dcomps);

            auto g_yy_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 551 * ccomps * dcomps);

            auto g_yy_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 552 * ccomps * dcomps);

            auto g_yy_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 553 * ccomps * dcomps);

            auto g_yy_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 554 * ccomps * dcomps);

            auto g_yy_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 555 * ccomps * dcomps);

            auto g_yy_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 556 * ccomps * dcomps);

            auto g_yy_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 557 * ccomps * dcomps);

            auto g_yy_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 558 * ccomps * dcomps);

            auto g_yy_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 559 * ccomps * dcomps);

            auto g_yy_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 560 * ccomps * dcomps);

            auto g_yy_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 561 * ccomps * dcomps);

            auto g_yy_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 562 * ccomps * dcomps);

            auto g_yy_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 563 * ccomps * dcomps);

            auto g_yy_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 564 * ccomps * dcomps);

            auto g_yy_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 565 * ccomps * dcomps);

            auto g_yy_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 566 * ccomps * dcomps);

            auto g_yy_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 567 * ccomps * dcomps);

            auto g_yy_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 568 * ccomps * dcomps);

            auto g_yy_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 569 * ccomps * dcomps);

            auto g_yy_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 570 * ccomps * dcomps);

            auto g_yy_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 571 * ccomps * dcomps);

            auto g_yy_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 572 * ccomps * dcomps);

            auto g_yy_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 573 * ccomps * dcomps);

            auto g_yy_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 574 * ccomps * dcomps);

            auto g_yy_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 575 * ccomps * dcomps);

            auto g_yy_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 576 * ccomps * dcomps);

            auto g_yy_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 577 * ccomps * dcomps);

            auto g_yy_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 578 * ccomps * dcomps);

            auto g_yy_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 579 * ccomps * dcomps);

            auto g_yy_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 580 * ccomps * dcomps);

            auto g_yy_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 581 * ccomps * dcomps);

            auto g_yy_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 582 * ccomps * dcomps);

            auto g_yy_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 583 * ccomps * dcomps);

            auto g_yy_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 584 * ccomps * dcomps);

            auto g_yy_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 585 * ccomps * dcomps);

            auto g_yy_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 586 * ccomps * dcomps);

            auto g_yy_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 587 * ccomps * dcomps);

            auto g_yy_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 588 * ccomps * dcomps);

            auto g_yy_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 589 * ccomps * dcomps);

            auto g_yy_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 590 * ccomps * dcomps);

            auto g_yy_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 591 * ccomps * dcomps);

            auto g_yy_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 592 * ccomps * dcomps);

            auto g_yy_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 593 * ccomps * dcomps);

            auto g_yy_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 594 * ccomps * dcomps);

            auto g_yy_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 595 * ccomps * dcomps);

            auto g_yy_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 596 * ccomps * dcomps);

            auto g_yy_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 597 * ccomps * dcomps);

            auto g_yy_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 598 * ccomps * dcomps);

            auto g_yy_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 599 * ccomps * dcomps);

            auto g_yz_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 600 * ccomps * dcomps);

            auto g_yz_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 601 * ccomps * dcomps);

            auto g_yz_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 602 * ccomps * dcomps);

            auto g_yz_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 603 * ccomps * dcomps);

            auto g_yz_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 604 * ccomps * dcomps);

            auto g_yz_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 605 * ccomps * dcomps);

            auto g_yz_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 606 * ccomps * dcomps);

            auto g_yz_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 607 * ccomps * dcomps);

            auto g_yz_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 608 * ccomps * dcomps);

            auto g_yz_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 609 * ccomps * dcomps);

            auto g_yz_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 610 * ccomps * dcomps);

            auto g_yz_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 611 * ccomps * dcomps);

            auto g_yz_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 612 * ccomps * dcomps);

            auto g_yz_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 613 * ccomps * dcomps);

            auto g_yz_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 614 * ccomps * dcomps);

            auto g_yz_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 615 * ccomps * dcomps);

            auto g_yz_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 616 * ccomps * dcomps);

            auto g_yz_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 617 * ccomps * dcomps);

            auto g_yz_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 618 * ccomps * dcomps);

            auto g_yz_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 619 * ccomps * dcomps);

            auto g_yz_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 620 * ccomps * dcomps);

            auto g_yz_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 621 * ccomps * dcomps);

            auto g_yz_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 622 * ccomps * dcomps);

            auto g_yz_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 623 * ccomps * dcomps);

            auto g_yz_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 624 * ccomps * dcomps);

            auto g_yz_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 625 * ccomps * dcomps);

            auto g_yz_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 626 * ccomps * dcomps);

            auto g_yz_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 627 * ccomps * dcomps);

            auto g_yz_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 628 * ccomps * dcomps);

            auto g_yz_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 629 * ccomps * dcomps);

            auto g_yz_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 630 * ccomps * dcomps);

            auto g_yz_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 631 * ccomps * dcomps);

            auto g_yz_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 632 * ccomps * dcomps);

            auto g_yz_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 633 * ccomps * dcomps);

            auto g_yz_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 634 * ccomps * dcomps);

            auto g_yz_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 635 * ccomps * dcomps);

            auto g_yz_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 636 * ccomps * dcomps);

            auto g_yz_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 637 * ccomps * dcomps);

            auto g_yz_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 638 * ccomps * dcomps);

            auto g_yz_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 639 * ccomps * dcomps);

            auto g_yz_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 640 * ccomps * dcomps);

            auto g_yz_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 641 * ccomps * dcomps);

            auto g_yz_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 642 * ccomps * dcomps);

            auto g_yz_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 643 * ccomps * dcomps);

            auto g_yz_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 644 * ccomps * dcomps);

            auto g_yz_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 645 * ccomps * dcomps);

            auto g_yz_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 646 * ccomps * dcomps);

            auto g_yz_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 647 * ccomps * dcomps);

            auto g_yz_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 648 * ccomps * dcomps);

            auto g_yz_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 649 * ccomps * dcomps);

            auto g_yz_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 650 * ccomps * dcomps);

            auto g_yz_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 651 * ccomps * dcomps);

            auto g_yz_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 652 * ccomps * dcomps);

            auto g_yz_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 653 * ccomps * dcomps);

            auto g_yz_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 654 * ccomps * dcomps);

            auto g_yz_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 655 * ccomps * dcomps);

            auto g_yz_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 656 * ccomps * dcomps);

            auto g_yz_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 657 * ccomps * dcomps);

            auto g_yz_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 658 * ccomps * dcomps);

            auto g_yz_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 659 * ccomps * dcomps);

            auto g_yz_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 660 * ccomps * dcomps);

            auto g_yz_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 661 * ccomps * dcomps);

            auto g_yz_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 662 * ccomps * dcomps);

            auto g_yz_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 663 * ccomps * dcomps);

            auto g_yz_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 664 * ccomps * dcomps);

            auto g_yz_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 665 * ccomps * dcomps);

            auto g_yz_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 666 * ccomps * dcomps);

            auto g_yz_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 667 * ccomps * dcomps);

            auto g_yz_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 668 * ccomps * dcomps);

            auto g_yz_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 669 * ccomps * dcomps);

            auto g_yz_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 670 * ccomps * dcomps);

            auto g_yz_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 671 * ccomps * dcomps);

            auto g_yz_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 672 * ccomps * dcomps);

            auto g_yz_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 673 * ccomps * dcomps);

            auto g_yz_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 674 * ccomps * dcomps);

            auto g_yz_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 675 * ccomps * dcomps);

            auto g_yz_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 676 * ccomps * dcomps);

            auto g_yz_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 677 * ccomps * dcomps);

            auto g_yz_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 678 * ccomps * dcomps);

            auto g_yz_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 679 * ccomps * dcomps);

            auto g_yz_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 680 * ccomps * dcomps);

            auto g_yz_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 681 * ccomps * dcomps);

            auto g_yz_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 682 * ccomps * dcomps);

            auto g_yz_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 683 * ccomps * dcomps);

            auto g_yz_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 684 * ccomps * dcomps);

            auto g_yz_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 685 * ccomps * dcomps);

            auto g_yz_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 686 * ccomps * dcomps);

            auto g_yz_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 687 * ccomps * dcomps);

            auto g_yz_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 688 * ccomps * dcomps);

            auto g_yz_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 689 * ccomps * dcomps);

            auto g_yz_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 690 * ccomps * dcomps);

            auto g_yz_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 691 * ccomps * dcomps);

            auto g_yz_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 692 * ccomps * dcomps);

            auto g_yz_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 693 * ccomps * dcomps);

            auto g_yz_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 694 * ccomps * dcomps);

            auto g_yz_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 695 * ccomps * dcomps);

            auto g_yz_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 696 * ccomps * dcomps);

            auto g_yz_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 697 * ccomps * dcomps);

            auto g_yz_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 698 * ccomps * dcomps);

            auto g_yz_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 699 * ccomps * dcomps);

            auto g_yz_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 700 * ccomps * dcomps);

            auto g_yz_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 701 * ccomps * dcomps);

            auto g_yz_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 702 * ccomps * dcomps);

            auto g_yz_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 703 * ccomps * dcomps);

            auto g_yz_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 704 * ccomps * dcomps);

            auto g_yz_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 705 * ccomps * dcomps);

            auto g_yz_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 706 * ccomps * dcomps);

            auto g_yz_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 707 * ccomps * dcomps);

            auto g_yz_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 708 * ccomps * dcomps);

            auto g_yz_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 709 * ccomps * dcomps);

            auto g_yz_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 710 * ccomps * dcomps);

            auto g_yz_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 711 * ccomps * dcomps);

            auto g_yz_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 712 * ccomps * dcomps);

            auto g_yz_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 713 * ccomps * dcomps);

            auto g_yz_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 714 * ccomps * dcomps);

            auto g_yz_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 715 * ccomps * dcomps);

            auto g_yz_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 716 * ccomps * dcomps);

            auto g_yz_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 717 * ccomps * dcomps);

            auto g_yz_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 718 * ccomps * dcomps);

            auto g_yz_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 719 * ccomps * dcomps);

            auto g_yz_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 720 * ccomps * dcomps);

            auto g_yz_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 721 * ccomps * dcomps);

            auto g_yz_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 722 * ccomps * dcomps);

            auto g_yz_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 723 * ccomps * dcomps);

            auto g_yz_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 724 * ccomps * dcomps);

            auto g_yz_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 725 * ccomps * dcomps);

            auto g_yz_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 726 * ccomps * dcomps);

            auto g_yz_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 727 * ccomps * dcomps);

            auto g_yz_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 728 * ccomps * dcomps);

            auto g_yz_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 729 * ccomps * dcomps);

            auto g_yz_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 730 * ccomps * dcomps);

            auto g_yz_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 731 * ccomps * dcomps);

            auto g_yz_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 732 * ccomps * dcomps);

            auto g_yz_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 733 * ccomps * dcomps);

            auto g_yz_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 734 * ccomps * dcomps);

            auto g_yz_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 735 * ccomps * dcomps);

            auto g_yz_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 736 * ccomps * dcomps);

            auto g_yz_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 737 * ccomps * dcomps);

            auto g_yz_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 738 * ccomps * dcomps);

            auto g_yz_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 739 * ccomps * dcomps);

            auto g_yz_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 740 * ccomps * dcomps);

            auto g_yz_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 741 * ccomps * dcomps);

            auto g_yz_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 742 * ccomps * dcomps);

            auto g_yz_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 743 * ccomps * dcomps);

            auto g_yz_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 744 * ccomps * dcomps);

            auto g_yz_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 745 * ccomps * dcomps);

            auto g_yz_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 746 * ccomps * dcomps);

            auto g_yz_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 747 * ccomps * dcomps);

            auto g_yz_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 748 * ccomps * dcomps);

            auto g_yz_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 749 * ccomps * dcomps);

            auto g_zz_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 750 * ccomps * dcomps);

            auto g_zz_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 751 * ccomps * dcomps);

            auto g_zz_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 752 * ccomps * dcomps);

            auto g_zz_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 753 * ccomps * dcomps);

            auto g_zz_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 754 * ccomps * dcomps);

            auto g_zz_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 755 * ccomps * dcomps);

            auto g_zz_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 756 * ccomps * dcomps);

            auto g_zz_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 757 * ccomps * dcomps);

            auto g_zz_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 758 * ccomps * dcomps);

            auto g_zz_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 759 * ccomps * dcomps);

            auto g_zz_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 760 * ccomps * dcomps);

            auto g_zz_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 761 * ccomps * dcomps);

            auto g_zz_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 762 * ccomps * dcomps);

            auto g_zz_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 763 * ccomps * dcomps);

            auto g_zz_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 764 * ccomps * dcomps);

            auto g_zz_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 765 * ccomps * dcomps);

            auto g_zz_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 766 * ccomps * dcomps);

            auto g_zz_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 767 * ccomps * dcomps);

            auto g_zz_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 768 * ccomps * dcomps);

            auto g_zz_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 769 * ccomps * dcomps);

            auto g_zz_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 770 * ccomps * dcomps);

            auto g_zz_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 771 * ccomps * dcomps);

            auto g_zz_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 772 * ccomps * dcomps);

            auto g_zz_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 773 * ccomps * dcomps);

            auto g_zz_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 774 * ccomps * dcomps);

            auto g_zz_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 775 * ccomps * dcomps);

            auto g_zz_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 776 * ccomps * dcomps);

            auto g_zz_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 777 * ccomps * dcomps);

            auto g_zz_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 778 * ccomps * dcomps);

            auto g_zz_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 779 * ccomps * dcomps);

            auto g_zz_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 780 * ccomps * dcomps);

            auto g_zz_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 781 * ccomps * dcomps);

            auto g_zz_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 782 * ccomps * dcomps);

            auto g_zz_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 783 * ccomps * dcomps);

            auto g_zz_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 784 * ccomps * dcomps);

            auto g_zz_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 785 * ccomps * dcomps);

            auto g_zz_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 786 * ccomps * dcomps);

            auto g_zz_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 787 * ccomps * dcomps);

            auto g_zz_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 788 * ccomps * dcomps);

            auto g_zz_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 789 * ccomps * dcomps);

            auto g_zz_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 790 * ccomps * dcomps);

            auto g_zz_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 791 * ccomps * dcomps);

            auto g_zz_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 792 * ccomps * dcomps);

            auto g_zz_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 793 * ccomps * dcomps);

            auto g_zz_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 794 * ccomps * dcomps);

            auto g_zz_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 795 * ccomps * dcomps);

            auto g_zz_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 796 * ccomps * dcomps);

            auto g_zz_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 797 * ccomps * dcomps);

            auto g_zz_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 798 * ccomps * dcomps);

            auto g_zz_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 799 * ccomps * dcomps);

            auto g_zz_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 800 * ccomps * dcomps);

            auto g_zz_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 801 * ccomps * dcomps);

            auto g_zz_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 802 * ccomps * dcomps);

            auto g_zz_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 803 * ccomps * dcomps);

            auto g_zz_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 804 * ccomps * dcomps);

            auto g_zz_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 805 * ccomps * dcomps);

            auto g_zz_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 806 * ccomps * dcomps);

            auto g_zz_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 807 * ccomps * dcomps);

            auto g_zz_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 808 * ccomps * dcomps);

            auto g_zz_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 809 * ccomps * dcomps);

            auto g_zz_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 810 * ccomps * dcomps);

            auto g_zz_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 811 * ccomps * dcomps);

            auto g_zz_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 812 * ccomps * dcomps);

            auto g_zz_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 813 * ccomps * dcomps);

            auto g_zz_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 814 * ccomps * dcomps);

            auto g_zz_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 815 * ccomps * dcomps);

            auto g_zz_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 816 * ccomps * dcomps);

            auto g_zz_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 817 * ccomps * dcomps);

            auto g_zz_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 818 * ccomps * dcomps);

            auto g_zz_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 819 * ccomps * dcomps);

            auto g_zz_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 820 * ccomps * dcomps);

            auto g_zz_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 821 * ccomps * dcomps);

            auto g_zz_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 822 * ccomps * dcomps);

            auto g_zz_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 823 * ccomps * dcomps);

            auto g_zz_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 824 * ccomps * dcomps);

            auto g_zz_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 825 * ccomps * dcomps);

            auto g_zz_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 826 * ccomps * dcomps);

            auto g_zz_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 827 * ccomps * dcomps);

            auto g_zz_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 828 * ccomps * dcomps);

            auto g_zz_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 829 * ccomps * dcomps);

            auto g_zz_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 830 * ccomps * dcomps);

            auto g_zz_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 831 * ccomps * dcomps);

            auto g_zz_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 832 * ccomps * dcomps);

            auto g_zz_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 833 * ccomps * dcomps);

            auto g_zz_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 834 * ccomps * dcomps);

            auto g_zz_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 835 * ccomps * dcomps);

            auto g_zz_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 836 * ccomps * dcomps);

            auto g_zz_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 837 * ccomps * dcomps);

            auto g_zz_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 838 * ccomps * dcomps);

            auto g_zz_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 839 * ccomps * dcomps);

            auto g_zz_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 840 * ccomps * dcomps);

            auto g_zz_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 841 * ccomps * dcomps);

            auto g_zz_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 842 * ccomps * dcomps);

            auto g_zz_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 843 * ccomps * dcomps);

            auto g_zz_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 844 * ccomps * dcomps);

            auto g_zz_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 845 * ccomps * dcomps);

            auto g_zz_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 846 * ccomps * dcomps);

            auto g_zz_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 847 * ccomps * dcomps);

            auto g_zz_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 848 * ccomps * dcomps);

            auto g_zz_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 849 * ccomps * dcomps);

            auto g_zz_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 850 * ccomps * dcomps);

            auto g_zz_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 851 * ccomps * dcomps);

            auto g_zz_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 852 * ccomps * dcomps);

            auto g_zz_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 853 * ccomps * dcomps);

            auto g_zz_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 854 * ccomps * dcomps);

            auto g_zz_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 855 * ccomps * dcomps);

            auto g_zz_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 856 * ccomps * dcomps);

            auto g_zz_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 857 * ccomps * dcomps);

            auto g_zz_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 858 * ccomps * dcomps);

            auto g_zz_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 859 * ccomps * dcomps);

            auto g_zz_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 860 * ccomps * dcomps);

            auto g_zz_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 861 * ccomps * dcomps);

            auto g_zz_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 862 * ccomps * dcomps);

            auto g_zz_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 863 * ccomps * dcomps);

            auto g_zz_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 864 * ccomps * dcomps);

            auto g_zz_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 865 * ccomps * dcomps);

            auto g_zz_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 866 * ccomps * dcomps);

            auto g_zz_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 867 * ccomps * dcomps);

            auto g_zz_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 868 * ccomps * dcomps);

            auto g_zz_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 869 * ccomps * dcomps);

            auto g_zz_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 870 * ccomps * dcomps);

            auto g_zz_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 871 * ccomps * dcomps);

            auto g_zz_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 872 * ccomps * dcomps);

            auto g_zz_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 873 * ccomps * dcomps);

            auto g_zz_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 874 * ccomps * dcomps);

            auto g_zz_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 875 * ccomps * dcomps);

            auto g_zz_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 876 * ccomps * dcomps);

            auto g_zz_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 877 * ccomps * dcomps);

            auto g_zz_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 878 * ccomps * dcomps);

            auto g_zz_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 879 * ccomps * dcomps);

            auto g_zz_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 880 * ccomps * dcomps);

            auto g_zz_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 881 * ccomps * dcomps);

            auto g_zz_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 882 * ccomps * dcomps);

            auto g_zz_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 883 * ccomps * dcomps);

            auto g_zz_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 884 * ccomps * dcomps);

            auto g_zz_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 885 * ccomps * dcomps);

            auto g_zz_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 886 * ccomps * dcomps);

            auto g_zz_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 887 * ccomps * dcomps);

            auto g_zz_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 888 * ccomps * dcomps);

            auto g_zz_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 889 * ccomps * dcomps);

            auto g_zz_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 890 * ccomps * dcomps);

            auto g_zz_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 891 * ccomps * dcomps);

            auto g_zz_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 892 * ccomps * dcomps);

            auto g_zz_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 893 * ccomps * dcomps);

            auto g_zz_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 894 * ccomps * dcomps);

            auto g_zz_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 895 * ccomps * dcomps);

            auto g_zz_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 896 * ccomps * dcomps);

            auto g_zz_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 897 * ccomps * dcomps);

            auto g_zz_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 898 * ccomps * dcomps);

            auto g_zz_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 899 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gfxx

            const auto gf_geom_20_off = idx_geom_20_gfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_xxx, g_x_0_xxx_xxy, g_x_0_xxx_xxz, g_x_0_xxx_xyy, g_x_0_xxx_xyz, g_x_0_xxx_xzz, g_x_0_xxx_yyy, g_x_0_xxx_yyz, g_x_0_xxx_yzz, g_x_0_xxx_zzz, g_xx_0_xxx_xxx, g_xx_0_xxx_xxxx, g_xx_0_xxx_xxxy, g_xx_0_xxx_xxxz, g_xx_0_xxx_xxy, g_xx_0_xxx_xxyy, g_xx_0_xxx_xxyz, g_xx_0_xxx_xxz, g_xx_0_xxx_xxzz, g_xx_0_xxx_xyy, g_xx_0_xxx_xyyy, g_xx_0_xxx_xyyz, g_xx_0_xxx_xyz, g_xx_0_xxx_xyzz, g_xx_0_xxx_xzz, g_xx_0_xxx_xzzz, g_xx_0_xxx_yyy, g_xx_0_xxx_yyz, g_xx_0_xxx_yzz, g_xx_0_xxx_zzz, g_xx_0_xxxx_xxx, g_xx_0_xxxx_xxy, g_xx_0_xxxx_xxz, g_xx_0_xxxx_xyy, g_xx_0_xxxx_xyz, g_xx_0_xxxx_xzz, g_xx_0_xxxx_yyy, g_xx_0_xxxx_yyz, g_xx_0_xxxx_yzz, g_xx_0_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxx_xxx[k] = -2.0 * g_x_0_xxx_xxx[k] - g_xx_0_xxx_xxx[k] * ab_x + g_xx_0_xxx_xxxx[k];

                g_xx_0_xxxx_xxy[k] = -2.0 * g_x_0_xxx_xxy[k] - g_xx_0_xxx_xxy[k] * ab_x + g_xx_0_xxx_xxxy[k];

                g_xx_0_xxxx_xxz[k] = -2.0 * g_x_0_xxx_xxz[k] - g_xx_0_xxx_xxz[k] * ab_x + g_xx_0_xxx_xxxz[k];

                g_xx_0_xxxx_xyy[k] = -2.0 * g_x_0_xxx_xyy[k] - g_xx_0_xxx_xyy[k] * ab_x + g_xx_0_xxx_xxyy[k];

                g_xx_0_xxxx_xyz[k] = -2.0 * g_x_0_xxx_xyz[k] - g_xx_0_xxx_xyz[k] * ab_x + g_xx_0_xxx_xxyz[k];

                g_xx_0_xxxx_xzz[k] = -2.0 * g_x_0_xxx_xzz[k] - g_xx_0_xxx_xzz[k] * ab_x + g_xx_0_xxx_xxzz[k];

                g_xx_0_xxxx_yyy[k] = -2.0 * g_x_0_xxx_yyy[k] - g_xx_0_xxx_yyy[k] * ab_x + g_xx_0_xxx_xyyy[k];

                g_xx_0_xxxx_yyz[k] = -2.0 * g_x_0_xxx_yyz[k] - g_xx_0_xxx_yyz[k] * ab_x + g_xx_0_xxx_xyyz[k];

                g_xx_0_xxxx_yzz[k] = -2.0 * g_x_0_xxx_yzz[k] - g_xx_0_xxx_yzz[k] * ab_x + g_xx_0_xxx_xyzz[k];

                g_xx_0_xxxx_zzz[k] = -2.0 * g_x_0_xxx_zzz[k] - g_xx_0_xxx_zzz[k] * ab_x + g_xx_0_xxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxx_xxx, g_xx_0_xxx_xxxy, g_xx_0_xxx_xxy, g_xx_0_xxx_xxyy, g_xx_0_xxx_xxyz, g_xx_0_xxx_xxz, g_xx_0_xxx_xyy, g_xx_0_xxx_xyyy, g_xx_0_xxx_xyyz, g_xx_0_xxx_xyz, g_xx_0_xxx_xyzz, g_xx_0_xxx_xzz, g_xx_0_xxx_yyy, g_xx_0_xxx_yyyy, g_xx_0_xxx_yyyz, g_xx_0_xxx_yyz, g_xx_0_xxx_yyzz, g_xx_0_xxx_yzz, g_xx_0_xxx_yzzz, g_xx_0_xxx_zzz, g_xx_0_xxxy_xxx, g_xx_0_xxxy_xxy, g_xx_0_xxxy_xxz, g_xx_0_xxxy_xyy, g_xx_0_xxxy_xyz, g_xx_0_xxxy_xzz, g_xx_0_xxxy_yyy, g_xx_0_xxxy_yyz, g_xx_0_xxxy_yzz, g_xx_0_xxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxy_xxx[k] = -g_xx_0_xxx_xxx[k] * ab_y + g_xx_0_xxx_xxxy[k];

                g_xx_0_xxxy_xxy[k] = -g_xx_0_xxx_xxy[k] * ab_y + g_xx_0_xxx_xxyy[k];

                g_xx_0_xxxy_xxz[k] = -g_xx_0_xxx_xxz[k] * ab_y + g_xx_0_xxx_xxyz[k];

                g_xx_0_xxxy_xyy[k] = -g_xx_0_xxx_xyy[k] * ab_y + g_xx_0_xxx_xyyy[k];

                g_xx_0_xxxy_xyz[k] = -g_xx_0_xxx_xyz[k] * ab_y + g_xx_0_xxx_xyyz[k];

                g_xx_0_xxxy_xzz[k] = -g_xx_0_xxx_xzz[k] * ab_y + g_xx_0_xxx_xyzz[k];

                g_xx_0_xxxy_yyy[k] = -g_xx_0_xxx_yyy[k] * ab_y + g_xx_0_xxx_yyyy[k];

                g_xx_0_xxxy_yyz[k] = -g_xx_0_xxx_yyz[k] * ab_y + g_xx_0_xxx_yyyz[k];

                g_xx_0_xxxy_yzz[k] = -g_xx_0_xxx_yzz[k] * ab_y + g_xx_0_xxx_yyzz[k];

                g_xx_0_xxxy_zzz[k] = -g_xx_0_xxx_zzz[k] * ab_y + g_xx_0_xxx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxx_xxx, g_xx_0_xxx_xxxz, g_xx_0_xxx_xxy, g_xx_0_xxx_xxyz, g_xx_0_xxx_xxz, g_xx_0_xxx_xxzz, g_xx_0_xxx_xyy, g_xx_0_xxx_xyyz, g_xx_0_xxx_xyz, g_xx_0_xxx_xyzz, g_xx_0_xxx_xzz, g_xx_0_xxx_xzzz, g_xx_0_xxx_yyy, g_xx_0_xxx_yyyz, g_xx_0_xxx_yyz, g_xx_0_xxx_yyzz, g_xx_0_xxx_yzz, g_xx_0_xxx_yzzz, g_xx_0_xxx_zzz, g_xx_0_xxx_zzzz, g_xx_0_xxxz_xxx, g_xx_0_xxxz_xxy, g_xx_0_xxxz_xxz, g_xx_0_xxxz_xyy, g_xx_0_xxxz_xyz, g_xx_0_xxxz_xzz, g_xx_0_xxxz_yyy, g_xx_0_xxxz_yyz, g_xx_0_xxxz_yzz, g_xx_0_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxz_xxx[k] = -g_xx_0_xxx_xxx[k] * ab_z + g_xx_0_xxx_xxxz[k];

                g_xx_0_xxxz_xxy[k] = -g_xx_0_xxx_xxy[k] * ab_z + g_xx_0_xxx_xxyz[k];

                g_xx_0_xxxz_xxz[k] = -g_xx_0_xxx_xxz[k] * ab_z + g_xx_0_xxx_xxzz[k];

                g_xx_0_xxxz_xyy[k] = -g_xx_0_xxx_xyy[k] * ab_z + g_xx_0_xxx_xyyz[k];

                g_xx_0_xxxz_xyz[k] = -g_xx_0_xxx_xyz[k] * ab_z + g_xx_0_xxx_xyzz[k];

                g_xx_0_xxxz_xzz[k] = -g_xx_0_xxx_xzz[k] * ab_z + g_xx_0_xxx_xzzz[k];

                g_xx_0_xxxz_yyy[k] = -g_xx_0_xxx_yyy[k] * ab_z + g_xx_0_xxx_yyyz[k];

                g_xx_0_xxxz_yyz[k] = -g_xx_0_xxx_yyz[k] * ab_z + g_xx_0_xxx_yyzz[k];

                g_xx_0_xxxz_yzz[k] = -g_xx_0_xxx_yzz[k] * ab_z + g_xx_0_xxx_yzzz[k];

                g_xx_0_xxxz_zzz[k] = -g_xx_0_xxx_zzz[k] * ab_z + g_xx_0_xxx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxy_xxx, g_xx_0_xxy_xxxy, g_xx_0_xxy_xxy, g_xx_0_xxy_xxyy, g_xx_0_xxy_xxyz, g_xx_0_xxy_xxz, g_xx_0_xxy_xyy, g_xx_0_xxy_xyyy, g_xx_0_xxy_xyyz, g_xx_0_xxy_xyz, g_xx_0_xxy_xyzz, g_xx_0_xxy_xzz, g_xx_0_xxy_yyy, g_xx_0_xxy_yyyy, g_xx_0_xxy_yyyz, g_xx_0_xxy_yyz, g_xx_0_xxy_yyzz, g_xx_0_xxy_yzz, g_xx_0_xxy_yzzz, g_xx_0_xxy_zzz, g_xx_0_xxyy_xxx, g_xx_0_xxyy_xxy, g_xx_0_xxyy_xxz, g_xx_0_xxyy_xyy, g_xx_0_xxyy_xyz, g_xx_0_xxyy_xzz, g_xx_0_xxyy_yyy, g_xx_0_xxyy_yyz, g_xx_0_xxyy_yzz, g_xx_0_xxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyy_xxx[k] = -g_xx_0_xxy_xxx[k] * ab_y + g_xx_0_xxy_xxxy[k];

                g_xx_0_xxyy_xxy[k] = -g_xx_0_xxy_xxy[k] * ab_y + g_xx_0_xxy_xxyy[k];

                g_xx_0_xxyy_xxz[k] = -g_xx_0_xxy_xxz[k] * ab_y + g_xx_0_xxy_xxyz[k];

                g_xx_0_xxyy_xyy[k] = -g_xx_0_xxy_xyy[k] * ab_y + g_xx_0_xxy_xyyy[k];

                g_xx_0_xxyy_xyz[k] = -g_xx_0_xxy_xyz[k] * ab_y + g_xx_0_xxy_xyyz[k];

                g_xx_0_xxyy_xzz[k] = -g_xx_0_xxy_xzz[k] * ab_y + g_xx_0_xxy_xyzz[k];

                g_xx_0_xxyy_yyy[k] = -g_xx_0_xxy_yyy[k] * ab_y + g_xx_0_xxy_yyyy[k];

                g_xx_0_xxyy_yyz[k] = -g_xx_0_xxy_yyz[k] * ab_y + g_xx_0_xxy_yyyz[k];

                g_xx_0_xxyy_yzz[k] = -g_xx_0_xxy_yzz[k] * ab_y + g_xx_0_xxy_yyzz[k];

                g_xx_0_xxyy_zzz[k] = -g_xx_0_xxy_zzz[k] * ab_y + g_xx_0_xxy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyz_xxx, g_xx_0_xxyz_xxy, g_xx_0_xxyz_xxz, g_xx_0_xxyz_xyy, g_xx_0_xxyz_xyz, g_xx_0_xxyz_xzz, g_xx_0_xxyz_yyy, g_xx_0_xxyz_yyz, g_xx_0_xxyz_yzz, g_xx_0_xxyz_zzz, g_xx_0_xxz_xxx, g_xx_0_xxz_xxxy, g_xx_0_xxz_xxy, g_xx_0_xxz_xxyy, g_xx_0_xxz_xxyz, g_xx_0_xxz_xxz, g_xx_0_xxz_xyy, g_xx_0_xxz_xyyy, g_xx_0_xxz_xyyz, g_xx_0_xxz_xyz, g_xx_0_xxz_xyzz, g_xx_0_xxz_xzz, g_xx_0_xxz_yyy, g_xx_0_xxz_yyyy, g_xx_0_xxz_yyyz, g_xx_0_xxz_yyz, g_xx_0_xxz_yyzz, g_xx_0_xxz_yzz, g_xx_0_xxz_yzzz, g_xx_0_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyz_xxx[k] = -g_xx_0_xxz_xxx[k] * ab_y + g_xx_0_xxz_xxxy[k];

                g_xx_0_xxyz_xxy[k] = -g_xx_0_xxz_xxy[k] * ab_y + g_xx_0_xxz_xxyy[k];

                g_xx_0_xxyz_xxz[k] = -g_xx_0_xxz_xxz[k] * ab_y + g_xx_0_xxz_xxyz[k];

                g_xx_0_xxyz_xyy[k] = -g_xx_0_xxz_xyy[k] * ab_y + g_xx_0_xxz_xyyy[k];

                g_xx_0_xxyz_xyz[k] = -g_xx_0_xxz_xyz[k] * ab_y + g_xx_0_xxz_xyyz[k];

                g_xx_0_xxyz_xzz[k] = -g_xx_0_xxz_xzz[k] * ab_y + g_xx_0_xxz_xyzz[k];

                g_xx_0_xxyz_yyy[k] = -g_xx_0_xxz_yyy[k] * ab_y + g_xx_0_xxz_yyyy[k];

                g_xx_0_xxyz_yyz[k] = -g_xx_0_xxz_yyz[k] * ab_y + g_xx_0_xxz_yyyz[k];

                g_xx_0_xxyz_yzz[k] = -g_xx_0_xxz_yzz[k] * ab_y + g_xx_0_xxz_yyzz[k];

                g_xx_0_xxyz_zzz[k] = -g_xx_0_xxz_zzz[k] * ab_y + g_xx_0_xxz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxz_xxx, g_xx_0_xxz_xxxz, g_xx_0_xxz_xxy, g_xx_0_xxz_xxyz, g_xx_0_xxz_xxz, g_xx_0_xxz_xxzz, g_xx_0_xxz_xyy, g_xx_0_xxz_xyyz, g_xx_0_xxz_xyz, g_xx_0_xxz_xyzz, g_xx_0_xxz_xzz, g_xx_0_xxz_xzzz, g_xx_0_xxz_yyy, g_xx_0_xxz_yyyz, g_xx_0_xxz_yyz, g_xx_0_xxz_yyzz, g_xx_0_xxz_yzz, g_xx_0_xxz_yzzz, g_xx_0_xxz_zzz, g_xx_0_xxz_zzzz, g_xx_0_xxzz_xxx, g_xx_0_xxzz_xxy, g_xx_0_xxzz_xxz, g_xx_0_xxzz_xyy, g_xx_0_xxzz_xyz, g_xx_0_xxzz_xzz, g_xx_0_xxzz_yyy, g_xx_0_xxzz_yyz, g_xx_0_xxzz_yzz, g_xx_0_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxzz_xxx[k] = -g_xx_0_xxz_xxx[k] * ab_z + g_xx_0_xxz_xxxz[k];

                g_xx_0_xxzz_xxy[k] = -g_xx_0_xxz_xxy[k] * ab_z + g_xx_0_xxz_xxyz[k];

                g_xx_0_xxzz_xxz[k] = -g_xx_0_xxz_xxz[k] * ab_z + g_xx_0_xxz_xxzz[k];

                g_xx_0_xxzz_xyy[k] = -g_xx_0_xxz_xyy[k] * ab_z + g_xx_0_xxz_xyyz[k];

                g_xx_0_xxzz_xyz[k] = -g_xx_0_xxz_xyz[k] * ab_z + g_xx_0_xxz_xyzz[k];

                g_xx_0_xxzz_xzz[k] = -g_xx_0_xxz_xzz[k] * ab_z + g_xx_0_xxz_xzzz[k];

                g_xx_0_xxzz_yyy[k] = -g_xx_0_xxz_yyy[k] * ab_z + g_xx_0_xxz_yyyz[k];

                g_xx_0_xxzz_yyz[k] = -g_xx_0_xxz_yyz[k] * ab_z + g_xx_0_xxz_yyzz[k];

                g_xx_0_xxzz_yzz[k] = -g_xx_0_xxz_yzz[k] * ab_z + g_xx_0_xxz_yzzz[k];

                g_xx_0_xxzz_zzz[k] = -g_xx_0_xxz_zzz[k] * ab_z + g_xx_0_xxz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyy_xxx, g_xx_0_xyy_xxxy, g_xx_0_xyy_xxy, g_xx_0_xyy_xxyy, g_xx_0_xyy_xxyz, g_xx_0_xyy_xxz, g_xx_0_xyy_xyy, g_xx_0_xyy_xyyy, g_xx_0_xyy_xyyz, g_xx_0_xyy_xyz, g_xx_0_xyy_xyzz, g_xx_0_xyy_xzz, g_xx_0_xyy_yyy, g_xx_0_xyy_yyyy, g_xx_0_xyy_yyyz, g_xx_0_xyy_yyz, g_xx_0_xyy_yyzz, g_xx_0_xyy_yzz, g_xx_0_xyy_yzzz, g_xx_0_xyy_zzz, g_xx_0_xyyy_xxx, g_xx_0_xyyy_xxy, g_xx_0_xyyy_xxz, g_xx_0_xyyy_xyy, g_xx_0_xyyy_xyz, g_xx_0_xyyy_xzz, g_xx_0_xyyy_yyy, g_xx_0_xyyy_yyz, g_xx_0_xyyy_yzz, g_xx_0_xyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyy_xxx[k] = -g_xx_0_xyy_xxx[k] * ab_y + g_xx_0_xyy_xxxy[k];

                g_xx_0_xyyy_xxy[k] = -g_xx_0_xyy_xxy[k] * ab_y + g_xx_0_xyy_xxyy[k];

                g_xx_0_xyyy_xxz[k] = -g_xx_0_xyy_xxz[k] * ab_y + g_xx_0_xyy_xxyz[k];

                g_xx_0_xyyy_xyy[k] = -g_xx_0_xyy_xyy[k] * ab_y + g_xx_0_xyy_xyyy[k];

                g_xx_0_xyyy_xyz[k] = -g_xx_0_xyy_xyz[k] * ab_y + g_xx_0_xyy_xyyz[k];

                g_xx_0_xyyy_xzz[k] = -g_xx_0_xyy_xzz[k] * ab_y + g_xx_0_xyy_xyzz[k];

                g_xx_0_xyyy_yyy[k] = -g_xx_0_xyy_yyy[k] * ab_y + g_xx_0_xyy_yyyy[k];

                g_xx_0_xyyy_yyz[k] = -g_xx_0_xyy_yyz[k] * ab_y + g_xx_0_xyy_yyyz[k];

                g_xx_0_xyyy_yzz[k] = -g_xx_0_xyy_yzz[k] * ab_y + g_xx_0_xyy_yyzz[k];

                g_xx_0_xyyy_zzz[k] = -g_xx_0_xyy_zzz[k] * ab_y + g_xx_0_xyy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyz_xxx, g_xx_0_xyyz_xxy, g_xx_0_xyyz_xxz, g_xx_0_xyyz_xyy, g_xx_0_xyyz_xyz, g_xx_0_xyyz_xzz, g_xx_0_xyyz_yyy, g_xx_0_xyyz_yyz, g_xx_0_xyyz_yzz, g_xx_0_xyyz_zzz, g_xx_0_xyz_xxx, g_xx_0_xyz_xxxy, g_xx_0_xyz_xxy, g_xx_0_xyz_xxyy, g_xx_0_xyz_xxyz, g_xx_0_xyz_xxz, g_xx_0_xyz_xyy, g_xx_0_xyz_xyyy, g_xx_0_xyz_xyyz, g_xx_0_xyz_xyz, g_xx_0_xyz_xyzz, g_xx_0_xyz_xzz, g_xx_0_xyz_yyy, g_xx_0_xyz_yyyy, g_xx_0_xyz_yyyz, g_xx_0_xyz_yyz, g_xx_0_xyz_yyzz, g_xx_0_xyz_yzz, g_xx_0_xyz_yzzz, g_xx_0_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyz_xxx[k] = -g_xx_0_xyz_xxx[k] * ab_y + g_xx_0_xyz_xxxy[k];

                g_xx_0_xyyz_xxy[k] = -g_xx_0_xyz_xxy[k] * ab_y + g_xx_0_xyz_xxyy[k];

                g_xx_0_xyyz_xxz[k] = -g_xx_0_xyz_xxz[k] * ab_y + g_xx_0_xyz_xxyz[k];

                g_xx_0_xyyz_xyy[k] = -g_xx_0_xyz_xyy[k] * ab_y + g_xx_0_xyz_xyyy[k];

                g_xx_0_xyyz_xyz[k] = -g_xx_0_xyz_xyz[k] * ab_y + g_xx_0_xyz_xyyz[k];

                g_xx_0_xyyz_xzz[k] = -g_xx_0_xyz_xzz[k] * ab_y + g_xx_0_xyz_xyzz[k];

                g_xx_0_xyyz_yyy[k] = -g_xx_0_xyz_yyy[k] * ab_y + g_xx_0_xyz_yyyy[k];

                g_xx_0_xyyz_yyz[k] = -g_xx_0_xyz_yyz[k] * ab_y + g_xx_0_xyz_yyyz[k];

                g_xx_0_xyyz_yzz[k] = -g_xx_0_xyz_yzz[k] * ab_y + g_xx_0_xyz_yyzz[k];

                g_xx_0_xyyz_zzz[k] = -g_xx_0_xyz_zzz[k] * ab_y + g_xx_0_xyz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyzz_xxx, g_xx_0_xyzz_xxy, g_xx_0_xyzz_xxz, g_xx_0_xyzz_xyy, g_xx_0_xyzz_xyz, g_xx_0_xyzz_xzz, g_xx_0_xyzz_yyy, g_xx_0_xyzz_yyz, g_xx_0_xyzz_yzz, g_xx_0_xyzz_zzz, g_xx_0_xzz_xxx, g_xx_0_xzz_xxxy, g_xx_0_xzz_xxy, g_xx_0_xzz_xxyy, g_xx_0_xzz_xxyz, g_xx_0_xzz_xxz, g_xx_0_xzz_xyy, g_xx_0_xzz_xyyy, g_xx_0_xzz_xyyz, g_xx_0_xzz_xyz, g_xx_0_xzz_xyzz, g_xx_0_xzz_xzz, g_xx_0_xzz_yyy, g_xx_0_xzz_yyyy, g_xx_0_xzz_yyyz, g_xx_0_xzz_yyz, g_xx_0_xzz_yyzz, g_xx_0_xzz_yzz, g_xx_0_xzz_yzzz, g_xx_0_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyzz_xxx[k] = -g_xx_0_xzz_xxx[k] * ab_y + g_xx_0_xzz_xxxy[k];

                g_xx_0_xyzz_xxy[k] = -g_xx_0_xzz_xxy[k] * ab_y + g_xx_0_xzz_xxyy[k];

                g_xx_0_xyzz_xxz[k] = -g_xx_0_xzz_xxz[k] * ab_y + g_xx_0_xzz_xxyz[k];

                g_xx_0_xyzz_xyy[k] = -g_xx_0_xzz_xyy[k] * ab_y + g_xx_0_xzz_xyyy[k];

                g_xx_0_xyzz_xyz[k] = -g_xx_0_xzz_xyz[k] * ab_y + g_xx_0_xzz_xyyz[k];

                g_xx_0_xyzz_xzz[k] = -g_xx_0_xzz_xzz[k] * ab_y + g_xx_0_xzz_xyzz[k];

                g_xx_0_xyzz_yyy[k] = -g_xx_0_xzz_yyy[k] * ab_y + g_xx_0_xzz_yyyy[k];

                g_xx_0_xyzz_yyz[k] = -g_xx_0_xzz_yyz[k] * ab_y + g_xx_0_xzz_yyyz[k];

                g_xx_0_xyzz_yzz[k] = -g_xx_0_xzz_yzz[k] * ab_y + g_xx_0_xzz_yyzz[k];

                g_xx_0_xyzz_zzz[k] = -g_xx_0_xzz_zzz[k] * ab_y + g_xx_0_xzz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xzz_xxx, g_xx_0_xzz_xxxz, g_xx_0_xzz_xxy, g_xx_0_xzz_xxyz, g_xx_0_xzz_xxz, g_xx_0_xzz_xxzz, g_xx_0_xzz_xyy, g_xx_0_xzz_xyyz, g_xx_0_xzz_xyz, g_xx_0_xzz_xyzz, g_xx_0_xzz_xzz, g_xx_0_xzz_xzzz, g_xx_0_xzz_yyy, g_xx_0_xzz_yyyz, g_xx_0_xzz_yyz, g_xx_0_xzz_yyzz, g_xx_0_xzz_yzz, g_xx_0_xzz_yzzz, g_xx_0_xzz_zzz, g_xx_0_xzz_zzzz, g_xx_0_xzzz_xxx, g_xx_0_xzzz_xxy, g_xx_0_xzzz_xxz, g_xx_0_xzzz_xyy, g_xx_0_xzzz_xyz, g_xx_0_xzzz_xzz, g_xx_0_xzzz_yyy, g_xx_0_xzzz_yyz, g_xx_0_xzzz_yzz, g_xx_0_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzzz_xxx[k] = -g_xx_0_xzz_xxx[k] * ab_z + g_xx_0_xzz_xxxz[k];

                g_xx_0_xzzz_xxy[k] = -g_xx_0_xzz_xxy[k] * ab_z + g_xx_0_xzz_xxyz[k];

                g_xx_0_xzzz_xxz[k] = -g_xx_0_xzz_xxz[k] * ab_z + g_xx_0_xzz_xxzz[k];

                g_xx_0_xzzz_xyy[k] = -g_xx_0_xzz_xyy[k] * ab_z + g_xx_0_xzz_xyyz[k];

                g_xx_0_xzzz_xyz[k] = -g_xx_0_xzz_xyz[k] * ab_z + g_xx_0_xzz_xyzz[k];

                g_xx_0_xzzz_xzz[k] = -g_xx_0_xzz_xzz[k] * ab_z + g_xx_0_xzz_xzzz[k];

                g_xx_0_xzzz_yyy[k] = -g_xx_0_xzz_yyy[k] * ab_z + g_xx_0_xzz_yyyz[k];

                g_xx_0_xzzz_yyz[k] = -g_xx_0_xzz_yyz[k] * ab_z + g_xx_0_xzz_yyzz[k];

                g_xx_0_xzzz_yzz[k] = -g_xx_0_xzz_yzz[k] * ab_z + g_xx_0_xzz_yzzz[k];

                g_xx_0_xzzz_zzz[k] = -g_xx_0_xzz_zzz[k] * ab_z + g_xx_0_xzz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyy_xxx, g_xx_0_yyy_xxxy, g_xx_0_yyy_xxy, g_xx_0_yyy_xxyy, g_xx_0_yyy_xxyz, g_xx_0_yyy_xxz, g_xx_0_yyy_xyy, g_xx_0_yyy_xyyy, g_xx_0_yyy_xyyz, g_xx_0_yyy_xyz, g_xx_0_yyy_xyzz, g_xx_0_yyy_xzz, g_xx_0_yyy_yyy, g_xx_0_yyy_yyyy, g_xx_0_yyy_yyyz, g_xx_0_yyy_yyz, g_xx_0_yyy_yyzz, g_xx_0_yyy_yzz, g_xx_0_yyy_yzzz, g_xx_0_yyy_zzz, g_xx_0_yyyy_xxx, g_xx_0_yyyy_xxy, g_xx_0_yyyy_xxz, g_xx_0_yyyy_xyy, g_xx_0_yyyy_xyz, g_xx_0_yyyy_xzz, g_xx_0_yyyy_yyy, g_xx_0_yyyy_yyz, g_xx_0_yyyy_yzz, g_xx_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyy_xxx[k] = -g_xx_0_yyy_xxx[k] * ab_y + g_xx_0_yyy_xxxy[k];

                g_xx_0_yyyy_xxy[k] = -g_xx_0_yyy_xxy[k] * ab_y + g_xx_0_yyy_xxyy[k];

                g_xx_0_yyyy_xxz[k] = -g_xx_0_yyy_xxz[k] * ab_y + g_xx_0_yyy_xxyz[k];

                g_xx_0_yyyy_xyy[k] = -g_xx_0_yyy_xyy[k] * ab_y + g_xx_0_yyy_xyyy[k];

                g_xx_0_yyyy_xyz[k] = -g_xx_0_yyy_xyz[k] * ab_y + g_xx_0_yyy_xyyz[k];

                g_xx_0_yyyy_xzz[k] = -g_xx_0_yyy_xzz[k] * ab_y + g_xx_0_yyy_xyzz[k];

                g_xx_0_yyyy_yyy[k] = -g_xx_0_yyy_yyy[k] * ab_y + g_xx_0_yyy_yyyy[k];

                g_xx_0_yyyy_yyz[k] = -g_xx_0_yyy_yyz[k] * ab_y + g_xx_0_yyy_yyyz[k];

                g_xx_0_yyyy_yzz[k] = -g_xx_0_yyy_yzz[k] * ab_y + g_xx_0_yyy_yyzz[k];

                g_xx_0_yyyy_zzz[k] = -g_xx_0_yyy_zzz[k] * ab_y + g_xx_0_yyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyz_xxx, g_xx_0_yyyz_xxy, g_xx_0_yyyz_xxz, g_xx_0_yyyz_xyy, g_xx_0_yyyz_xyz, g_xx_0_yyyz_xzz, g_xx_0_yyyz_yyy, g_xx_0_yyyz_yyz, g_xx_0_yyyz_yzz, g_xx_0_yyyz_zzz, g_xx_0_yyz_xxx, g_xx_0_yyz_xxxy, g_xx_0_yyz_xxy, g_xx_0_yyz_xxyy, g_xx_0_yyz_xxyz, g_xx_0_yyz_xxz, g_xx_0_yyz_xyy, g_xx_0_yyz_xyyy, g_xx_0_yyz_xyyz, g_xx_0_yyz_xyz, g_xx_0_yyz_xyzz, g_xx_0_yyz_xzz, g_xx_0_yyz_yyy, g_xx_0_yyz_yyyy, g_xx_0_yyz_yyyz, g_xx_0_yyz_yyz, g_xx_0_yyz_yyzz, g_xx_0_yyz_yzz, g_xx_0_yyz_yzzz, g_xx_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyz_xxx[k] = -g_xx_0_yyz_xxx[k] * ab_y + g_xx_0_yyz_xxxy[k];

                g_xx_0_yyyz_xxy[k] = -g_xx_0_yyz_xxy[k] * ab_y + g_xx_0_yyz_xxyy[k];

                g_xx_0_yyyz_xxz[k] = -g_xx_0_yyz_xxz[k] * ab_y + g_xx_0_yyz_xxyz[k];

                g_xx_0_yyyz_xyy[k] = -g_xx_0_yyz_xyy[k] * ab_y + g_xx_0_yyz_xyyy[k];

                g_xx_0_yyyz_xyz[k] = -g_xx_0_yyz_xyz[k] * ab_y + g_xx_0_yyz_xyyz[k];

                g_xx_0_yyyz_xzz[k] = -g_xx_0_yyz_xzz[k] * ab_y + g_xx_0_yyz_xyzz[k];

                g_xx_0_yyyz_yyy[k] = -g_xx_0_yyz_yyy[k] * ab_y + g_xx_0_yyz_yyyy[k];

                g_xx_0_yyyz_yyz[k] = -g_xx_0_yyz_yyz[k] * ab_y + g_xx_0_yyz_yyyz[k];

                g_xx_0_yyyz_yzz[k] = -g_xx_0_yyz_yzz[k] * ab_y + g_xx_0_yyz_yyzz[k];

                g_xx_0_yyyz_zzz[k] = -g_xx_0_yyz_zzz[k] * ab_y + g_xx_0_yyz_yzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 125 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyzz_xxx, g_xx_0_yyzz_xxy, g_xx_0_yyzz_xxz, g_xx_0_yyzz_xyy, g_xx_0_yyzz_xyz, g_xx_0_yyzz_xzz, g_xx_0_yyzz_yyy, g_xx_0_yyzz_yyz, g_xx_0_yyzz_yzz, g_xx_0_yyzz_zzz, g_xx_0_yzz_xxx, g_xx_0_yzz_xxxy, g_xx_0_yzz_xxy, g_xx_0_yzz_xxyy, g_xx_0_yzz_xxyz, g_xx_0_yzz_xxz, g_xx_0_yzz_xyy, g_xx_0_yzz_xyyy, g_xx_0_yzz_xyyz, g_xx_0_yzz_xyz, g_xx_0_yzz_xyzz, g_xx_0_yzz_xzz, g_xx_0_yzz_yyy, g_xx_0_yzz_yyyy, g_xx_0_yzz_yyyz, g_xx_0_yzz_yyz, g_xx_0_yzz_yyzz, g_xx_0_yzz_yzz, g_xx_0_yzz_yzzz, g_xx_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyzz_xxx[k] = -g_xx_0_yzz_xxx[k] * ab_y + g_xx_0_yzz_xxxy[k];

                g_xx_0_yyzz_xxy[k] = -g_xx_0_yzz_xxy[k] * ab_y + g_xx_0_yzz_xxyy[k];

                g_xx_0_yyzz_xxz[k] = -g_xx_0_yzz_xxz[k] * ab_y + g_xx_0_yzz_xxyz[k];

                g_xx_0_yyzz_xyy[k] = -g_xx_0_yzz_xyy[k] * ab_y + g_xx_0_yzz_xyyy[k];

                g_xx_0_yyzz_xyz[k] = -g_xx_0_yzz_xyz[k] * ab_y + g_xx_0_yzz_xyyz[k];

                g_xx_0_yyzz_xzz[k] = -g_xx_0_yzz_xzz[k] * ab_y + g_xx_0_yzz_xyzz[k];

                g_xx_0_yyzz_yyy[k] = -g_xx_0_yzz_yyy[k] * ab_y + g_xx_0_yzz_yyyy[k];

                g_xx_0_yyzz_yyz[k] = -g_xx_0_yzz_yyz[k] * ab_y + g_xx_0_yzz_yyyz[k];

                g_xx_0_yyzz_yzz[k] = -g_xx_0_yzz_yzz[k] * ab_y + g_xx_0_yzz_yyzz[k];

                g_xx_0_yyzz_zzz[k] = -g_xx_0_yzz_zzz[k] * ab_y + g_xx_0_yzz_yzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 131 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 134 * ccomps * dcomps);

            auto g_xx_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 137 * ccomps * dcomps);

            auto g_xx_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzzz_xxx, g_xx_0_yzzz_xxy, g_xx_0_yzzz_xxz, g_xx_0_yzzz_xyy, g_xx_0_yzzz_xyz, g_xx_0_yzzz_xzz, g_xx_0_yzzz_yyy, g_xx_0_yzzz_yyz, g_xx_0_yzzz_yzz, g_xx_0_yzzz_zzz, g_xx_0_zzz_xxx, g_xx_0_zzz_xxxy, g_xx_0_zzz_xxy, g_xx_0_zzz_xxyy, g_xx_0_zzz_xxyz, g_xx_0_zzz_xxz, g_xx_0_zzz_xyy, g_xx_0_zzz_xyyy, g_xx_0_zzz_xyyz, g_xx_0_zzz_xyz, g_xx_0_zzz_xyzz, g_xx_0_zzz_xzz, g_xx_0_zzz_yyy, g_xx_0_zzz_yyyy, g_xx_0_zzz_yyyz, g_xx_0_zzz_yyz, g_xx_0_zzz_yyzz, g_xx_0_zzz_yzz, g_xx_0_zzz_yzzz, g_xx_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzzz_xxx[k] = -g_xx_0_zzz_xxx[k] * ab_y + g_xx_0_zzz_xxxy[k];

                g_xx_0_yzzz_xxy[k] = -g_xx_0_zzz_xxy[k] * ab_y + g_xx_0_zzz_xxyy[k];

                g_xx_0_yzzz_xxz[k] = -g_xx_0_zzz_xxz[k] * ab_y + g_xx_0_zzz_xxyz[k];

                g_xx_0_yzzz_xyy[k] = -g_xx_0_zzz_xyy[k] * ab_y + g_xx_0_zzz_xyyy[k];

                g_xx_0_yzzz_xyz[k] = -g_xx_0_zzz_xyz[k] * ab_y + g_xx_0_zzz_xyyz[k];

                g_xx_0_yzzz_xzz[k] = -g_xx_0_zzz_xzz[k] * ab_y + g_xx_0_zzz_xyzz[k];

                g_xx_0_yzzz_yyy[k] = -g_xx_0_zzz_yyy[k] * ab_y + g_xx_0_zzz_yyyy[k];

                g_xx_0_yzzz_yyz[k] = -g_xx_0_zzz_yyz[k] * ab_y + g_xx_0_zzz_yyyz[k];

                g_xx_0_yzzz_yzz[k] = -g_xx_0_zzz_yzz[k] * ab_y + g_xx_0_zzz_yyzz[k];

                g_xx_0_yzzz_zzz[k] = -g_xx_0_zzz_zzz[k] * ab_y + g_xx_0_zzz_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 143 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zzz_xxx, g_xx_0_zzz_xxxz, g_xx_0_zzz_xxy, g_xx_0_zzz_xxyz, g_xx_0_zzz_xxz, g_xx_0_zzz_xxzz, g_xx_0_zzz_xyy, g_xx_0_zzz_xyyz, g_xx_0_zzz_xyz, g_xx_0_zzz_xyzz, g_xx_0_zzz_xzz, g_xx_0_zzz_xzzz, g_xx_0_zzz_yyy, g_xx_0_zzz_yyyz, g_xx_0_zzz_yyz, g_xx_0_zzz_yyzz, g_xx_0_zzz_yzz, g_xx_0_zzz_yzzz, g_xx_0_zzz_zzz, g_xx_0_zzz_zzzz, g_xx_0_zzzz_xxx, g_xx_0_zzzz_xxy, g_xx_0_zzzz_xxz, g_xx_0_zzzz_xyy, g_xx_0_zzzz_xyz, g_xx_0_zzzz_xzz, g_xx_0_zzzz_yyy, g_xx_0_zzzz_yyz, g_xx_0_zzzz_yzz, g_xx_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzzz_xxx[k] = -g_xx_0_zzz_xxx[k] * ab_z + g_xx_0_zzz_xxxz[k];

                g_xx_0_zzzz_xxy[k] = -g_xx_0_zzz_xxy[k] * ab_z + g_xx_0_zzz_xxyz[k];

                g_xx_0_zzzz_xxz[k] = -g_xx_0_zzz_xxz[k] * ab_z + g_xx_0_zzz_xxzz[k];

                g_xx_0_zzzz_xyy[k] = -g_xx_0_zzz_xyy[k] * ab_z + g_xx_0_zzz_xyyz[k];

                g_xx_0_zzzz_xyz[k] = -g_xx_0_zzz_xyz[k] * ab_z + g_xx_0_zzz_xyzz[k];

                g_xx_0_zzzz_xzz[k] = -g_xx_0_zzz_xzz[k] * ab_z + g_xx_0_zzz_xzzz[k];

                g_xx_0_zzzz_yyy[k] = -g_xx_0_zzz_yyy[k] * ab_z + g_xx_0_zzz_yyyz[k];

                g_xx_0_zzzz_yyz[k] = -g_xx_0_zzz_yyz[k] * ab_z + g_xx_0_zzz_yyzz[k];

                g_xx_0_zzzz_yzz[k] = -g_xx_0_zzz_yzz[k] * ab_z + g_xx_0_zzz_yzzz[k];

                g_xx_0_zzzz_zzz[k] = -g_xx_0_zzz_zzz[k] * ab_z + g_xx_0_zzz_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxx_xxx, g_xy_0_xxx_xxxx, g_xy_0_xxx_xxxy, g_xy_0_xxx_xxxz, g_xy_0_xxx_xxy, g_xy_0_xxx_xxyy, g_xy_0_xxx_xxyz, g_xy_0_xxx_xxz, g_xy_0_xxx_xxzz, g_xy_0_xxx_xyy, g_xy_0_xxx_xyyy, g_xy_0_xxx_xyyz, g_xy_0_xxx_xyz, g_xy_0_xxx_xyzz, g_xy_0_xxx_xzz, g_xy_0_xxx_xzzz, g_xy_0_xxx_yyy, g_xy_0_xxx_yyz, g_xy_0_xxx_yzz, g_xy_0_xxx_zzz, g_xy_0_xxxx_xxx, g_xy_0_xxxx_xxy, g_xy_0_xxxx_xxz, g_xy_0_xxxx_xyy, g_xy_0_xxxx_xyz, g_xy_0_xxxx_xzz, g_xy_0_xxxx_yyy, g_xy_0_xxxx_yyz, g_xy_0_xxxx_yzz, g_xy_0_xxxx_zzz, g_y_0_xxx_xxx, g_y_0_xxx_xxy, g_y_0_xxx_xxz, g_y_0_xxx_xyy, g_y_0_xxx_xyz, g_y_0_xxx_xzz, g_y_0_xxx_yyy, g_y_0_xxx_yyz, g_y_0_xxx_yzz, g_y_0_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxx_xxx[k] = -g_y_0_xxx_xxx[k] - g_xy_0_xxx_xxx[k] * ab_x + g_xy_0_xxx_xxxx[k];

                g_xy_0_xxxx_xxy[k] = -g_y_0_xxx_xxy[k] - g_xy_0_xxx_xxy[k] * ab_x + g_xy_0_xxx_xxxy[k];

                g_xy_0_xxxx_xxz[k] = -g_y_0_xxx_xxz[k] - g_xy_0_xxx_xxz[k] * ab_x + g_xy_0_xxx_xxxz[k];

                g_xy_0_xxxx_xyy[k] = -g_y_0_xxx_xyy[k] - g_xy_0_xxx_xyy[k] * ab_x + g_xy_0_xxx_xxyy[k];

                g_xy_0_xxxx_xyz[k] = -g_y_0_xxx_xyz[k] - g_xy_0_xxx_xyz[k] * ab_x + g_xy_0_xxx_xxyz[k];

                g_xy_0_xxxx_xzz[k] = -g_y_0_xxx_xzz[k] - g_xy_0_xxx_xzz[k] * ab_x + g_xy_0_xxx_xxzz[k];

                g_xy_0_xxxx_yyy[k] = -g_y_0_xxx_yyy[k] - g_xy_0_xxx_yyy[k] * ab_x + g_xy_0_xxx_xyyy[k];

                g_xy_0_xxxx_yyz[k] = -g_y_0_xxx_yyz[k] - g_xy_0_xxx_yyz[k] * ab_x + g_xy_0_xxx_xyyz[k];

                g_xy_0_xxxx_yzz[k] = -g_y_0_xxx_yzz[k] - g_xy_0_xxx_yzz[k] * ab_x + g_xy_0_xxx_xyzz[k];

                g_xy_0_xxxx_zzz[k] = -g_y_0_xxx_zzz[k] - g_xy_0_xxx_zzz[k] * ab_x + g_xy_0_xxx_xzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxy_xxx, g_xy_0_xxxy_xxy, g_xy_0_xxxy_xxz, g_xy_0_xxxy_xyy, g_xy_0_xxxy_xyz, g_xy_0_xxxy_xzz, g_xy_0_xxxy_yyy, g_xy_0_xxxy_yyz, g_xy_0_xxxy_yzz, g_xy_0_xxxy_zzz, g_xy_0_xxy_xxx, g_xy_0_xxy_xxxx, g_xy_0_xxy_xxxy, g_xy_0_xxy_xxxz, g_xy_0_xxy_xxy, g_xy_0_xxy_xxyy, g_xy_0_xxy_xxyz, g_xy_0_xxy_xxz, g_xy_0_xxy_xxzz, g_xy_0_xxy_xyy, g_xy_0_xxy_xyyy, g_xy_0_xxy_xyyz, g_xy_0_xxy_xyz, g_xy_0_xxy_xyzz, g_xy_0_xxy_xzz, g_xy_0_xxy_xzzz, g_xy_0_xxy_yyy, g_xy_0_xxy_yyz, g_xy_0_xxy_yzz, g_xy_0_xxy_zzz, g_y_0_xxy_xxx, g_y_0_xxy_xxy, g_y_0_xxy_xxz, g_y_0_xxy_xyy, g_y_0_xxy_xyz, g_y_0_xxy_xzz, g_y_0_xxy_yyy, g_y_0_xxy_yyz, g_y_0_xxy_yzz, g_y_0_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxy_xxx[k] = -g_y_0_xxy_xxx[k] - g_xy_0_xxy_xxx[k] * ab_x + g_xy_0_xxy_xxxx[k];

                g_xy_0_xxxy_xxy[k] = -g_y_0_xxy_xxy[k] - g_xy_0_xxy_xxy[k] * ab_x + g_xy_0_xxy_xxxy[k];

                g_xy_0_xxxy_xxz[k] = -g_y_0_xxy_xxz[k] - g_xy_0_xxy_xxz[k] * ab_x + g_xy_0_xxy_xxxz[k];

                g_xy_0_xxxy_xyy[k] = -g_y_0_xxy_xyy[k] - g_xy_0_xxy_xyy[k] * ab_x + g_xy_0_xxy_xxyy[k];

                g_xy_0_xxxy_xyz[k] = -g_y_0_xxy_xyz[k] - g_xy_0_xxy_xyz[k] * ab_x + g_xy_0_xxy_xxyz[k];

                g_xy_0_xxxy_xzz[k] = -g_y_0_xxy_xzz[k] - g_xy_0_xxy_xzz[k] * ab_x + g_xy_0_xxy_xxzz[k];

                g_xy_0_xxxy_yyy[k] = -g_y_0_xxy_yyy[k] - g_xy_0_xxy_yyy[k] * ab_x + g_xy_0_xxy_xyyy[k];

                g_xy_0_xxxy_yyz[k] = -g_y_0_xxy_yyz[k] - g_xy_0_xxy_yyz[k] * ab_x + g_xy_0_xxy_xyyz[k];

                g_xy_0_xxxy_yzz[k] = -g_y_0_xxy_yzz[k] - g_xy_0_xxy_yzz[k] * ab_x + g_xy_0_xxy_xyzz[k];

                g_xy_0_xxxy_zzz[k] = -g_y_0_xxy_zzz[k] - g_xy_0_xxy_zzz[k] * ab_x + g_xy_0_xxy_xzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxx_xxx, g_xy_0_xxx_xxxz, g_xy_0_xxx_xxy, g_xy_0_xxx_xxyz, g_xy_0_xxx_xxz, g_xy_0_xxx_xxzz, g_xy_0_xxx_xyy, g_xy_0_xxx_xyyz, g_xy_0_xxx_xyz, g_xy_0_xxx_xyzz, g_xy_0_xxx_xzz, g_xy_0_xxx_xzzz, g_xy_0_xxx_yyy, g_xy_0_xxx_yyyz, g_xy_0_xxx_yyz, g_xy_0_xxx_yyzz, g_xy_0_xxx_yzz, g_xy_0_xxx_yzzz, g_xy_0_xxx_zzz, g_xy_0_xxx_zzzz, g_xy_0_xxxz_xxx, g_xy_0_xxxz_xxy, g_xy_0_xxxz_xxz, g_xy_0_xxxz_xyy, g_xy_0_xxxz_xyz, g_xy_0_xxxz_xzz, g_xy_0_xxxz_yyy, g_xy_0_xxxz_yyz, g_xy_0_xxxz_yzz, g_xy_0_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxz_xxx[k] = -g_xy_0_xxx_xxx[k] * ab_z + g_xy_0_xxx_xxxz[k];

                g_xy_0_xxxz_xxy[k] = -g_xy_0_xxx_xxy[k] * ab_z + g_xy_0_xxx_xxyz[k];

                g_xy_0_xxxz_xxz[k] = -g_xy_0_xxx_xxz[k] * ab_z + g_xy_0_xxx_xxzz[k];

                g_xy_0_xxxz_xyy[k] = -g_xy_0_xxx_xyy[k] * ab_z + g_xy_0_xxx_xyyz[k];

                g_xy_0_xxxz_xyz[k] = -g_xy_0_xxx_xyz[k] * ab_z + g_xy_0_xxx_xyzz[k];

                g_xy_0_xxxz_xzz[k] = -g_xy_0_xxx_xzz[k] * ab_z + g_xy_0_xxx_xzzz[k];

                g_xy_0_xxxz_yyy[k] = -g_xy_0_xxx_yyy[k] * ab_z + g_xy_0_xxx_yyyz[k];

                g_xy_0_xxxz_yyz[k] = -g_xy_0_xxx_yyz[k] * ab_z + g_xy_0_xxx_yyzz[k];

                g_xy_0_xxxz_yzz[k] = -g_xy_0_xxx_yzz[k] * ab_z + g_xy_0_xxx_yzzz[k];

                g_xy_0_xxxz_zzz[k] = -g_xy_0_xxx_zzz[k] * ab_z + g_xy_0_xxx_zzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyy_xxx, g_xy_0_xxyy_xxy, g_xy_0_xxyy_xxz, g_xy_0_xxyy_xyy, g_xy_0_xxyy_xyz, g_xy_0_xxyy_xzz, g_xy_0_xxyy_yyy, g_xy_0_xxyy_yyz, g_xy_0_xxyy_yzz, g_xy_0_xxyy_zzz, g_xy_0_xyy_xxx, g_xy_0_xyy_xxxx, g_xy_0_xyy_xxxy, g_xy_0_xyy_xxxz, g_xy_0_xyy_xxy, g_xy_0_xyy_xxyy, g_xy_0_xyy_xxyz, g_xy_0_xyy_xxz, g_xy_0_xyy_xxzz, g_xy_0_xyy_xyy, g_xy_0_xyy_xyyy, g_xy_0_xyy_xyyz, g_xy_0_xyy_xyz, g_xy_0_xyy_xyzz, g_xy_0_xyy_xzz, g_xy_0_xyy_xzzz, g_xy_0_xyy_yyy, g_xy_0_xyy_yyz, g_xy_0_xyy_yzz, g_xy_0_xyy_zzz, g_y_0_xyy_xxx, g_y_0_xyy_xxy, g_y_0_xyy_xxz, g_y_0_xyy_xyy, g_y_0_xyy_xyz, g_y_0_xyy_xzz, g_y_0_xyy_yyy, g_y_0_xyy_yyz, g_y_0_xyy_yzz, g_y_0_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyy_xxx[k] = -g_y_0_xyy_xxx[k] - g_xy_0_xyy_xxx[k] * ab_x + g_xy_0_xyy_xxxx[k];

                g_xy_0_xxyy_xxy[k] = -g_y_0_xyy_xxy[k] - g_xy_0_xyy_xxy[k] * ab_x + g_xy_0_xyy_xxxy[k];

                g_xy_0_xxyy_xxz[k] = -g_y_0_xyy_xxz[k] - g_xy_0_xyy_xxz[k] * ab_x + g_xy_0_xyy_xxxz[k];

                g_xy_0_xxyy_xyy[k] = -g_y_0_xyy_xyy[k] - g_xy_0_xyy_xyy[k] * ab_x + g_xy_0_xyy_xxyy[k];

                g_xy_0_xxyy_xyz[k] = -g_y_0_xyy_xyz[k] - g_xy_0_xyy_xyz[k] * ab_x + g_xy_0_xyy_xxyz[k];

                g_xy_0_xxyy_xzz[k] = -g_y_0_xyy_xzz[k] - g_xy_0_xyy_xzz[k] * ab_x + g_xy_0_xyy_xxzz[k];

                g_xy_0_xxyy_yyy[k] = -g_y_0_xyy_yyy[k] - g_xy_0_xyy_yyy[k] * ab_x + g_xy_0_xyy_xyyy[k];

                g_xy_0_xxyy_yyz[k] = -g_y_0_xyy_yyz[k] - g_xy_0_xyy_yyz[k] * ab_x + g_xy_0_xyy_xyyz[k];

                g_xy_0_xxyy_yzz[k] = -g_y_0_xyy_yzz[k] - g_xy_0_xyy_yzz[k] * ab_x + g_xy_0_xyy_xyzz[k];

                g_xy_0_xxyy_zzz[k] = -g_y_0_xyy_zzz[k] - g_xy_0_xyy_zzz[k] * ab_x + g_xy_0_xyy_xzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxy_xxx, g_xy_0_xxy_xxxz, g_xy_0_xxy_xxy, g_xy_0_xxy_xxyz, g_xy_0_xxy_xxz, g_xy_0_xxy_xxzz, g_xy_0_xxy_xyy, g_xy_0_xxy_xyyz, g_xy_0_xxy_xyz, g_xy_0_xxy_xyzz, g_xy_0_xxy_xzz, g_xy_0_xxy_xzzz, g_xy_0_xxy_yyy, g_xy_0_xxy_yyyz, g_xy_0_xxy_yyz, g_xy_0_xxy_yyzz, g_xy_0_xxy_yzz, g_xy_0_xxy_yzzz, g_xy_0_xxy_zzz, g_xy_0_xxy_zzzz, g_xy_0_xxyz_xxx, g_xy_0_xxyz_xxy, g_xy_0_xxyz_xxz, g_xy_0_xxyz_xyy, g_xy_0_xxyz_xyz, g_xy_0_xxyz_xzz, g_xy_0_xxyz_yyy, g_xy_0_xxyz_yyz, g_xy_0_xxyz_yzz, g_xy_0_xxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyz_xxx[k] = -g_xy_0_xxy_xxx[k] * ab_z + g_xy_0_xxy_xxxz[k];

                g_xy_0_xxyz_xxy[k] = -g_xy_0_xxy_xxy[k] * ab_z + g_xy_0_xxy_xxyz[k];

                g_xy_0_xxyz_xxz[k] = -g_xy_0_xxy_xxz[k] * ab_z + g_xy_0_xxy_xxzz[k];

                g_xy_0_xxyz_xyy[k] = -g_xy_0_xxy_xyy[k] * ab_z + g_xy_0_xxy_xyyz[k];

                g_xy_0_xxyz_xyz[k] = -g_xy_0_xxy_xyz[k] * ab_z + g_xy_0_xxy_xyzz[k];

                g_xy_0_xxyz_xzz[k] = -g_xy_0_xxy_xzz[k] * ab_z + g_xy_0_xxy_xzzz[k];

                g_xy_0_xxyz_yyy[k] = -g_xy_0_xxy_yyy[k] * ab_z + g_xy_0_xxy_yyyz[k];

                g_xy_0_xxyz_yyz[k] = -g_xy_0_xxy_yyz[k] * ab_z + g_xy_0_xxy_yyzz[k];

                g_xy_0_xxyz_yzz[k] = -g_xy_0_xxy_yzz[k] * ab_z + g_xy_0_xxy_yzzz[k];

                g_xy_0_xxyz_zzz[k] = -g_xy_0_xxy_zzz[k] * ab_z + g_xy_0_xxy_zzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxz_xxx, g_xy_0_xxz_xxxz, g_xy_0_xxz_xxy, g_xy_0_xxz_xxyz, g_xy_0_xxz_xxz, g_xy_0_xxz_xxzz, g_xy_0_xxz_xyy, g_xy_0_xxz_xyyz, g_xy_0_xxz_xyz, g_xy_0_xxz_xyzz, g_xy_0_xxz_xzz, g_xy_0_xxz_xzzz, g_xy_0_xxz_yyy, g_xy_0_xxz_yyyz, g_xy_0_xxz_yyz, g_xy_0_xxz_yyzz, g_xy_0_xxz_yzz, g_xy_0_xxz_yzzz, g_xy_0_xxz_zzz, g_xy_0_xxz_zzzz, g_xy_0_xxzz_xxx, g_xy_0_xxzz_xxy, g_xy_0_xxzz_xxz, g_xy_0_xxzz_xyy, g_xy_0_xxzz_xyz, g_xy_0_xxzz_xzz, g_xy_0_xxzz_yyy, g_xy_0_xxzz_yyz, g_xy_0_xxzz_yzz, g_xy_0_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxzz_xxx[k] = -g_xy_0_xxz_xxx[k] * ab_z + g_xy_0_xxz_xxxz[k];

                g_xy_0_xxzz_xxy[k] = -g_xy_0_xxz_xxy[k] * ab_z + g_xy_0_xxz_xxyz[k];

                g_xy_0_xxzz_xxz[k] = -g_xy_0_xxz_xxz[k] * ab_z + g_xy_0_xxz_xxzz[k];

                g_xy_0_xxzz_xyy[k] = -g_xy_0_xxz_xyy[k] * ab_z + g_xy_0_xxz_xyyz[k];

                g_xy_0_xxzz_xyz[k] = -g_xy_0_xxz_xyz[k] * ab_z + g_xy_0_xxz_xyzz[k];

                g_xy_0_xxzz_xzz[k] = -g_xy_0_xxz_xzz[k] * ab_z + g_xy_0_xxz_xzzz[k];

                g_xy_0_xxzz_yyy[k] = -g_xy_0_xxz_yyy[k] * ab_z + g_xy_0_xxz_yyyz[k];

                g_xy_0_xxzz_yyz[k] = -g_xy_0_xxz_yyz[k] * ab_z + g_xy_0_xxz_yyzz[k];

                g_xy_0_xxzz_yzz[k] = -g_xy_0_xxz_yzz[k] * ab_z + g_xy_0_xxz_yzzz[k];

                g_xy_0_xxzz_zzz[k] = -g_xy_0_xxz_zzz[k] * ab_z + g_xy_0_xxz_zzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyy_xxx, g_xy_0_xyyy_xxy, g_xy_0_xyyy_xxz, g_xy_0_xyyy_xyy, g_xy_0_xyyy_xyz, g_xy_0_xyyy_xzz, g_xy_0_xyyy_yyy, g_xy_0_xyyy_yyz, g_xy_0_xyyy_yzz, g_xy_0_xyyy_zzz, g_xy_0_yyy_xxx, g_xy_0_yyy_xxxx, g_xy_0_yyy_xxxy, g_xy_0_yyy_xxxz, g_xy_0_yyy_xxy, g_xy_0_yyy_xxyy, g_xy_0_yyy_xxyz, g_xy_0_yyy_xxz, g_xy_0_yyy_xxzz, g_xy_0_yyy_xyy, g_xy_0_yyy_xyyy, g_xy_0_yyy_xyyz, g_xy_0_yyy_xyz, g_xy_0_yyy_xyzz, g_xy_0_yyy_xzz, g_xy_0_yyy_xzzz, g_xy_0_yyy_yyy, g_xy_0_yyy_yyz, g_xy_0_yyy_yzz, g_xy_0_yyy_zzz, g_y_0_yyy_xxx, g_y_0_yyy_xxy, g_y_0_yyy_xxz, g_y_0_yyy_xyy, g_y_0_yyy_xyz, g_y_0_yyy_xzz, g_y_0_yyy_yyy, g_y_0_yyy_yyz, g_y_0_yyy_yzz, g_y_0_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyy_xxx[k] = -g_y_0_yyy_xxx[k] - g_xy_0_yyy_xxx[k] * ab_x + g_xy_0_yyy_xxxx[k];

                g_xy_0_xyyy_xxy[k] = -g_y_0_yyy_xxy[k] - g_xy_0_yyy_xxy[k] * ab_x + g_xy_0_yyy_xxxy[k];

                g_xy_0_xyyy_xxz[k] = -g_y_0_yyy_xxz[k] - g_xy_0_yyy_xxz[k] * ab_x + g_xy_0_yyy_xxxz[k];

                g_xy_0_xyyy_xyy[k] = -g_y_0_yyy_xyy[k] - g_xy_0_yyy_xyy[k] * ab_x + g_xy_0_yyy_xxyy[k];

                g_xy_0_xyyy_xyz[k] = -g_y_0_yyy_xyz[k] - g_xy_0_yyy_xyz[k] * ab_x + g_xy_0_yyy_xxyz[k];

                g_xy_0_xyyy_xzz[k] = -g_y_0_yyy_xzz[k] - g_xy_0_yyy_xzz[k] * ab_x + g_xy_0_yyy_xxzz[k];

                g_xy_0_xyyy_yyy[k] = -g_y_0_yyy_yyy[k] - g_xy_0_yyy_yyy[k] * ab_x + g_xy_0_yyy_xyyy[k];

                g_xy_0_xyyy_yyz[k] = -g_y_0_yyy_yyz[k] - g_xy_0_yyy_yyz[k] * ab_x + g_xy_0_yyy_xyyz[k];

                g_xy_0_xyyy_yzz[k] = -g_y_0_yyy_yzz[k] - g_xy_0_yyy_yzz[k] * ab_x + g_xy_0_yyy_xyzz[k];

                g_xy_0_xyyy_zzz[k] = -g_y_0_yyy_zzz[k] - g_xy_0_yyy_zzz[k] * ab_x + g_xy_0_yyy_xzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyy_xxx, g_xy_0_xyy_xxxz, g_xy_0_xyy_xxy, g_xy_0_xyy_xxyz, g_xy_0_xyy_xxz, g_xy_0_xyy_xxzz, g_xy_0_xyy_xyy, g_xy_0_xyy_xyyz, g_xy_0_xyy_xyz, g_xy_0_xyy_xyzz, g_xy_0_xyy_xzz, g_xy_0_xyy_xzzz, g_xy_0_xyy_yyy, g_xy_0_xyy_yyyz, g_xy_0_xyy_yyz, g_xy_0_xyy_yyzz, g_xy_0_xyy_yzz, g_xy_0_xyy_yzzz, g_xy_0_xyy_zzz, g_xy_0_xyy_zzzz, g_xy_0_xyyz_xxx, g_xy_0_xyyz_xxy, g_xy_0_xyyz_xxz, g_xy_0_xyyz_xyy, g_xy_0_xyyz_xyz, g_xy_0_xyyz_xzz, g_xy_0_xyyz_yyy, g_xy_0_xyyz_yyz, g_xy_0_xyyz_yzz, g_xy_0_xyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyz_xxx[k] = -g_xy_0_xyy_xxx[k] * ab_z + g_xy_0_xyy_xxxz[k];

                g_xy_0_xyyz_xxy[k] = -g_xy_0_xyy_xxy[k] * ab_z + g_xy_0_xyy_xxyz[k];

                g_xy_0_xyyz_xxz[k] = -g_xy_0_xyy_xxz[k] * ab_z + g_xy_0_xyy_xxzz[k];

                g_xy_0_xyyz_xyy[k] = -g_xy_0_xyy_xyy[k] * ab_z + g_xy_0_xyy_xyyz[k];

                g_xy_0_xyyz_xyz[k] = -g_xy_0_xyy_xyz[k] * ab_z + g_xy_0_xyy_xyzz[k];

                g_xy_0_xyyz_xzz[k] = -g_xy_0_xyy_xzz[k] * ab_z + g_xy_0_xyy_xzzz[k];

                g_xy_0_xyyz_yyy[k] = -g_xy_0_xyy_yyy[k] * ab_z + g_xy_0_xyy_yyyz[k];

                g_xy_0_xyyz_yyz[k] = -g_xy_0_xyy_yyz[k] * ab_z + g_xy_0_xyy_yyzz[k];

                g_xy_0_xyyz_yzz[k] = -g_xy_0_xyy_yzz[k] * ab_z + g_xy_0_xyy_yzzz[k];

                g_xy_0_xyyz_zzz[k] = -g_xy_0_xyy_zzz[k] * ab_z + g_xy_0_xyy_zzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyz_xxx, g_xy_0_xyz_xxxz, g_xy_0_xyz_xxy, g_xy_0_xyz_xxyz, g_xy_0_xyz_xxz, g_xy_0_xyz_xxzz, g_xy_0_xyz_xyy, g_xy_0_xyz_xyyz, g_xy_0_xyz_xyz, g_xy_0_xyz_xyzz, g_xy_0_xyz_xzz, g_xy_0_xyz_xzzz, g_xy_0_xyz_yyy, g_xy_0_xyz_yyyz, g_xy_0_xyz_yyz, g_xy_0_xyz_yyzz, g_xy_0_xyz_yzz, g_xy_0_xyz_yzzz, g_xy_0_xyz_zzz, g_xy_0_xyz_zzzz, g_xy_0_xyzz_xxx, g_xy_0_xyzz_xxy, g_xy_0_xyzz_xxz, g_xy_0_xyzz_xyy, g_xy_0_xyzz_xyz, g_xy_0_xyzz_xzz, g_xy_0_xyzz_yyy, g_xy_0_xyzz_yyz, g_xy_0_xyzz_yzz, g_xy_0_xyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyzz_xxx[k] = -g_xy_0_xyz_xxx[k] * ab_z + g_xy_0_xyz_xxxz[k];

                g_xy_0_xyzz_xxy[k] = -g_xy_0_xyz_xxy[k] * ab_z + g_xy_0_xyz_xxyz[k];

                g_xy_0_xyzz_xxz[k] = -g_xy_0_xyz_xxz[k] * ab_z + g_xy_0_xyz_xxzz[k];

                g_xy_0_xyzz_xyy[k] = -g_xy_0_xyz_xyy[k] * ab_z + g_xy_0_xyz_xyyz[k];

                g_xy_0_xyzz_xyz[k] = -g_xy_0_xyz_xyz[k] * ab_z + g_xy_0_xyz_xyzz[k];

                g_xy_0_xyzz_xzz[k] = -g_xy_0_xyz_xzz[k] * ab_z + g_xy_0_xyz_xzzz[k];

                g_xy_0_xyzz_yyy[k] = -g_xy_0_xyz_yyy[k] * ab_z + g_xy_0_xyz_yyyz[k];

                g_xy_0_xyzz_yyz[k] = -g_xy_0_xyz_yyz[k] * ab_z + g_xy_0_xyz_yyzz[k];

                g_xy_0_xyzz_yzz[k] = -g_xy_0_xyz_yzz[k] * ab_z + g_xy_0_xyz_yzzz[k];

                g_xy_0_xyzz_zzz[k] = -g_xy_0_xyz_zzz[k] * ab_z + g_xy_0_xyz_zzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xzz_xxx, g_xy_0_xzz_xxxz, g_xy_0_xzz_xxy, g_xy_0_xzz_xxyz, g_xy_0_xzz_xxz, g_xy_0_xzz_xxzz, g_xy_0_xzz_xyy, g_xy_0_xzz_xyyz, g_xy_0_xzz_xyz, g_xy_0_xzz_xyzz, g_xy_0_xzz_xzz, g_xy_0_xzz_xzzz, g_xy_0_xzz_yyy, g_xy_0_xzz_yyyz, g_xy_0_xzz_yyz, g_xy_0_xzz_yyzz, g_xy_0_xzz_yzz, g_xy_0_xzz_yzzz, g_xy_0_xzz_zzz, g_xy_0_xzz_zzzz, g_xy_0_xzzz_xxx, g_xy_0_xzzz_xxy, g_xy_0_xzzz_xxz, g_xy_0_xzzz_xyy, g_xy_0_xzzz_xyz, g_xy_0_xzzz_xzz, g_xy_0_xzzz_yyy, g_xy_0_xzzz_yyz, g_xy_0_xzzz_yzz, g_xy_0_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzzz_xxx[k] = -g_xy_0_xzz_xxx[k] * ab_z + g_xy_0_xzz_xxxz[k];

                g_xy_0_xzzz_xxy[k] = -g_xy_0_xzz_xxy[k] * ab_z + g_xy_0_xzz_xxyz[k];

                g_xy_0_xzzz_xxz[k] = -g_xy_0_xzz_xxz[k] * ab_z + g_xy_0_xzz_xxzz[k];

                g_xy_0_xzzz_xyy[k] = -g_xy_0_xzz_xyy[k] * ab_z + g_xy_0_xzz_xyyz[k];

                g_xy_0_xzzz_xyz[k] = -g_xy_0_xzz_xyz[k] * ab_z + g_xy_0_xzz_xyzz[k];

                g_xy_0_xzzz_xzz[k] = -g_xy_0_xzz_xzz[k] * ab_z + g_xy_0_xzz_xzzz[k];

                g_xy_0_xzzz_yyy[k] = -g_xy_0_xzz_yyy[k] * ab_z + g_xy_0_xzz_yyyz[k];

                g_xy_0_xzzz_yyz[k] = -g_xy_0_xzz_yyz[k] * ab_z + g_xy_0_xzz_yyzz[k];

                g_xy_0_xzzz_yzz[k] = -g_xy_0_xzz_yzz[k] * ab_z + g_xy_0_xzz_yzzz[k];

                g_xy_0_xzzz_zzz[k] = -g_xy_0_xzz_zzz[k] * ab_z + g_xy_0_xzz_zzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 251 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 252 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 253 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 254 * ccomps * dcomps);

            auto g_xy_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 255 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 256 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 257 * ccomps * dcomps);

            auto g_xy_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 258 * ccomps * dcomps);

            auto g_xy_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 259 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyy_xxx, g_x_0_yyy_xxy, g_x_0_yyy_xxz, g_x_0_yyy_xyy, g_x_0_yyy_xyz, g_x_0_yyy_xzz, g_x_0_yyy_yyy, g_x_0_yyy_yyz, g_x_0_yyy_yzz, g_x_0_yyy_zzz, g_xy_0_yyy_xxx, g_xy_0_yyy_xxxy, g_xy_0_yyy_xxy, g_xy_0_yyy_xxyy, g_xy_0_yyy_xxyz, g_xy_0_yyy_xxz, g_xy_0_yyy_xyy, g_xy_0_yyy_xyyy, g_xy_0_yyy_xyyz, g_xy_0_yyy_xyz, g_xy_0_yyy_xyzz, g_xy_0_yyy_xzz, g_xy_0_yyy_yyy, g_xy_0_yyy_yyyy, g_xy_0_yyy_yyyz, g_xy_0_yyy_yyz, g_xy_0_yyy_yyzz, g_xy_0_yyy_yzz, g_xy_0_yyy_yzzz, g_xy_0_yyy_zzz, g_xy_0_yyyy_xxx, g_xy_0_yyyy_xxy, g_xy_0_yyyy_xxz, g_xy_0_yyyy_xyy, g_xy_0_yyyy_xyz, g_xy_0_yyyy_xzz, g_xy_0_yyyy_yyy, g_xy_0_yyyy_yyz, g_xy_0_yyyy_yzz, g_xy_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyy_xxx[k] = -g_x_0_yyy_xxx[k] - g_xy_0_yyy_xxx[k] * ab_y + g_xy_0_yyy_xxxy[k];

                g_xy_0_yyyy_xxy[k] = -g_x_0_yyy_xxy[k] - g_xy_0_yyy_xxy[k] * ab_y + g_xy_0_yyy_xxyy[k];

                g_xy_0_yyyy_xxz[k] = -g_x_0_yyy_xxz[k] - g_xy_0_yyy_xxz[k] * ab_y + g_xy_0_yyy_xxyz[k];

                g_xy_0_yyyy_xyy[k] = -g_x_0_yyy_xyy[k] - g_xy_0_yyy_xyy[k] * ab_y + g_xy_0_yyy_xyyy[k];

                g_xy_0_yyyy_xyz[k] = -g_x_0_yyy_xyz[k] - g_xy_0_yyy_xyz[k] * ab_y + g_xy_0_yyy_xyyz[k];

                g_xy_0_yyyy_xzz[k] = -g_x_0_yyy_xzz[k] - g_xy_0_yyy_xzz[k] * ab_y + g_xy_0_yyy_xyzz[k];

                g_xy_0_yyyy_yyy[k] = -g_x_0_yyy_yyy[k] - g_xy_0_yyy_yyy[k] * ab_y + g_xy_0_yyy_yyyy[k];

                g_xy_0_yyyy_yyz[k] = -g_x_0_yyy_yyz[k] - g_xy_0_yyy_yyz[k] * ab_y + g_xy_0_yyy_yyyz[k];

                g_xy_0_yyyy_yzz[k] = -g_x_0_yyy_yzz[k] - g_xy_0_yyy_yzz[k] * ab_y + g_xy_0_yyy_yyzz[k];

                g_xy_0_yyyy_zzz[k] = -g_x_0_yyy_zzz[k] - g_xy_0_yyy_zzz[k] * ab_y + g_xy_0_yyy_yzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 260 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 261 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 262 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 263 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 264 * ccomps * dcomps);

            auto g_xy_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 265 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 266 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 267 * ccomps * dcomps);

            auto g_xy_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 268 * ccomps * dcomps);

            auto g_xy_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyy_xxx, g_xy_0_yyy_xxxz, g_xy_0_yyy_xxy, g_xy_0_yyy_xxyz, g_xy_0_yyy_xxz, g_xy_0_yyy_xxzz, g_xy_0_yyy_xyy, g_xy_0_yyy_xyyz, g_xy_0_yyy_xyz, g_xy_0_yyy_xyzz, g_xy_0_yyy_xzz, g_xy_0_yyy_xzzz, g_xy_0_yyy_yyy, g_xy_0_yyy_yyyz, g_xy_0_yyy_yyz, g_xy_0_yyy_yyzz, g_xy_0_yyy_yzz, g_xy_0_yyy_yzzz, g_xy_0_yyy_zzz, g_xy_0_yyy_zzzz, g_xy_0_yyyz_xxx, g_xy_0_yyyz_xxy, g_xy_0_yyyz_xxz, g_xy_0_yyyz_xyy, g_xy_0_yyyz_xyz, g_xy_0_yyyz_xzz, g_xy_0_yyyz_yyy, g_xy_0_yyyz_yyz, g_xy_0_yyyz_yzz, g_xy_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyz_xxx[k] = -g_xy_0_yyy_xxx[k] * ab_z + g_xy_0_yyy_xxxz[k];

                g_xy_0_yyyz_xxy[k] = -g_xy_0_yyy_xxy[k] * ab_z + g_xy_0_yyy_xxyz[k];

                g_xy_0_yyyz_xxz[k] = -g_xy_0_yyy_xxz[k] * ab_z + g_xy_0_yyy_xxzz[k];

                g_xy_0_yyyz_xyy[k] = -g_xy_0_yyy_xyy[k] * ab_z + g_xy_0_yyy_xyyz[k];

                g_xy_0_yyyz_xyz[k] = -g_xy_0_yyy_xyz[k] * ab_z + g_xy_0_yyy_xyzz[k];

                g_xy_0_yyyz_xzz[k] = -g_xy_0_yyy_xzz[k] * ab_z + g_xy_0_yyy_xzzz[k];

                g_xy_0_yyyz_yyy[k] = -g_xy_0_yyy_yyy[k] * ab_z + g_xy_0_yyy_yyyz[k];

                g_xy_0_yyyz_yyz[k] = -g_xy_0_yyy_yyz[k] * ab_z + g_xy_0_yyy_yyzz[k];

                g_xy_0_yyyz_yzz[k] = -g_xy_0_yyy_yzz[k] * ab_z + g_xy_0_yyy_yzzz[k];

                g_xy_0_yyyz_zzz[k] = -g_xy_0_yyy_zzz[k] * ab_z + g_xy_0_yyy_zzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 270 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 271 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 272 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 273 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 274 * ccomps * dcomps);

            auto g_xy_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 275 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 276 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 277 * ccomps * dcomps);

            auto g_xy_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 278 * ccomps * dcomps);

            auto g_xy_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyz_xxx, g_xy_0_yyz_xxxz, g_xy_0_yyz_xxy, g_xy_0_yyz_xxyz, g_xy_0_yyz_xxz, g_xy_0_yyz_xxzz, g_xy_0_yyz_xyy, g_xy_0_yyz_xyyz, g_xy_0_yyz_xyz, g_xy_0_yyz_xyzz, g_xy_0_yyz_xzz, g_xy_0_yyz_xzzz, g_xy_0_yyz_yyy, g_xy_0_yyz_yyyz, g_xy_0_yyz_yyz, g_xy_0_yyz_yyzz, g_xy_0_yyz_yzz, g_xy_0_yyz_yzzz, g_xy_0_yyz_zzz, g_xy_0_yyz_zzzz, g_xy_0_yyzz_xxx, g_xy_0_yyzz_xxy, g_xy_0_yyzz_xxz, g_xy_0_yyzz_xyy, g_xy_0_yyzz_xyz, g_xy_0_yyzz_xzz, g_xy_0_yyzz_yyy, g_xy_0_yyzz_yyz, g_xy_0_yyzz_yzz, g_xy_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyzz_xxx[k] = -g_xy_0_yyz_xxx[k] * ab_z + g_xy_0_yyz_xxxz[k];

                g_xy_0_yyzz_xxy[k] = -g_xy_0_yyz_xxy[k] * ab_z + g_xy_0_yyz_xxyz[k];

                g_xy_0_yyzz_xxz[k] = -g_xy_0_yyz_xxz[k] * ab_z + g_xy_0_yyz_xxzz[k];

                g_xy_0_yyzz_xyy[k] = -g_xy_0_yyz_xyy[k] * ab_z + g_xy_0_yyz_xyyz[k];

                g_xy_0_yyzz_xyz[k] = -g_xy_0_yyz_xyz[k] * ab_z + g_xy_0_yyz_xyzz[k];

                g_xy_0_yyzz_xzz[k] = -g_xy_0_yyz_xzz[k] * ab_z + g_xy_0_yyz_xzzz[k];

                g_xy_0_yyzz_yyy[k] = -g_xy_0_yyz_yyy[k] * ab_z + g_xy_0_yyz_yyyz[k];

                g_xy_0_yyzz_yyz[k] = -g_xy_0_yyz_yyz[k] * ab_z + g_xy_0_yyz_yyzz[k];

                g_xy_0_yyzz_yzz[k] = -g_xy_0_yyz_yzz[k] * ab_z + g_xy_0_yyz_yzzz[k];

                g_xy_0_yyzz_zzz[k] = -g_xy_0_yyz_zzz[k] * ab_z + g_xy_0_yyz_zzzz[k];
            }

            /// Set up 280-290 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 280 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 281 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 282 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 283 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 284 * ccomps * dcomps);

            auto g_xy_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 285 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 286 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 287 * ccomps * dcomps);

            auto g_xy_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 288 * ccomps * dcomps);

            auto g_xy_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 289 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yzz_xxx, g_xy_0_yzz_xxxz, g_xy_0_yzz_xxy, g_xy_0_yzz_xxyz, g_xy_0_yzz_xxz, g_xy_0_yzz_xxzz, g_xy_0_yzz_xyy, g_xy_0_yzz_xyyz, g_xy_0_yzz_xyz, g_xy_0_yzz_xyzz, g_xy_0_yzz_xzz, g_xy_0_yzz_xzzz, g_xy_0_yzz_yyy, g_xy_0_yzz_yyyz, g_xy_0_yzz_yyz, g_xy_0_yzz_yyzz, g_xy_0_yzz_yzz, g_xy_0_yzz_yzzz, g_xy_0_yzz_zzz, g_xy_0_yzz_zzzz, g_xy_0_yzzz_xxx, g_xy_0_yzzz_xxy, g_xy_0_yzzz_xxz, g_xy_0_yzzz_xyy, g_xy_0_yzzz_xyz, g_xy_0_yzzz_xzz, g_xy_0_yzzz_yyy, g_xy_0_yzzz_yyz, g_xy_0_yzzz_yzz, g_xy_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzzz_xxx[k] = -g_xy_0_yzz_xxx[k] * ab_z + g_xy_0_yzz_xxxz[k];

                g_xy_0_yzzz_xxy[k] = -g_xy_0_yzz_xxy[k] * ab_z + g_xy_0_yzz_xxyz[k];

                g_xy_0_yzzz_xxz[k] = -g_xy_0_yzz_xxz[k] * ab_z + g_xy_0_yzz_xxzz[k];

                g_xy_0_yzzz_xyy[k] = -g_xy_0_yzz_xyy[k] * ab_z + g_xy_0_yzz_xyyz[k];

                g_xy_0_yzzz_xyz[k] = -g_xy_0_yzz_xyz[k] * ab_z + g_xy_0_yzz_xyzz[k];

                g_xy_0_yzzz_xzz[k] = -g_xy_0_yzz_xzz[k] * ab_z + g_xy_0_yzz_xzzz[k];

                g_xy_0_yzzz_yyy[k] = -g_xy_0_yzz_yyy[k] * ab_z + g_xy_0_yzz_yyyz[k];

                g_xy_0_yzzz_yyz[k] = -g_xy_0_yzz_yyz[k] * ab_z + g_xy_0_yzz_yyzz[k];

                g_xy_0_yzzz_yzz[k] = -g_xy_0_yzz_yzz[k] * ab_z + g_xy_0_yzz_yzzz[k];

                g_xy_0_yzzz_zzz[k] = -g_xy_0_yzz_zzz[k] * ab_z + g_xy_0_yzz_zzzz[k];
            }

            /// Set up 290-300 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 290 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 291 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 292 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 293 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 294 * ccomps * dcomps);

            auto g_xy_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 295 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 296 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 297 * ccomps * dcomps);

            auto g_xy_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 298 * ccomps * dcomps);

            auto g_xy_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zzz_xxx, g_xy_0_zzz_xxxz, g_xy_0_zzz_xxy, g_xy_0_zzz_xxyz, g_xy_0_zzz_xxz, g_xy_0_zzz_xxzz, g_xy_0_zzz_xyy, g_xy_0_zzz_xyyz, g_xy_0_zzz_xyz, g_xy_0_zzz_xyzz, g_xy_0_zzz_xzz, g_xy_0_zzz_xzzz, g_xy_0_zzz_yyy, g_xy_0_zzz_yyyz, g_xy_0_zzz_yyz, g_xy_0_zzz_yyzz, g_xy_0_zzz_yzz, g_xy_0_zzz_yzzz, g_xy_0_zzz_zzz, g_xy_0_zzz_zzzz, g_xy_0_zzzz_xxx, g_xy_0_zzzz_xxy, g_xy_0_zzzz_xxz, g_xy_0_zzzz_xyy, g_xy_0_zzzz_xyz, g_xy_0_zzzz_xzz, g_xy_0_zzzz_yyy, g_xy_0_zzzz_yyz, g_xy_0_zzzz_yzz, g_xy_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzzz_xxx[k] = -g_xy_0_zzz_xxx[k] * ab_z + g_xy_0_zzz_xxxz[k];

                g_xy_0_zzzz_xxy[k] = -g_xy_0_zzz_xxy[k] * ab_z + g_xy_0_zzz_xxyz[k];

                g_xy_0_zzzz_xxz[k] = -g_xy_0_zzz_xxz[k] * ab_z + g_xy_0_zzz_xxzz[k];

                g_xy_0_zzzz_xyy[k] = -g_xy_0_zzz_xyy[k] * ab_z + g_xy_0_zzz_xyyz[k];

                g_xy_0_zzzz_xyz[k] = -g_xy_0_zzz_xyz[k] * ab_z + g_xy_0_zzz_xyzz[k];

                g_xy_0_zzzz_xzz[k] = -g_xy_0_zzz_xzz[k] * ab_z + g_xy_0_zzz_xzzz[k];

                g_xy_0_zzzz_yyy[k] = -g_xy_0_zzz_yyy[k] * ab_z + g_xy_0_zzz_yyyz[k];

                g_xy_0_zzzz_yyz[k] = -g_xy_0_zzz_yyz[k] * ab_z + g_xy_0_zzz_yyzz[k];

                g_xy_0_zzzz_yzz[k] = -g_xy_0_zzz_yzz[k] * ab_z + g_xy_0_zzz_yzzz[k];

                g_xy_0_zzzz_zzz[k] = -g_xy_0_zzz_zzz[k] * ab_z + g_xy_0_zzz_zzzz[k];
            }

            /// Set up 300-310 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 305 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 309 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxx_xxx, g_xz_0_xxx_xxxx, g_xz_0_xxx_xxxy, g_xz_0_xxx_xxxz, g_xz_0_xxx_xxy, g_xz_0_xxx_xxyy, g_xz_0_xxx_xxyz, g_xz_0_xxx_xxz, g_xz_0_xxx_xxzz, g_xz_0_xxx_xyy, g_xz_0_xxx_xyyy, g_xz_0_xxx_xyyz, g_xz_0_xxx_xyz, g_xz_0_xxx_xyzz, g_xz_0_xxx_xzz, g_xz_0_xxx_xzzz, g_xz_0_xxx_yyy, g_xz_0_xxx_yyz, g_xz_0_xxx_yzz, g_xz_0_xxx_zzz, g_xz_0_xxxx_xxx, g_xz_0_xxxx_xxy, g_xz_0_xxxx_xxz, g_xz_0_xxxx_xyy, g_xz_0_xxxx_xyz, g_xz_0_xxxx_xzz, g_xz_0_xxxx_yyy, g_xz_0_xxxx_yyz, g_xz_0_xxxx_yzz, g_xz_0_xxxx_zzz, g_z_0_xxx_xxx, g_z_0_xxx_xxy, g_z_0_xxx_xxz, g_z_0_xxx_xyy, g_z_0_xxx_xyz, g_z_0_xxx_xzz, g_z_0_xxx_yyy, g_z_0_xxx_yyz, g_z_0_xxx_yzz, g_z_0_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxx_xxx[k] = -g_z_0_xxx_xxx[k] - g_xz_0_xxx_xxx[k] * ab_x + g_xz_0_xxx_xxxx[k];

                g_xz_0_xxxx_xxy[k] = -g_z_0_xxx_xxy[k] - g_xz_0_xxx_xxy[k] * ab_x + g_xz_0_xxx_xxxy[k];

                g_xz_0_xxxx_xxz[k] = -g_z_0_xxx_xxz[k] - g_xz_0_xxx_xxz[k] * ab_x + g_xz_0_xxx_xxxz[k];

                g_xz_0_xxxx_xyy[k] = -g_z_0_xxx_xyy[k] - g_xz_0_xxx_xyy[k] * ab_x + g_xz_0_xxx_xxyy[k];

                g_xz_0_xxxx_xyz[k] = -g_z_0_xxx_xyz[k] - g_xz_0_xxx_xyz[k] * ab_x + g_xz_0_xxx_xxyz[k];

                g_xz_0_xxxx_xzz[k] = -g_z_0_xxx_xzz[k] - g_xz_0_xxx_xzz[k] * ab_x + g_xz_0_xxx_xxzz[k];

                g_xz_0_xxxx_yyy[k] = -g_z_0_xxx_yyy[k] - g_xz_0_xxx_yyy[k] * ab_x + g_xz_0_xxx_xyyy[k];

                g_xz_0_xxxx_yyz[k] = -g_z_0_xxx_yyz[k] - g_xz_0_xxx_yyz[k] * ab_x + g_xz_0_xxx_xyyz[k];

                g_xz_0_xxxx_yzz[k] = -g_z_0_xxx_yzz[k] - g_xz_0_xxx_yzz[k] * ab_x + g_xz_0_xxx_xyzz[k];

                g_xz_0_xxxx_zzz[k] = -g_z_0_xxx_zzz[k] - g_xz_0_xxx_zzz[k] * ab_x + g_xz_0_xxx_xzzz[k];
            }

            /// Set up 310-320 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 311 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 314 * ccomps * dcomps);

            auto g_xz_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 317 * ccomps * dcomps);

            auto g_xz_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 319 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxx_xxx, g_xz_0_xxx_xxxy, g_xz_0_xxx_xxy, g_xz_0_xxx_xxyy, g_xz_0_xxx_xxyz, g_xz_0_xxx_xxz, g_xz_0_xxx_xyy, g_xz_0_xxx_xyyy, g_xz_0_xxx_xyyz, g_xz_0_xxx_xyz, g_xz_0_xxx_xyzz, g_xz_0_xxx_xzz, g_xz_0_xxx_yyy, g_xz_0_xxx_yyyy, g_xz_0_xxx_yyyz, g_xz_0_xxx_yyz, g_xz_0_xxx_yyzz, g_xz_0_xxx_yzz, g_xz_0_xxx_yzzz, g_xz_0_xxx_zzz, g_xz_0_xxxy_xxx, g_xz_0_xxxy_xxy, g_xz_0_xxxy_xxz, g_xz_0_xxxy_xyy, g_xz_0_xxxy_xyz, g_xz_0_xxxy_xzz, g_xz_0_xxxy_yyy, g_xz_0_xxxy_yyz, g_xz_0_xxxy_yzz, g_xz_0_xxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxy_xxx[k] = -g_xz_0_xxx_xxx[k] * ab_y + g_xz_0_xxx_xxxy[k];

                g_xz_0_xxxy_xxy[k] = -g_xz_0_xxx_xxy[k] * ab_y + g_xz_0_xxx_xxyy[k];

                g_xz_0_xxxy_xxz[k] = -g_xz_0_xxx_xxz[k] * ab_y + g_xz_0_xxx_xxyz[k];

                g_xz_0_xxxy_xyy[k] = -g_xz_0_xxx_xyy[k] * ab_y + g_xz_0_xxx_xyyy[k];

                g_xz_0_xxxy_xyz[k] = -g_xz_0_xxx_xyz[k] * ab_y + g_xz_0_xxx_xyyz[k];

                g_xz_0_xxxy_xzz[k] = -g_xz_0_xxx_xzz[k] * ab_y + g_xz_0_xxx_xyzz[k];

                g_xz_0_xxxy_yyy[k] = -g_xz_0_xxx_yyy[k] * ab_y + g_xz_0_xxx_yyyy[k];

                g_xz_0_xxxy_yyz[k] = -g_xz_0_xxx_yyz[k] * ab_y + g_xz_0_xxx_yyyz[k];

                g_xz_0_xxxy_yzz[k] = -g_xz_0_xxx_yzz[k] * ab_y + g_xz_0_xxx_yyzz[k];

                g_xz_0_xxxy_zzz[k] = -g_xz_0_xxx_zzz[k] * ab_y + g_xz_0_xxx_yzzz[k];
            }

            /// Set up 320-330 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 323 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 324 * ccomps * dcomps);

            auto g_xz_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 325 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 326 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 327 * ccomps * dcomps);

            auto g_xz_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 328 * ccomps * dcomps);

            auto g_xz_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxz_xxx, g_xz_0_xxxz_xxy, g_xz_0_xxxz_xxz, g_xz_0_xxxz_xyy, g_xz_0_xxxz_xyz, g_xz_0_xxxz_xzz, g_xz_0_xxxz_yyy, g_xz_0_xxxz_yyz, g_xz_0_xxxz_yzz, g_xz_0_xxxz_zzz, g_xz_0_xxz_xxx, g_xz_0_xxz_xxxx, g_xz_0_xxz_xxxy, g_xz_0_xxz_xxxz, g_xz_0_xxz_xxy, g_xz_0_xxz_xxyy, g_xz_0_xxz_xxyz, g_xz_0_xxz_xxz, g_xz_0_xxz_xxzz, g_xz_0_xxz_xyy, g_xz_0_xxz_xyyy, g_xz_0_xxz_xyyz, g_xz_0_xxz_xyz, g_xz_0_xxz_xyzz, g_xz_0_xxz_xzz, g_xz_0_xxz_xzzz, g_xz_0_xxz_yyy, g_xz_0_xxz_yyz, g_xz_0_xxz_yzz, g_xz_0_xxz_zzz, g_z_0_xxz_xxx, g_z_0_xxz_xxy, g_z_0_xxz_xxz, g_z_0_xxz_xyy, g_z_0_xxz_xyz, g_z_0_xxz_xzz, g_z_0_xxz_yyy, g_z_0_xxz_yyz, g_z_0_xxz_yzz, g_z_0_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxz_xxx[k] = -g_z_0_xxz_xxx[k] - g_xz_0_xxz_xxx[k] * ab_x + g_xz_0_xxz_xxxx[k];

                g_xz_0_xxxz_xxy[k] = -g_z_0_xxz_xxy[k] - g_xz_0_xxz_xxy[k] * ab_x + g_xz_0_xxz_xxxy[k];

                g_xz_0_xxxz_xxz[k] = -g_z_0_xxz_xxz[k] - g_xz_0_xxz_xxz[k] * ab_x + g_xz_0_xxz_xxxz[k];

                g_xz_0_xxxz_xyy[k] = -g_z_0_xxz_xyy[k] - g_xz_0_xxz_xyy[k] * ab_x + g_xz_0_xxz_xxyy[k];

                g_xz_0_xxxz_xyz[k] = -g_z_0_xxz_xyz[k] - g_xz_0_xxz_xyz[k] * ab_x + g_xz_0_xxz_xxyz[k];

                g_xz_0_xxxz_xzz[k] = -g_z_0_xxz_xzz[k] - g_xz_0_xxz_xzz[k] * ab_x + g_xz_0_xxz_xxzz[k];

                g_xz_0_xxxz_yyy[k] = -g_z_0_xxz_yyy[k] - g_xz_0_xxz_yyy[k] * ab_x + g_xz_0_xxz_xyyy[k];

                g_xz_0_xxxz_yyz[k] = -g_z_0_xxz_yyz[k] - g_xz_0_xxz_yyz[k] * ab_x + g_xz_0_xxz_xyyz[k];

                g_xz_0_xxxz_yzz[k] = -g_z_0_xxz_yzz[k] - g_xz_0_xxz_yzz[k] * ab_x + g_xz_0_xxz_xyzz[k];

                g_xz_0_xxxz_zzz[k] = -g_z_0_xxz_zzz[k] - g_xz_0_xxz_zzz[k] * ab_x + g_xz_0_xxz_xzzz[k];
            }

            /// Set up 330-340 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 330 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 331 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 332 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 333 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 334 * ccomps * dcomps);

            auto g_xz_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 335 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 339 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxy_xxx, g_xz_0_xxy_xxxy, g_xz_0_xxy_xxy, g_xz_0_xxy_xxyy, g_xz_0_xxy_xxyz, g_xz_0_xxy_xxz, g_xz_0_xxy_xyy, g_xz_0_xxy_xyyy, g_xz_0_xxy_xyyz, g_xz_0_xxy_xyz, g_xz_0_xxy_xyzz, g_xz_0_xxy_xzz, g_xz_0_xxy_yyy, g_xz_0_xxy_yyyy, g_xz_0_xxy_yyyz, g_xz_0_xxy_yyz, g_xz_0_xxy_yyzz, g_xz_0_xxy_yzz, g_xz_0_xxy_yzzz, g_xz_0_xxy_zzz, g_xz_0_xxyy_xxx, g_xz_0_xxyy_xxy, g_xz_0_xxyy_xxz, g_xz_0_xxyy_xyy, g_xz_0_xxyy_xyz, g_xz_0_xxyy_xzz, g_xz_0_xxyy_yyy, g_xz_0_xxyy_yyz, g_xz_0_xxyy_yzz, g_xz_0_xxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyy_xxx[k] = -g_xz_0_xxy_xxx[k] * ab_y + g_xz_0_xxy_xxxy[k];

                g_xz_0_xxyy_xxy[k] = -g_xz_0_xxy_xxy[k] * ab_y + g_xz_0_xxy_xxyy[k];

                g_xz_0_xxyy_xxz[k] = -g_xz_0_xxy_xxz[k] * ab_y + g_xz_0_xxy_xxyz[k];

                g_xz_0_xxyy_xyy[k] = -g_xz_0_xxy_xyy[k] * ab_y + g_xz_0_xxy_xyyy[k];

                g_xz_0_xxyy_xyz[k] = -g_xz_0_xxy_xyz[k] * ab_y + g_xz_0_xxy_xyyz[k];

                g_xz_0_xxyy_xzz[k] = -g_xz_0_xxy_xzz[k] * ab_y + g_xz_0_xxy_xyzz[k];

                g_xz_0_xxyy_yyy[k] = -g_xz_0_xxy_yyy[k] * ab_y + g_xz_0_xxy_yyyy[k];

                g_xz_0_xxyy_yyz[k] = -g_xz_0_xxy_yyz[k] * ab_y + g_xz_0_xxy_yyyz[k];

                g_xz_0_xxyy_yzz[k] = -g_xz_0_xxy_yzz[k] * ab_y + g_xz_0_xxy_yyzz[k];

                g_xz_0_xxyy_zzz[k] = -g_xz_0_xxy_zzz[k] * ab_y + g_xz_0_xxy_yzzz[k];
            }

            /// Set up 340-350 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 341 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 347 * ccomps * dcomps);

            auto g_xz_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 349 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyz_xxx, g_xz_0_xxyz_xxy, g_xz_0_xxyz_xxz, g_xz_0_xxyz_xyy, g_xz_0_xxyz_xyz, g_xz_0_xxyz_xzz, g_xz_0_xxyz_yyy, g_xz_0_xxyz_yyz, g_xz_0_xxyz_yzz, g_xz_0_xxyz_zzz, g_xz_0_xxz_xxx, g_xz_0_xxz_xxxy, g_xz_0_xxz_xxy, g_xz_0_xxz_xxyy, g_xz_0_xxz_xxyz, g_xz_0_xxz_xxz, g_xz_0_xxz_xyy, g_xz_0_xxz_xyyy, g_xz_0_xxz_xyyz, g_xz_0_xxz_xyz, g_xz_0_xxz_xyzz, g_xz_0_xxz_xzz, g_xz_0_xxz_yyy, g_xz_0_xxz_yyyy, g_xz_0_xxz_yyyz, g_xz_0_xxz_yyz, g_xz_0_xxz_yyzz, g_xz_0_xxz_yzz, g_xz_0_xxz_yzzz, g_xz_0_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyz_xxx[k] = -g_xz_0_xxz_xxx[k] * ab_y + g_xz_0_xxz_xxxy[k];

                g_xz_0_xxyz_xxy[k] = -g_xz_0_xxz_xxy[k] * ab_y + g_xz_0_xxz_xxyy[k];

                g_xz_0_xxyz_xxz[k] = -g_xz_0_xxz_xxz[k] * ab_y + g_xz_0_xxz_xxyz[k];

                g_xz_0_xxyz_xyy[k] = -g_xz_0_xxz_xyy[k] * ab_y + g_xz_0_xxz_xyyy[k];

                g_xz_0_xxyz_xyz[k] = -g_xz_0_xxz_xyz[k] * ab_y + g_xz_0_xxz_xyyz[k];

                g_xz_0_xxyz_xzz[k] = -g_xz_0_xxz_xzz[k] * ab_y + g_xz_0_xxz_xyzz[k];

                g_xz_0_xxyz_yyy[k] = -g_xz_0_xxz_yyy[k] * ab_y + g_xz_0_xxz_yyyy[k];

                g_xz_0_xxyz_yyz[k] = -g_xz_0_xxz_yyz[k] * ab_y + g_xz_0_xxz_yyyz[k];

                g_xz_0_xxyz_yzz[k] = -g_xz_0_xxz_yzz[k] * ab_y + g_xz_0_xxz_yyzz[k];

                g_xz_0_xxyz_zzz[k] = -g_xz_0_xxz_zzz[k] * ab_y + g_xz_0_xxz_yzzz[k];
            }

            /// Set up 350-360 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 353 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxzz_xxx, g_xz_0_xxzz_xxy, g_xz_0_xxzz_xxz, g_xz_0_xxzz_xyy, g_xz_0_xxzz_xyz, g_xz_0_xxzz_xzz, g_xz_0_xxzz_yyy, g_xz_0_xxzz_yyz, g_xz_0_xxzz_yzz, g_xz_0_xxzz_zzz, g_xz_0_xzz_xxx, g_xz_0_xzz_xxxx, g_xz_0_xzz_xxxy, g_xz_0_xzz_xxxz, g_xz_0_xzz_xxy, g_xz_0_xzz_xxyy, g_xz_0_xzz_xxyz, g_xz_0_xzz_xxz, g_xz_0_xzz_xxzz, g_xz_0_xzz_xyy, g_xz_0_xzz_xyyy, g_xz_0_xzz_xyyz, g_xz_0_xzz_xyz, g_xz_0_xzz_xyzz, g_xz_0_xzz_xzz, g_xz_0_xzz_xzzz, g_xz_0_xzz_yyy, g_xz_0_xzz_yyz, g_xz_0_xzz_yzz, g_xz_0_xzz_zzz, g_z_0_xzz_xxx, g_z_0_xzz_xxy, g_z_0_xzz_xxz, g_z_0_xzz_xyy, g_z_0_xzz_xyz, g_z_0_xzz_xzz, g_z_0_xzz_yyy, g_z_0_xzz_yyz, g_z_0_xzz_yzz, g_z_0_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxzz_xxx[k] = -g_z_0_xzz_xxx[k] - g_xz_0_xzz_xxx[k] * ab_x + g_xz_0_xzz_xxxx[k];

                g_xz_0_xxzz_xxy[k] = -g_z_0_xzz_xxy[k] - g_xz_0_xzz_xxy[k] * ab_x + g_xz_0_xzz_xxxy[k];

                g_xz_0_xxzz_xxz[k] = -g_z_0_xzz_xxz[k] - g_xz_0_xzz_xxz[k] * ab_x + g_xz_0_xzz_xxxz[k];

                g_xz_0_xxzz_xyy[k] = -g_z_0_xzz_xyy[k] - g_xz_0_xzz_xyy[k] * ab_x + g_xz_0_xzz_xxyy[k];

                g_xz_0_xxzz_xyz[k] = -g_z_0_xzz_xyz[k] - g_xz_0_xzz_xyz[k] * ab_x + g_xz_0_xzz_xxyz[k];

                g_xz_0_xxzz_xzz[k] = -g_z_0_xzz_xzz[k] - g_xz_0_xzz_xzz[k] * ab_x + g_xz_0_xzz_xxzz[k];

                g_xz_0_xxzz_yyy[k] = -g_z_0_xzz_yyy[k] - g_xz_0_xzz_yyy[k] * ab_x + g_xz_0_xzz_xyyy[k];

                g_xz_0_xxzz_yyz[k] = -g_z_0_xzz_yyz[k] - g_xz_0_xzz_yyz[k] * ab_x + g_xz_0_xzz_xyyz[k];

                g_xz_0_xxzz_yzz[k] = -g_z_0_xzz_yzz[k] - g_xz_0_xzz_yzz[k] * ab_x + g_xz_0_xzz_xyzz[k];

                g_xz_0_xxzz_zzz[k] = -g_z_0_xzz_zzz[k] - g_xz_0_xzz_zzz[k] * ab_x + g_xz_0_xzz_xzzz[k];
            }

            /// Set up 360-370 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 365 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 369 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyy_xxx, g_xz_0_xyy_xxxy, g_xz_0_xyy_xxy, g_xz_0_xyy_xxyy, g_xz_0_xyy_xxyz, g_xz_0_xyy_xxz, g_xz_0_xyy_xyy, g_xz_0_xyy_xyyy, g_xz_0_xyy_xyyz, g_xz_0_xyy_xyz, g_xz_0_xyy_xyzz, g_xz_0_xyy_xzz, g_xz_0_xyy_yyy, g_xz_0_xyy_yyyy, g_xz_0_xyy_yyyz, g_xz_0_xyy_yyz, g_xz_0_xyy_yyzz, g_xz_0_xyy_yzz, g_xz_0_xyy_yzzz, g_xz_0_xyy_zzz, g_xz_0_xyyy_xxx, g_xz_0_xyyy_xxy, g_xz_0_xyyy_xxz, g_xz_0_xyyy_xyy, g_xz_0_xyyy_xyz, g_xz_0_xyyy_xzz, g_xz_0_xyyy_yyy, g_xz_0_xyyy_yyz, g_xz_0_xyyy_yzz, g_xz_0_xyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyy_xxx[k] = -g_xz_0_xyy_xxx[k] * ab_y + g_xz_0_xyy_xxxy[k];

                g_xz_0_xyyy_xxy[k] = -g_xz_0_xyy_xxy[k] * ab_y + g_xz_0_xyy_xxyy[k];

                g_xz_0_xyyy_xxz[k] = -g_xz_0_xyy_xxz[k] * ab_y + g_xz_0_xyy_xxyz[k];

                g_xz_0_xyyy_xyy[k] = -g_xz_0_xyy_xyy[k] * ab_y + g_xz_0_xyy_xyyy[k];

                g_xz_0_xyyy_xyz[k] = -g_xz_0_xyy_xyz[k] * ab_y + g_xz_0_xyy_xyyz[k];

                g_xz_0_xyyy_xzz[k] = -g_xz_0_xyy_xzz[k] * ab_y + g_xz_0_xyy_xyzz[k];

                g_xz_0_xyyy_yyy[k] = -g_xz_0_xyy_yyy[k] * ab_y + g_xz_0_xyy_yyyy[k];

                g_xz_0_xyyy_yyz[k] = -g_xz_0_xyy_yyz[k] * ab_y + g_xz_0_xyy_yyyz[k];

                g_xz_0_xyyy_yzz[k] = -g_xz_0_xyy_yzz[k] * ab_y + g_xz_0_xyy_yyzz[k];

                g_xz_0_xyyy_zzz[k] = -g_xz_0_xyy_zzz[k] * ab_y + g_xz_0_xyy_yzzz[k];
            }

            /// Set up 370-380 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 371 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 377 * ccomps * dcomps);

            auto g_xz_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 378 * ccomps * dcomps);

            auto g_xz_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 379 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyz_xxx, g_xz_0_xyyz_xxy, g_xz_0_xyyz_xxz, g_xz_0_xyyz_xyy, g_xz_0_xyyz_xyz, g_xz_0_xyyz_xzz, g_xz_0_xyyz_yyy, g_xz_0_xyyz_yyz, g_xz_0_xyyz_yzz, g_xz_0_xyyz_zzz, g_xz_0_xyz_xxx, g_xz_0_xyz_xxxy, g_xz_0_xyz_xxy, g_xz_0_xyz_xxyy, g_xz_0_xyz_xxyz, g_xz_0_xyz_xxz, g_xz_0_xyz_xyy, g_xz_0_xyz_xyyy, g_xz_0_xyz_xyyz, g_xz_0_xyz_xyz, g_xz_0_xyz_xyzz, g_xz_0_xyz_xzz, g_xz_0_xyz_yyy, g_xz_0_xyz_yyyy, g_xz_0_xyz_yyyz, g_xz_0_xyz_yyz, g_xz_0_xyz_yyzz, g_xz_0_xyz_yzz, g_xz_0_xyz_yzzz, g_xz_0_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyz_xxx[k] = -g_xz_0_xyz_xxx[k] * ab_y + g_xz_0_xyz_xxxy[k];

                g_xz_0_xyyz_xxy[k] = -g_xz_0_xyz_xxy[k] * ab_y + g_xz_0_xyz_xxyy[k];

                g_xz_0_xyyz_xxz[k] = -g_xz_0_xyz_xxz[k] * ab_y + g_xz_0_xyz_xxyz[k];

                g_xz_0_xyyz_xyy[k] = -g_xz_0_xyz_xyy[k] * ab_y + g_xz_0_xyz_xyyy[k];

                g_xz_0_xyyz_xyz[k] = -g_xz_0_xyz_xyz[k] * ab_y + g_xz_0_xyz_xyyz[k];

                g_xz_0_xyyz_xzz[k] = -g_xz_0_xyz_xzz[k] * ab_y + g_xz_0_xyz_xyzz[k];

                g_xz_0_xyyz_yyy[k] = -g_xz_0_xyz_yyy[k] * ab_y + g_xz_0_xyz_yyyy[k];

                g_xz_0_xyyz_yyz[k] = -g_xz_0_xyz_yyz[k] * ab_y + g_xz_0_xyz_yyyz[k];

                g_xz_0_xyyz_yzz[k] = -g_xz_0_xyz_yzz[k] * ab_y + g_xz_0_xyz_yyzz[k];

                g_xz_0_xyyz_zzz[k] = -g_xz_0_xyz_zzz[k] * ab_y + g_xz_0_xyz_yzzz[k];
            }

            /// Set up 380-390 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 380 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 381 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 382 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 383 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 384 * ccomps * dcomps);

            auto g_xz_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 385 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 386 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 387 * ccomps * dcomps);

            auto g_xz_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 388 * ccomps * dcomps);

            auto g_xz_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyzz_xxx, g_xz_0_xyzz_xxy, g_xz_0_xyzz_xxz, g_xz_0_xyzz_xyy, g_xz_0_xyzz_xyz, g_xz_0_xyzz_xzz, g_xz_0_xyzz_yyy, g_xz_0_xyzz_yyz, g_xz_0_xyzz_yzz, g_xz_0_xyzz_zzz, g_xz_0_xzz_xxx, g_xz_0_xzz_xxxy, g_xz_0_xzz_xxy, g_xz_0_xzz_xxyy, g_xz_0_xzz_xxyz, g_xz_0_xzz_xxz, g_xz_0_xzz_xyy, g_xz_0_xzz_xyyy, g_xz_0_xzz_xyyz, g_xz_0_xzz_xyz, g_xz_0_xzz_xyzz, g_xz_0_xzz_xzz, g_xz_0_xzz_yyy, g_xz_0_xzz_yyyy, g_xz_0_xzz_yyyz, g_xz_0_xzz_yyz, g_xz_0_xzz_yyzz, g_xz_0_xzz_yzz, g_xz_0_xzz_yzzz, g_xz_0_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyzz_xxx[k] = -g_xz_0_xzz_xxx[k] * ab_y + g_xz_0_xzz_xxxy[k];

                g_xz_0_xyzz_xxy[k] = -g_xz_0_xzz_xxy[k] * ab_y + g_xz_0_xzz_xxyy[k];

                g_xz_0_xyzz_xxz[k] = -g_xz_0_xzz_xxz[k] * ab_y + g_xz_0_xzz_xxyz[k];

                g_xz_0_xyzz_xyy[k] = -g_xz_0_xzz_xyy[k] * ab_y + g_xz_0_xzz_xyyy[k];

                g_xz_0_xyzz_xyz[k] = -g_xz_0_xzz_xyz[k] * ab_y + g_xz_0_xzz_xyyz[k];

                g_xz_0_xyzz_xzz[k] = -g_xz_0_xzz_xzz[k] * ab_y + g_xz_0_xzz_xyzz[k];

                g_xz_0_xyzz_yyy[k] = -g_xz_0_xzz_yyy[k] * ab_y + g_xz_0_xzz_yyyy[k];

                g_xz_0_xyzz_yyz[k] = -g_xz_0_xzz_yyz[k] * ab_y + g_xz_0_xzz_yyyz[k];

                g_xz_0_xyzz_yzz[k] = -g_xz_0_xzz_yzz[k] * ab_y + g_xz_0_xzz_yyzz[k];

                g_xz_0_xyzz_zzz[k] = -g_xz_0_xzz_zzz[k] * ab_y + g_xz_0_xzz_yzzz[k];
            }

            /// Set up 390-400 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 390 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 391 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 392 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 393 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 394 * ccomps * dcomps);

            auto g_xz_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 395 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 396 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 397 * ccomps * dcomps);

            auto g_xz_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 398 * ccomps * dcomps);

            auto g_xz_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 399 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzzz_xxx, g_xz_0_xzzz_xxy, g_xz_0_xzzz_xxz, g_xz_0_xzzz_xyy, g_xz_0_xzzz_xyz, g_xz_0_xzzz_xzz, g_xz_0_xzzz_yyy, g_xz_0_xzzz_yyz, g_xz_0_xzzz_yzz, g_xz_0_xzzz_zzz, g_xz_0_zzz_xxx, g_xz_0_zzz_xxxx, g_xz_0_zzz_xxxy, g_xz_0_zzz_xxxz, g_xz_0_zzz_xxy, g_xz_0_zzz_xxyy, g_xz_0_zzz_xxyz, g_xz_0_zzz_xxz, g_xz_0_zzz_xxzz, g_xz_0_zzz_xyy, g_xz_0_zzz_xyyy, g_xz_0_zzz_xyyz, g_xz_0_zzz_xyz, g_xz_0_zzz_xyzz, g_xz_0_zzz_xzz, g_xz_0_zzz_xzzz, g_xz_0_zzz_yyy, g_xz_0_zzz_yyz, g_xz_0_zzz_yzz, g_xz_0_zzz_zzz, g_z_0_zzz_xxx, g_z_0_zzz_xxy, g_z_0_zzz_xxz, g_z_0_zzz_xyy, g_z_0_zzz_xyz, g_z_0_zzz_xzz, g_z_0_zzz_yyy, g_z_0_zzz_yyz, g_z_0_zzz_yzz, g_z_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzzz_xxx[k] = -g_z_0_zzz_xxx[k] - g_xz_0_zzz_xxx[k] * ab_x + g_xz_0_zzz_xxxx[k];

                g_xz_0_xzzz_xxy[k] = -g_z_0_zzz_xxy[k] - g_xz_0_zzz_xxy[k] * ab_x + g_xz_0_zzz_xxxy[k];

                g_xz_0_xzzz_xxz[k] = -g_z_0_zzz_xxz[k] - g_xz_0_zzz_xxz[k] * ab_x + g_xz_0_zzz_xxxz[k];

                g_xz_0_xzzz_xyy[k] = -g_z_0_zzz_xyy[k] - g_xz_0_zzz_xyy[k] * ab_x + g_xz_0_zzz_xxyy[k];

                g_xz_0_xzzz_xyz[k] = -g_z_0_zzz_xyz[k] - g_xz_0_zzz_xyz[k] * ab_x + g_xz_0_zzz_xxyz[k];

                g_xz_0_xzzz_xzz[k] = -g_z_0_zzz_xzz[k] - g_xz_0_zzz_xzz[k] * ab_x + g_xz_0_zzz_xxzz[k];

                g_xz_0_xzzz_yyy[k] = -g_z_0_zzz_yyy[k] - g_xz_0_zzz_yyy[k] * ab_x + g_xz_0_zzz_xyyy[k];

                g_xz_0_xzzz_yyz[k] = -g_z_0_zzz_yyz[k] - g_xz_0_zzz_yyz[k] * ab_x + g_xz_0_zzz_xyyz[k];

                g_xz_0_xzzz_yzz[k] = -g_z_0_zzz_yzz[k] - g_xz_0_zzz_yzz[k] * ab_x + g_xz_0_zzz_xyzz[k];

                g_xz_0_xzzz_zzz[k] = -g_z_0_zzz_zzz[k] - g_xz_0_zzz_zzz[k] * ab_x + g_xz_0_zzz_xzzz[k];
            }

            /// Set up 400-410 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 400 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 401 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 402 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 403 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 404 * ccomps * dcomps);

            auto g_xz_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 405 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 406 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 407 * ccomps * dcomps);

            auto g_xz_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 408 * ccomps * dcomps);

            auto g_xz_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 409 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyy_xxx, g_xz_0_yyy_xxxy, g_xz_0_yyy_xxy, g_xz_0_yyy_xxyy, g_xz_0_yyy_xxyz, g_xz_0_yyy_xxz, g_xz_0_yyy_xyy, g_xz_0_yyy_xyyy, g_xz_0_yyy_xyyz, g_xz_0_yyy_xyz, g_xz_0_yyy_xyzz, g_xz_0_yyy_xzz, g_xz_0_yyy_yyy, g_xz_0_yyy_yyyy, g_xz_0_yyy_yyyz, g_xz_0_yyy_yyz, g_xz_0_yyy_yyzz, g_xz_0_yyy_yzz, g_xz_0_yyy_yzzz, g_xz_0_yyy_zzz, g_xz_0_yyyy_xxx, g_xz_0_yyyy_xxy, g_xz_0_yyyy_xxz, g_xz_0_yyyy_xyy, g_xz_0_yyyy_xyz, g_xz_0_yyyy_xzz, g_xz_0_yyyy_yyy, g_xz_0_yyyy_yyz, g_xz_0_yyyy_yzz, g_xz_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyy_xxx[k] = -g_xz_0_yyy_xxx[k] * ab_y + g_xz_0_yyy_xxxy[k];

                g_xz_0_yyyy_xxy[k] = -g_xz_0_yyy_xxy[k] * ab_y + g_xz_0_yyy_xxyy[k];

                g_xz_0_yyyy_xxz[k] = -g_xz_0_yyy_xxz[k] * ab_y + g_xz_0_yyy_xxyz[k];

                g_xz_0_yyyy_xyy[k] = -g_xz_0_yyy_xyy[k] * ab_y + g_xz_0_yyy_xyyy[k];

                g_xz_0_yyyy_xyz[k] = -g_xz_0_yyy_xyz[k] * ab_y + g_xz_0_yyy_xyyz[k];

                g_xz_0_yyyy_xzz[k] = -g_xz_0_yyy_xzz[k] * ab_y + g_xz_0_yyy_xyzz[k];

                g_xz_0_yyyy_yyy[k] = -g_xz_0_yyy_yyy[k] * ab_y + g_xz_0_yyy_yyyy[k];

                g_xz_0_yyyy_yyz[k] = -g_xz_0_yyy_yyz[k] * ab_y + g_xz_0_yyy_yyyz[k];

                g_xz_0_yyyy_yzz[k] = -g_xz_0_yyy_yzz[k] * ab_y + g_xz_0_yyy_yyzz[k];

                g_xz_0_yyyy_zzz[k] = -g_xz_0_yyy_zzz[k] * ab_y + g_xz_0_yyy_yzzz[k];
            }

            /// Set up 410-420 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 410 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 411 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 412 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 413 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 414 * ccomps * dcomps);

            auto g_xz_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 415 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 416 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 417 * ccomps * dcomps);

            auto g_xz_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 418 * ccomps * dcomps);

            auto g_xz_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyz_xxx, g_xz_0_yyyz_xxy, g_xz_0_yyyz_xxz, g_xz_0_yyyz_xyy, g_xz_0_yyyz_xyz, g_xz_0_yyyz_xzz, g_xz_0_yyyz_yyy, g_xz_0_yyyz_yyz, g_xz_0_yyyz_yzz, g_xz_0_yyyz_zzz, g_xz_0_yyz_xxx, g_xz_0_yyz_xxxy, g_xz_0_yyz_xxy, g_xz_0_yyz_xxyy, g_xz_0_yyz_xxyz, g_xz_0_yyz_xxz, g_xz_0_yyz_xyy, g_xz_0_yyz_xyyy, g_xz_0_yyz_xyyz, g_xz_0_yyz_xyz, g_xz_0_yyz_xyzz, g_xz_0_yyz_xzz, g_xz_0_yyz_yyy, g_xz_0_yyz_yyyy, g_xz_0_yyz_yyyz, g_xz_0_yyz_yyz, g_xz_0_yyz_yyzz, g_xz_0_yyz_yzz, g_xz_0_yyz_yzzz, g_xz_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyz_xxx[k] = -g_xz_0_yyz_xxx[k] * ab_y + g_xz_0_yyz_xxxy[k];

                g_xz_0_yyyz_xxy[k] = -g_xz_0_yyz_xxy[k] * ab_y + g_xz_0_yyz_xxyy[k];

                g_xz_0_yyyz_xxz[k] = -g_xz_0_yyz_xxz[k] * ab_y + g_xz_0_yyz_xxyz[k];

                g_xz_0_yyyz_xyy[k] = -g_xz_0_yyz_xyy[k] * ab_y + g_xz_0_yyz_xyyy[k];

                g_xz_0_yyyz_xyz[k] = -g_xz_0_yyz_xyz[k] * ab_y + g_xz_0_yyz_xyyz[k];

                g_xz_0_yyyz_xzz[k] = -g_xz_0_yyz_xzz[k] * ab_y + g_xz_0_yyz_xyzz[k];

                g_xz_0_yyyz_yyy[k] = -g_xz_0_yyz_yyy[k] * ab_y + g_xz_0_yyz_yyyy[k];

                g_xz_0_yyyz_yyz[k] = -g_xz_0_yyz_yyz[k] * ab_y + g_xz_0_yyz_yyyz[k];

                g_xz_0_yyyz_yzz[k] = -g_xz_0_yyz_yzz[k] * ab_y + g_xz_0_yyz_yyzz[k];

                g_xz_0_yyyz_zzz[k] = -g_xz_0_yyz_zzz[k] * ab_y + g_xz_0_yyz_yzzz[k];
            }

            /// Set up 420-430 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 420 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 421 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 422 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 423 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 424 * ccomps * dcomps);

            auto g_xz_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 425 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 426 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 427 * ccomps * dcomps);

            auto g_xz_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 428 * ccomps * dcomps);

            auto g_xz_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 429 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyzz_xxx, g_xz_0_yyzz_xxy, g_xz_0_yyzz_xxz, g_xz_0_yyzz_xyy, g_xz_0_yyzz_xyz, g_xz_0_yyzz_xzz, g_xz_0_yyzz_yyy, g_xz_0_yyzz_yyz, g_xz_0_yyzz_yzz, g_xz_0_yyzz_zzz, g_xz_0_yzz_xxx, g_xz_0_yzz_xxxy, g_xz_0_yzz_xxy, g_xz_0_yzz_xxyy, g_xz_0_yzz_xxyz, g_xz_0_yzz_xxz, g_xz_0_yzz_xyy, g_xz_0_yzz_xyyy, g_xz_0_yzz_xyyz, g_xz_0_yzz_xyz, g_xz_0_yzz_xyzz, g_xz_0_yzz_xzz, g_xz_0_yzz_yyy, g_xz_0_yzz_yyyy, g_xz_0_yzz_yyyz, g_xz_0_yzz_yyz, g_xz_0_yzz_yyzz, g_xz_0_yzz_yzz, g_xz_0_yzz_yzzz, g_xz_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyzz_xxx[k] = -g_xz_0_yzz_xxx[k] * ab_y + g_xz_0_yzz_xxxy[k];

                g_xz_0_yyzz_xxy[k] = -g_xz_0_yzz_xxy[k] * ab_y + g_xz_0_yzz_xxyy[k];

                g_xz_0_yyzz_xxz[k] = -g_xz_0_yzz_xxz[k] * ab_y + g_xz_0_yzz_xxyz[k];

                g_xz_0_yyzz_xyy[k] = -g_xz_0_yzz_xyy[k] * ab_y + g_xz_0_yzz_xyyy[k];

                g_xz_0_yyzz_xyz[k] = -g_xz_0_yzz_xyz[k] * ab_y + g_xz_0_yzz_xyyz[k];

                g_xz_0_yyzz_xzz[k] = -g_xz_0_yzz_xzz[k] * ab_y + g_xz_0_yzz_xyzz[k];

                g_xz_0_yyzz_yyy[k] = -g_xz_0_yzz_yyy[k] * ab_y + g_xz_0_yzz_yyyy[k];

                g_xz_0_yyzz_yyz[k] = -g_xz_0_yzz_yyz[k] * ab_y + g_xz_0_yzz_yyyz[k];

                g_xz_0_yyzz_yzz[k] = -g_xz_0_yzz_yzz[k] * ab_y + g_xz_0_yzz_yyzz[k];

                g_xz_0_yyzz_zzz[k] = -g_xz_0_yzz_zzz[k] * ab_y + g_xz_0_yzz_yzzz[k];
            }

            /// Set up 430-440 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 430 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 431 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 432 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 433 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 434 * ccomps * dcomps);

            auto g_xz_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 435 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 436 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 437 * ccomps * dcomps);

            auto g_xz_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 438 * ccomps * dcomps);

            auto g_xz_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 439 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzzz_xxx, g_xz_0_yzzz_xxy, g_xz_0_yzzz_xxz, g_xz_0_yzzz_xyy, g_xz_0_yzzz_xyz, g_xz_0_yzzz_xzz, g_xz_0_yzzz_yyy, g_xz_0_yzzz_yyz, g_xz_0_yzzz_yzz, g_xz_0_yzzz_zzz, g_xz_0_zzz_xxx, g_xz_0_zzz_xxxy, g_xz_0_zzz_xxy, g_xz_0_zzz_xxyy, g_xz_0_zzz_xxyz, g_xz_0_zzz_xxz, g_xz_0_zzz_xyy, g_xz_0_zzz_xyyy, g_xz_0_zzz_xyyz, g_xz_0_zzz_xyz, g_xz_0_zzz_xyzz, g_xz_0_zzz_xzz, g_xz_0_zzz_yyy, g_xz_0_zzz_yyyy, g_xz_0_zzz_yyyz, g_xz_0_zzz_yyz, g_xz_0_zzz_yyzz, g_xz_0_zzz_yzz, g_xz_0_zzz_yzzz, g_xz_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzzz_xxx[k] = -g_xz_0_zzz_xxx[k] * ab_y + g_xz_0_zzz_xxxy[k];

                g_xz_0_yzzz_xxy[k] = -g_xz_0_zzz_xxy[k] * ab_y + g_xz_0_zzz_xxyy[k];

                g_xz_0_yzzz_xxz[k] = -g_xz_0_zzz_xxz[k] * ab_y + g_xz_0_zzz_xxyz[k];

                g_xz_0_yzzz_xyy[k] = -g_xz_0_zzz_xyy[k] * ab_y + g_xz_0_zzz_xyyy[k];

                g_xz_0_yzzz_xyz[k] = -g_xz_0_zzz_xyz[k] * ab_y + g_xz_0_zzz_xyyz[k];

                g_xz_0_yzzz_xzz[k] = -g_xz_0_zzz_xzz[k] * ab_y + g_xz_0_zzz_xyzz[k];

                g_xz_0_yzzz_yyy[k] = -g_xz_0_zzz_yyy[k] * ab_y + g_xz_0_zzz_yyyy[k];

                g_xz_0_yzzz_yyz[k] = -g_xz_0_zzz_yyz[k] * ab_y + g_xz_0_zzz_yyyz[k];

                g_xz_0_yzzz_yzz[k] = -g_xz_0_zzz_yzz[k] * ab_y + g_xz_0_zzz_yyzz[k];

                g_xz_0_yzzz_zzz[k] = -g_xz_0_zzz_zzz[k] * ab_y + g_xz_0_zzz_yzzz[k];
            }

            /// Set up 440-450 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 440 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 441 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 442 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 443 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 444 * ccomps * dcomps);

            auto g_xz_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 445 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 446 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 447 * ccomps * dcomps);

            auto g_xz_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 448 * ccomps * dcomps);

            auto g_xz_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzz_xxx, g_x_0_zzz_xxy, g_x_0_zzz_xxz, g_x_0_zzz_xyy, g_x_0_zzz_xyz, g_x_0_zzz_xzz, g_x_0_zzz_yyy, g_x_0_zzz_yyz, g_x_0_zzz_yzz, g_x_0_zzz_zzz, g_xz_0_zzz_xxx, g_xz_0_zzz_xxxz, g_xz_0_zzz_xxy, g_xz_0_zzz_xxyz, g_xz_0_zzz_xxz, g_xz_0_zzz_xxzz, g_xz_0_zzz_xyy, g_xz_0_zzz_xyyz, g_xz_0_zzz_xyz, g_xz_0_zzz_xyzz, g_xz_0_zzz_xzz, g_xz_0_zzz_xzzz, g_xz_0_zzz_yyy, g_xz_0_zzz_yyyz, g_xz_0_zzz_yyz, g_xz_0_zzz_yyzz, g_xz_0_zzz_yzz, g_xz_0_zzz_yzzz, g_xz_0_zzz_zzz, g_xz_0_zzz_zzzz, g_xz_0_zzzz_xxx, g_xz_0_zzzz_xxy, g_xz_0_zzzz_xxz, g_xz_0_zzzz_xyy, g_xz_0_zzzz_xyz, g_xz_0_zzzz_xzz, g_xz_0_zzzz_yyy, g_xz_0_zzzz_yyz, g_xz_0_zzzz_yzz, g_xz_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzzz_xxx[k] = -g_x_0_zzz_xxx[k] - g_xz_0_zzz_xxx[k] * ab_z + g_xz_0_zzz_xxxz[k];

                g_xz_0_zzzz_xxy[k] = -g_x_0_zzz_xxy[k] - g_xz_0_zzz_xxy[k] * ab_z + g_xz_0_zzz_xxyz[k];

                g_xz_0_zzzz_xxz[k] = -g_x_0_zzz_xxz[k] - g_xz_0_zzz_xxz[k] * ab_z + g_xz_0_zzz_xxzz[k];

                g_xz_0_zzzz_xyy[k] = -g_x_0_zzz_xyy[k] - g_xz_0_zzz_xyy[k] * ab_z + g_xz_0_zzz_xyyz[k];

                g_xz_0_zzzz_xyz[k] = -g_x_0_zzz_xyz[k] - g_xz_0_zzz_xyz[k] * ab_z + g_xz_0_zzz_xyzz[k];

                g_xz_0_zzzz_xzz[k] = -g_x_0_zzz_xzz[k] - g_xz_0_zzz_xzz[k] * ab_z + g_xz_0_zzz_xzzz[k];

                g_xz_0_zzzz_yyy[k] = -g_x_0_zzz_yyy[k] - g_xz_0_zzz_yyy[k] * ab_z + g_xz_0_zzz_yyyz[k];

                g_xz_0_zzzz_yyz[k] = -g_x_0_zzz_yyz[k] - g_xz_0_zzz_yyz[k] * ab_z + g_xz_0_zzz_yyzz[k];

                g_xz_0_zzzz_yzz[k] = -g_x_0_zzz_yzz[k] - g_xz_0_zzz_yzz[k] * ab_z + g_xz_0_zzz_yzzz[k];

                g_xz_0_zzzz_zzz[k] = -g_x_0_zzz_zzz[k] - g_xz_0_zzz_zzz[k] * ab_z + g_xz_0_zzz_zzzz[k];
            }

            /// Set up 450-460 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 450 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 451 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 452 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 453 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 454 * ccomps * dcomps);

            auto g_yy_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 455 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 456 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 457 * ccomps * dcomps);

            auto g_yy_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 458 * ccomps * dcomps);

            auto g_yy_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 459 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxx_xxx, g_yy_0_xxx_xxxx, g_yy_0_xxx_xxxy, g_yy_0_xxx_xxxz, g_yy_0_xxx_xxy, g_yy_0_xxx_xxyy, g_yy_0_xxx_xxyz, g_yy_0_xxx_xxz, g_yy_0_xxx_xxzz, g_yy_0_xxx_xyy, g_yy_0_xxx_xyyy, g_yy_0_xxx_xyyz, g_yy_0_xxx_xyz, g_yy_0_xxx_xyzz, g_yy_0_xxx_xzz, g_yy_0_xxx_xzzz, g_yy_0_xxx_yyy, g_yy_0_xxx_yyz, g_yy_0_xxx_yzz, g_yy_0_xxx_zzz, g_yy_0_xxxx_xxx, g_yy_0_xxxx_xxy, g_yy_0_xxxx_xxz, g_yy_0_xxxx_xyy, g_yy_0_xxxx_xyz, g_yy_0_xxxx_xzz, g_yy_0_xxxx_yyy, g_yy_0_xxxx_yyz, g_yy_0_xxxx_yzz, g_yy_0_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxx_xxx[k] = -g_yy_0_xxx_xxx[k] * ab_x + g_yy_0_xxx_xxxx[k];

                g_yy_0_xxxx_xxy[k] = -g_yy_0_xxx_xxy[k] * ab_x + g_yy_0_xxx_xxxy[k];

                g_yy_0_xxxx_xxz[k] = -g_yy_0_xxx_xxz[k] * ab_x + g_yy_0_xxx_xxxz[k];

                g_yy_0_xxxx_xyy[k] = -g_yy_0_xxx_xyy[k] * ab_x + g_yy_0_xxx_xxyy[k];

                g_yy_0_xxxx_xyz[k] = -g_yy_0_xxx_xyz[k] * ab_x + g_yy_0_xxx_xxyz[k];

                g_yy_0_xxxx_xzz[k] = -g_yy_0_xxx_xzz[k] * ab_x + g_yy_0_xxx_xxzz[k];

                g_yy_0_xxxx_yyy[k] = -g_yy_0_xxx_yyy[k] * ab_x + g_yy_0_xxx_xyyy[k];

                g_yy_0_xxxx_yyz[k] = -g_yy_0_xxx_yyz[k] * ab_x + g_yy_0_xxx_xyyz[k];

                g_yy_0_xxxx_yzz[k] = -g_yy_0_xxx_yzz[k] * ab_x + g_yy_0_xxx_xyzz[k];

                g_yy_0_xxxx_zzz[k] = -g_yy_0_xxx_zzz[k] * ab_x + g_yy_0_xxx_xzzz[k];
            }

            /// Set up 460-470 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 460 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 461 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 462 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 463 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 464 * ccomps * dcomps);

            auto g_yy_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 465 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 466 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 467 * ccomps * dcomps);

            auto g_yy_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 468 * ccomps * dcomps);

            auto g_yy_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 469 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxy_xxx, g_yy_0_xxxy_xxy, g_yy_0_xxxy_xxz, g_yy_0_xxxy_xyy, g_yy_0_xxxy_xyz, g_yy_0_xxxy_xzz, g_yy_0_xxxy_yyy, g_yy_0_xxxy_yyz, g_yy_0_xxxy_yzz, g_yy_0_xxxy_zzz, g_yy_0_xxy_xxx, g_yy_0_xxy_xxxx, g_yy_0_xxy_xxxy, g_yy_0_xxy_xxxz, g_yy_0_xxy_xxy, g_yy_0_xxy_xxyy, g_yy_0_xxy_xxyz, g_yy_0_xxy_xxz, g_yy_0_xxy_xxzz, g_yy_0_xxy_xyy, g_yy_0_xxy_xyyy, g_yy_0_xxy_xyyz, g_yy_0_xxy_xyz, g_yy_0_xxy_xyzz, g_yy_0_xxy_xzz, g_yy_0_xxy_xzzz, g_yy_0_xxy_yyy, g_yy_0_xxy_yyz, g_yy_0_xxy_yzz, g_yy_0_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxy_xxx[k] = -g_yy_0_xxy_xxx[k] * ab_x + g_yy_0_xxy_xxxx[k];

                g_yy_0_xxxy_xxy[k] = -g_yy_0_xxy_xxy[k] * ab_x + g_yy_0_xxy_xxxy[k];

                g_yy_0_xxxy_xxz[k] = -g_yy_0_xxy_xxz[k] * ab_x + g_yy_0_xxy_xxxz[k];

                g_yy_0_xxxy_xyy[k] = -g_yy_0_xxy_xyy[k] * ab_x + g_yy_0_xxy_xxyy[k];

                g_yy_0_xxxy_xyz[k] = -g_yy_0_xxy_xyz[k] * ab_x + g_yy_0_xxy_xxyz[k];

                g_yy_0_xxxy_xzz[k] = -g_yy_0_xxy_xzz[k] * ab_x + g_yy_0_xxy_xxzz[k];

                g_yy_0_xxxy_yyy[k] = -g_yy_0_xxy_yyy[k] * ab_x + g_yy_0_xxy_xyyy[k];

                g_yy_0_xxxy_yyz[k] = -g_yy_0_xxy_yyz[k] * ab_x + g_yy_0_xxy_xyyz[k];

                g_yy_0_xxxy_yzz[k] = -g_yy_0_xxy_yzz[k] * ab_x + g_yy_0_xxy_xyzz[k];

                g_yy_0_xxxy_zzz[k] = -g_yy_0_xxy_zzz[k] * ab_x + g_yy_0_xxy_xzzz[k];
            }

            /// Set up 470-480 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 470 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 471 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 472 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 473 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 474 * ccomps * dcomps);

            auto g_yy_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 475 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 476 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 477 * ccomps * dcomps);

            auto g_yy_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 478 * ccomps * dcomps);

            auto g_yy_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxz_xxx, g_yy_0_xxxz_xxy, g_yy_0_xxxz_xxz, g_yy_0_xxxz_xyy, g_yy_0_xxxz_xyz, g_yy_0_xxxz_xzz, g_yy_0_xxxz_yyy, g_yy_0_xxxz_yyz, g_yy_0_xxxz_yzz, g_yy_0_xxxz_zzz, g_yy_0_xxz_xxx, g_yy_0_xxz_xxxx, g_yy_0_xxz_xxxy, g_yy_0_xxz_xxxz, g_yy_0_xxz_xxy, g_yy_0_xxz_xxyy, g_yy_0_xxz_xxyz, g_yy_0_xxz_xxz, g_yy_0_xxz_xxzz, g_yy_0_xxz_xyy, g_yy_0_xxz_xyyy, g_yy_0_xxz_xyyz, g_yy_0_xxz_xyz, g_yy_0_xxz_xyzz, g_yy_0_xxz_xzz, g_yy_0_xxz_xzzz, g_yy_0_xxz_yyy, g_yy_0_xxz_yyz, g_yy_0_xxz_yzz, g_yy_0_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxz_xxx[k] = -g_yy_0_xxz_xxx[k] * ab_x + g_yy_0_xxz_xxxx[k];

                g_yy_0_xxxz_xxy[k] = -g_yy_0_xxz_xxy[k] * ab_x + g_yy_0_xxz_xxxy[k];

                g_yy_0_xxxz_xxz[k] = -g_yy_0_xxz_xxz[k] * ab_x + g_yy_0_xxz_xxxz[k];

                g_yy_0_xxxz_xyy[k] = -g_yy_0_xxz_xyy[k] * ab_x + g_yy_0_xxz_xxyy[k];

                g_yy_0_xxxz_xyz[k] = -g_yy_0_xxz_xyz[k] * ab_x + g_yy_0_xxz_xxyz[k];

                g_yy_0_xxxz_xzz[k] = -g_yy_0_xxz_xzz[k] * ab_x + g_yy_0_xxz_xxzz[k];

                g_yy_0_xxxz_yyy[k] = -g_yy_0_xxz_yyy[k] * ab_x + g_yy_0_xxz_xyyy[k];

                g_yy_0_xxxz_yyz[k] = -g_yy_0_xxz_yyz[k] * ab_x + g_yy_0_xxz_xyyz[k];

                g_yy_0_xxxz_yzz[k] = -g_yy_0_xxz_yzz[k] * ab_x + g_yy_0_xxz_xyzz[k];

                g_yy_0_xxxz_zzz[k] = -g_yy_0_xxz_zzz[k] * ab_x + g_yy_0_xxz_xzzz[k];
            }

            /// Set up 480-490 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 480 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 481 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 482 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 483 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 484 * ccomps * dcomps);

            auto g_yy_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 485 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 486 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 487 * ccomps * dcomps);

            auto g_yy_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 488 * ccomps * dcomps);

            auto g_yy_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 489 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyy_xxx, g_yy_0_xxyy_xxy, g_yy_0_xxyy_xxz, g_yy_0_xxyy_xyy, g_yy_0_xxyy_xyz, g_yy_0_xxyy_xzz, g_yy_0_xxyy_yyy, g_yy_0_xxyy_yyz, g_yy_0_xxyy_yzz, g_yy_0_xxyy_zzz, g_yy_0_xyy_xxx, g_yy_0_xyy_xxxx, g_yy_0_xyy_xxxy, g_yy_0_xyy_xxxz, g_yy_0_xyy_xxy, g_yy_0_xyy_xxyy, g_yy_0_xyy_xxyz, g_yy_0_xyy_xxz, g_yy_0_xyy_xxzz, g_yy_0_xyy_xyy, g_yy_0_xyy_xyyy, g_yy_0_xyy_xyyz, g_yy_0_xyy_xyz, g_yy_0_xyy_xyzz, g_yy_0_xyy_xzz, g_yy_0_xyy_xzzz, g_yy_0_xyy_yyy, g_yy_0_xyy_yyz, g_yy_0_xyy_yzz, g_yy_0_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyy_xxx[k] = -g_yy_0_xyy_xxx[k] * ab_x + g_yy_0_xyy_xxxx[k];

                g_yy_0_xxyy_xxy[k] = -g_yy_0_xyy_xxy[k] * ab_x + g_yy_0_xyy_xxxy[k];

                g_yy_0_xxyy_xxz[k] = -g_yy_0_xyy_xxz[k] * ab_x + g_yy_0_xyy_xxxz[k];

                g_yy_0_xxyy_xyy[k] = -g_yy_0_xyy_xyy[k] * ab_x + g_yy_0_xyy_xxyy[k];

                g_yy_0_xxyy_xyz[k] = -g_yy_0_xyy_xyz[k] * ab_x + g_yy_0_xyy_xxyz[k];

                g_yy_0_xxyy_xzz[k] = -g_yy_0_xyy_xzz[k] * ab_x + g_yy_0_xyy_xxzz[k];

                g_yy_0_xxyy_yyy[k] = -g_yy_0_xyy_yyy[k] * ab_x + g_yy_0_xyy_xyyy[k];

                g_yy_0_xxyy_yyz[k] = -g_yy_0_xyy_yyz[k] * ab_x + g_yy_0_xyy_xyyz[k];

                g_yy_0_xxyy_yzz[k] = -g_yy_0_xyy_yzz[k] * ab_x + g_yy_0_xyy_xyzz[k];

                g_yy_0_xxyy_zzz[k] = -g_yy_0_xyy_zzz[k] * ab_x + g_yy_0_xyy_xzzz[k];
            }

            /// Set up 490-500 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 490 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 491 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 492 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 493 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 494 * ccomps * dcomps);

            auto g_yy_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 495 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 496 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 497 * ccomps * dcomps);

            auto g_yy_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 498 * ccomps * dcomps);

            auto g_yy_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 499 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyz_xxx, g_yy_0_xxyz_xxy, g_yy_0_xxyz_xxz, g_yy_0_xxyz_xyy, g_yy_0_xxyz_xyz, g_yy_0_xxyz_xzz, g_yy_0_xxyz_yyy, g_yy_0_xxyz_yyz, g_yy_0_xxyz_yzz, g_yy_0_xxyz_zzz, g_yy_0_xyz_xxx, g_yy_0_xyz_xxxx, g_yy_0_xyz_xxxy, g_yy_0_xyz_xxxz, g_yy_0_xyz_xxy, g_yy_0_xyz_xxyy, g_yy_0_xyz_xxyz, g_yy_0_xyz_xxz, g_yy_0_xyz_xxzz, g_yy_0_xyz_xyy, g_yy_0_xyz_xyyy, g_yy_0_xyz_xyyz, g_yy_0_xyz_xyz, g_yy_0_xyz_xyzz, g_yy_0_xyz_xzz, g_yy_0_xyz_xzzz, g_yy_0_xyz_yyy, g_yy_0_xyz_yyz, g_yy_0_xyz_yzz, g_yy_0_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyz_xxx[k] = -g_yy_0_xyz_xxx[k] * ab_x + g_yy_0_xyz_xxxx[k];

                g_yy_0_xxyz_xxy[k] = -g_yy_0_xyz_xxy[k] * ab_x + g_yy_0_xyz_xxxy[k];

                g_yy_0_xxyz_xxz[k] = -g_yy_0_xyz_xxz[k] * ab_x + g_yy_0_xyz_xxxz[k];

                g_yy_0_xxyz_xyy[k] = -g_yy_0_xyz_xyy[k] * ab_x + g_yy_0_xyz_xxyy[k];

                g_yy_0_xxyz_xyz[k] = -g_yy_0_xyz_xyz[k] * ab_x + g_yy_0_xyz_xxyz[k];

                g_yy_0_xxyz_xzz[k] = -g_yy_0_xyz_xzz[k] * ab_x + g_yy_0_xyz_xxzz[k];

                g_yy_0_xxyz_yyy[k] = -g_yy_0_xyz_yyy[k] * ab_x + g_yy_0_xyz_xyyy[k];

                g_yy_0_xxyz_yyz[k] = -g_yy_0_xyz_yyz[k] * ab_x + g_yy_0_xyz_xyyz[k];

                g_yy_0_xxyz_yzz[k] = -g_yy_0_xyz_yzz[k] * ab_x + g_yy_0_xyz_xyzz[k];

                g_yy_0_xxyz_zzz[k] = -g_yy_0_xyz_zzz[k] * ab_x + g_yy_0_xyz_xzzz[k];
            }

            /// Set up 500-510 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 500 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 501 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 502 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 503 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 504 * ccomps * dcomps);

            auto g_yy_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 505 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 506 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 507 * ccomps * dcomps);

            auto g_yy_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 508 * ccomps * dcomps);

            auto g_yy_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxzz_xxx, g_yy_0_xxzz_xxy, g_yy_0_xxzz_xxz, g_yy_0_xxzz_xyy, g_yy_0_xxzz_xyz, g_yy_0_xxzz_xzz, g_yy_0_xxzz_yyy, g_yy_0_xxzz_yyz, g_yy_0_xxzz_yzz, g_yy_0_xxzz_zzz, g_yy_0_xzz_xxx, g_yy_0_xzz_xxxx, g_yy_0_xzz_xxxy, g_yy_0_xzz_xxxz, g_yy_0_xzz_xxy, g_yy_0_xzz_xxyy, g_yy_0_xzz_xxyz, g_yy_0_xzz_xxz, g_yy_0_xzz_xxzz, g_yy_0_xzz_xyy, g_yy_0_xzz_xyyy, g_yy_0_xzz_xyyz, g_yy_0_xzz_xyz, g_yy_0_xzz_xyzz, g_yy_0_xzz_xzz, g_yy_0_xzz_xzzz, g_yy_0_xzz_yyy, g_yy_0_xzz_yyz, g_yy_0_xzz_yzz, g_yy_0_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxzz_xxx[k] = -g_yy_0_xzz_xxx[k] * ab_x + g_yy_0_xzz_xxxx[k];

                g_yy_0_xxzz_xxy[k] = -g_yy_0_xzz_xxy[k] * ab_x + g_yy_0_xzz_xxxy[k];

                g_yy_0_xxzz_xxz[k] = -g_yy_0_xzz_xxz[k] * ab_x + g_yy_0_xzz_xxxz[k];

                g_yy_0_xxzz_xyy[k] = -g_yy_0_xzz_xyy[k] * ab_x + g_yy_0_xzz_xxyy[k];

                g_yy_0_xxzz_xyz[k] = -g_yy_0_xzz_xyz[k] * ab_x + g_yy_0_xzz_xxyz[k];

                g_yy_0_xxzz_xzz[k] = -g_yy_0_xzz_xzz[k] * ab_x + g_yy_0_xzz_xxzz[k];

                g_yy_0_xxzz_yyy[k] = -g_yy_0_xzz_yyy[k] * ab_x + g_yy_0_xzz_xyyy[k];

                g_yy_0_xxzz_yyz[k] = -g_yy_0_xzz_yyz[k] * ab_x + g_yy_0_xzz_xyyz[k];

                g_yy_0_xxzz_yzz[k] = -g_yy_0_xzz_yzz[k] * ab_x + g_yy_0_xzz_xyzz[k];

                g_yy_0_xxzz_zzz[k] = -g_yy_0_xzz_zzz[k] * ab_x + g_yy_0_xzz_xzzz[k];
            }

            /// Set up 510-520 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 510 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 511 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 512 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 513 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 514 * ccomps * dcomps);

            auto g_yy_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 515 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 516 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 517 * ccomps * dcomps);

            auto g_yy_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 518 * ccomps * dcomps);

            auto g_yy_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 519 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyy_xxx, g_yy_0_xyyy_xxy, g_yy_0_xyyy_xxz, g_yy_0_xyyy_xyy, g_yy_0_xyyy_xyz, g_yy_0_xyyy_xzz, g_yy_0_xyyy_yyy, g_yy_0_xyyy_yyz, g_yy_0_xyyy_yzz, g_yy_0_xyyy_zzz, g_yy_0_yyy_xxx, g_yy_0_yyy_xxxx, g_yy_0_yyy_xxxy, g_yy_0_yyy_xxxz, g_yy_0_yyy_xxy, g_yy_0_yyy_xxyy, g_yy_0_yyy_xxyz, g_yy_0_yyy_xxz, g_yy_0_yyy_xxzz, g_yy_0_yyy_xyy, g_yy_0_yyy_xyyy, g_yy_0_yyy_xyyz, g_yy_0_yyy_xyz, g_yy_0_yyy_xyzz, g_yy_0_yyy_xzz, g_yy_0_yyy_xzzz, g_yy_0_yyy_yyy, g_yy_0_yyy_yyz, g_yy_0_yyy_yzz, g_yy_0_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyy_xxx[k] = -g_yy_0_yyy_xxx[k] * ab_x + g_yy_0_yyy_xxxx[k];

                g_yy_0_xyyy_xxy[k] = -g_yy_0_yyy_xxy[k] * ab_x + g_yy_0_yyy_xxxy[k];

                g_yy_0_xyyy_xxz[k] = -g_yy_0_yyy_xxz[k] * ab_x + g_yy_0_yyy_xxxz[k];

                g_yy_0_xyyy_xyy[k] = -g_yy_0_yyy_xyy[k] * ab_x + g_yy_0_yyy_xxyy[k];

                g_yy_0_xyyy_xyz[k] = -g_yy_0_yyy_xyz[k] * ab_x + g_yy_0_yyy_xxyz[k];

                g_yy_0_xyyy_xzz[k] = -g_yy_0_yyy_xzz[k] * ab_x + g_yy_0_yyy_xxzz[k];

                g_yy_0_xyyy_yyy[k] = -g_yy_0_yyy_yyy[k] * ab_x + g_yy_0_yyy_xyyy[k];

                g_yy_0_xyyy_yyz[k] = -g_yy_0_yyy_yyz[k] * ab_x + g_yy_0_yyy_xyyz[k];

                g_yy_0_xyyy_yzz[k] = -g_yy_0_yyy_yzz[k] * ab_x + g_yy_0_yyy_xyzz[k];

                g_yy_0_xyyy_zzz[k] = -g_yy_0_yyy_zzz[k] * ab_x + g_yy_0_yyy_xzzz[k];
            }

            /// Set up 520-530 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 520 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 521 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 522 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 523 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 524 * ccomps * dcomps);

            auto g_yy_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 525 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 526 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 527 * ccomps * dcomps);

            auto g_yy_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 528 * ccomps * dcomps);

            auto g_yy_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 529 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyz_xxx, g_yy_0_xyyz_xxy, g_yy_0_xyyz_xxz, g_yy_0_xyyz_xyy, g_yy_0_xyyz_xyz, g_yy_0_xyyz_xzz, g_yy_0_xyyz_yyy, g_yy_0_xyyz_yyz, g_yy_0_xyyz_yzz, g_yy_0_xyyz_zzz, g_yy_0_yyz_xxx, g_yy_0_yyz_xxxx, g_yy_0_yyz_xxxy, g_yy_0_yyz_xxxz, g_yy_0_yyz_xxy, g_yy_0_yyz_xxyy, g_yy_0_yyz_xxyz, g_yy_0_yyz_xxz, g_yy_0_yyz_xxzz, g_yy_0_yyz_xyy, g_yy_0_yyz_xyyy, g_yy_0_yyz_xyyz, g_yy_0_yyz_xyz, g_yy_0_yyz_xyzz, g_yy_0_yyz_xzz, g_yy_0_yyz_xzzz, g_yy_0_yyz_yyy, g_yy_0_yyz_yyz, g_yy_0_yyz_yzz, g_yy_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyz_xxx[k] = -g_yy_0_yyz_xxx[k] * ab_x + g_yy_0_yyz_xxxx[k];

                g_yy_0_xyyz_xxy[k] = -g_yy_0_yyz_xxy[k] * ab_x + g_yy_0_yyz_xxxy[k];

                g_yy_0_xyyz_xxz[k] = -g_yy_0_yyz_xxz[k] * ab_x + g_yy_0_yyz_xxxz[k];

                g_yy_0_xyyz_xyy[k] = -g_yy_0_yyz_xyy[k] * ab_x + g_yy_0_yyz_xxyy[k];

                g_yy_0_xyyz_xyz[k] = -g_yy_0_yyz_xyz[k] * ab_x + g_yy_0_yyz_xxyz[k];

                g_yy_0_xyyz_xzz[k] = -g_yy_0_yyz_xzz[k] * ab_x + g_yy_0_yyz_xxzz[k];

                g_yy_0_xyyz_yyy[k] = -g_yy_0_yyz_yyy[k] * ab_x + g_yy_0_yyz_xyyy[k];

                g_yy_0_xyyz_yyz[k] = -g_yy_0_yyz_yyz[k] * ab_x + g_yy_0_yyz_xyyz[k];

                g_yy_0_xyyz_yzz[k] = -g_yy_0_yyz_yzz[k] * ab_x + g_yy_0_yyz_xyzz[k];

                g_yy_0_xyyz_zzz[k] = -g_yy_0_yyz_zzz[k] * ab_x + g_yy_0_yyz_xzzz[k];
            }

            /// Set up 530-540 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 530 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 531 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 532 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 533 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 534 * ccomps * dcomps);

            auto g_yy_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 535 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 536 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 537 * ccomps * dcomps);

            auto g_yy_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 538 * ccomps * dcomps);

            auto g_yy_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyzz_xxx, g_yy_0_xyzz_xxy, g_yy_0_xyzz_xxz, g_yy_0_xyzz_xyy, g_yy_0_xyzz_xyz, g_yy_0_xyzz_xzz, g_yy_0_xyzz_yyy, g_yy_0_xyzz_yyz, g_yy_0_xyzz_yzz, g_yy_0_xyzz_zzz, g_yy_0_yzz_xxx, g_yy_0_yzz_xxxx, g_yy_0_yzz_xxxy, g_yy_0_yzz_xxxz, g_yy_0_yzz_xxy, g_yy_0_yzz_xxyy, g_yy_0_yzz_xxyz, g_yy_0_yzz_xxz, g_yy_0_yzz_xxzz, g_yy_0_yzz_xyy, g_yy_0_yzz_xyyy, g_yy_0_yzz_xyyz, g_yy_0_yzz_xyz, g_yy_0_yzz_xyzz, g_yy_0_yzz_xzz, g_yy_0_yzz_xzzz, g_yy_0_yzz_yyy, g_yy_0_yzz_yyz, g_yy_0_yzz_yzz, g_yy_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyzz_xxx[k] = -g_yy_0_yzz_xxx[k] * ab_x + g_yy_0_yzz_xxxx[k];

                g_yy_0_xyzz_xxy[k] = -g_yy_0_yzz_xxy[k] * ab_x + g_yy_0_yzz_xxxy[k];

                g_yy_0_xyzz_xxz[k] = -g_yy_0_yzz_xxz[k] * ab_x + g_yy_0_yzz_xxxz[k];

                g_yy_0_xyzz_xyy[k] = -g_yy_0_yzz_xyy[k] * ab_x + g_yy_0_yzz_xxyy[k];

                g_yy_0_xyzz_xyz[k] = -g_yy_0_yzz_xyz[k] * ab_x + g_yy_0_yzz_xxyz[k];

                g_yy_0_xyzz_xzz[k] = -g_yy_0_yzz_xzz[k] * ab_x + g_yy_0_yzz_xxzz[k];

                g_yy_0_xyzz_yyy[k] = -g_yy_0_yzz_yyy[k] * ab_x + g_yy_0_yzz_xyyy[k];

                g_yy_0_xyzz_yyz[k] = -g_yy_0_yzz_yyz[k] * ab_x + g_yy_0_yzz_xyyz[k];

                g_yy_0_xyzz_yzz[k] = -g_yy_0_yzz_yzz[k] * ab_x + g_yy_0_yzz_xyzz[k];

                g_yy_0_xyzz_zzz[k] = -g_yy_0_yzz_zzz[k] * ab_x + g_yy_0_yzz_xzzz[k];
            }

            /// Set up 540-550 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 540 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 541 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 542 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 543 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 544 * ccomps * dcomps);

            auto g_yy_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 545 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 546 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 547 * ccomps * dcomps);

            auto g_yy_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 548 * ccomps * dcomps);

            auto g_yy_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 549 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzzz_xxx, g_yy_0_xzzz_xxy, g_yy_0_xzzz_xxz, g_yy_0_xzzz_xyy, g_yy_0_xzzz_xyz, g_yy_0_xzzz_xzz, g_yy_0_xzzz_yyy, g_yy_0_xzzz_yyz, g_yy_0_xzzz_yzz, g_yy_0_xzzz_zzz, g_yy_0_zzz_xxx, g_yy_0_zzz_xxxx, g_yy_0_zzz_xxxy, g_yy_0_zzz_xxxz, g_yy_0_zzz_xxy, g_yy_0_zzz_xxyy, g_yy_0_zzz_xxyz, g_yy_0_zzz_xxz, g_yy_0_zzz_xxzz, g_yy_0_zzz_xyy, g_yy_0_zzz_xyyy, g_yy_0_zzz_xyyz, g_yy_0_zzz_xyz, g_yy_0_zzz_xyzz, g_yy_0_zzz_xzz, g_yy_0_zzz_xzzz, g_yy_0_zzz_yyy, g_yy_0_zzz_yyz, g_yy_0_zzz_yzz, g_yy_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzzz_xxx[k] = -g_yy_0_zzz_xxx[k] * ab_x + g_yy_0_zzz_xxxx[k];

                g_yy_0_xzzz_xxy[k] = -g_yy_0_zzz_xxy[k] * ab_x + g_yy_0_zzz_xxxy[k];

                g_yy_0_xzzz_xxz[k] = -g_yy_0_zzz_xxz[k] * ab_x + g_yy_0_zzz_xxxz[k];

                g_yy_0_xzzz_xyy[k] = -g_yy_0_zzz_xyy[k] * ab_x + g_yy_0_zzz_xxyy[k];

                g_yy_0_xzzz_xyz[k] = -g_yy_0_zzz_xyz[k] * ab_x + g_yy_0_zzz_xxyz[k];

                g_yy_0_xzzz_xzz[k] = -g_yy_0_zzz_xzz[k] * ab_x + g_yy_0_zzz_xxzz[k];

                g_yy_0_xzzz_yyy[k] = -g_yy_0_zzz_yyy[k] * ab_x + g_yy_0_zzz_xyyy[k];

                g_yy_0_xzzz_yyz[k] = -g_yy_0_zzz_yyz[k] * ab_x + g_yy_0_zzz_xyyz[k];

                g_yy_0_xzzz_yzz[k] = -g_yy_0_zzz_yzz[k] * ab_x + g_yy_0_zzz_xyzz[k];

                g_yy_0_xzzz_zzz[k] = -g_yy_0_zzz_zzz[k] * ab_x + g_yy_0_zzz_xzzz[k];
            }

            /// Set up 550-560 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 550 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 551 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 552 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 553 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 554 * ccomps * dcomps);

            auto g_yy_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 555 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 556 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 557 * ccomps * dcomps);

            auto g_yy_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 558 * ccomps * dcomps);

            auto g_yy_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_xxx, g_y_0_yyy_xxy, g_y_0_yyy_xxz, g_y_0_yyy_xyy, g_y_0_yyy_xyz, g_y_0_yyy_xzz, g_y_0_yyy_yyy, g_y_0_yyy_yyz, g_y_0_yyy_yzz, g_y_0_yyy_zzz, g_yy_0_yyy_xxx, g_yy_0_yyy_xxxy, g_yy_0_yyy_xxy, g_yy_0_yyy_xxyy, g_yy_0_yyy_xxyz, g_yy_0_yyy_xxz, g_yy_0_yyy_xyy, g_yy_0_yyy_xyyy, g_yy_0_yyy_xyyz, g_yy_0_yyy_xyz, g_yy_0_yyy_xyzz, g_yy_0_yyy_xzz, g_yy_0_yyy_yyy, g_yy_0_yyy_yyyy, g_yy_0_yyy_yyyz, g_yy_0_yyy_yyz, g_yy_0_yyy_yyzz, g_yy_0_yyy_yzz, g_yy_0_yyy_yzzz, g_yy_0_yyy_zzz, g_yy_0_yyyy_xxx, g_yy_0_yyyy_xxy, g_yy_0_yyyy_xxz, g_yy_0_yyyy_xyy, g_yy_0_yyyy_xyz, g_yy_0_yyyy_xzz, g_yy_0_yyyy_yyy, g_yy_0_yyyy_yyz, g_yy_0_yyyy_yzz, g_yy_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyy_xxx[k] = -2.0 * g_y_0_yyy_xxx[k] - g_yy_0_yyy_xxx[k] * ab_y + g_yy_0_yyy_xxxy[k];

                g_yy_0_yyyy_xxy[k] = -2.0 * g_y_0_yyy_xxy[k] - g_yy_0_yyy_xxy[k] * ab_y + g_yy_0_yyy_xxyy[k];

                g_yy_0_yyyy_xxz[k] = -2.0 * g_y_0_yyy_xxz[k] - g_yy_0_yyy_xxz[k] * ab_y + g_yy_0_yyy_xxyz[k];

                g_yy_0_yyyy_xyy[k] = -2.0 * g_y_0_yyy_xyy[k] - g_yy_0_yyy_xyy[k] * ab_y + g_yy_0_yyy_xyyy[k];

                g_yy_0_yyyy_xyz[k] = -2.0 * g_y_0_yyy_xyz[k] - g_yy_0_yyy_xyz[k] * ab_y + g_yy_0_yyy_xyyz[k];

                g_yy_0_yyyy_xzz[k] = -2.0 * g_y_0_yyy_xzz[k] - g_yy_0_yyy_xzz[k] * ab_y + g_yy_0_yyy_xyzz[k];

                g_yy_0_yyyy_yyy[k] = -2.0 * g_y_0_yyy_yyy[k] - g_yy_0_yyy_yyy[k] * ab_y + g_yy_0_yyy_yyyy[k];

                g_yy_0_yyyy_yyz[k] = -2.0 * g_y_0_yyy_yyz[k] - g_yy_0_yyy_yyz[k] * ab_y + g_yy_0_yyy_yyyz[k];

                g_yy_0_yyyy_yzz[k] = -2.0 * g_y_0_yyy_yzz[k] - g_yy_0_yyy_yzz[k] * ab_y + g_yy_0_yyy_yyzz[k];

                g_yy_0_yyyy_zzz[k] = -2.0 * g_y_0_yyy_zzz[k] - g_yy_0_yyy_zzz[k] * ab_y + g_yy_0_yyy_yzzz[k];
            }

            /// Set up 560-570 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 560 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 561 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 562 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 563 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 564 * ccomps * dcomps);

            auto g_yy_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 565 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 566 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 567 * ccomps * dcomps);

            auto g_yy_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 568 * ccomps * dcomps);

            auto g_yy_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyy_xxx, g_yy_0_yyy_xxxz, g_yy_0_yyy_xxy, g_yy_0_yyy_xxyz, g_yy_0_yyy_xxz, g_yy_0_yyy_xxzz, g_yy_0_yyy_xyy, g_yy_0_yyy_xyyz, g_yy_0_yyy_xyz, g_yy_0_yyy_xyzz, g_yy_0_yyy_xzz, g_yy_0_yyy_xzzz, g_yy_0_yyy_yyy, g_yy_0_yyy_yyyz, g_yy_0_yyy_yyz, g_yy_0_yyy_yyzz, g_yy_0_yyy_yzz, g_yy_0_yyy_yzzz, g_yy_0_yyy_zzz, g_yy_0_yyy_zzzz, g_yy_0_yyyz_xxx, g_yy_0_yyyz_xxy, g_yy_0_yyyz_xxz, g_yy_0_yyyz_xyy, g_yy_0_yyyz_xyz, g_yy_0_yyyz_xzz, g_yy_0_yyyz_yyy, g_yy_0_yyyz_yyz, g_yy_0_yyyz_yzz, g_yy_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyz_xxx[k] = -g_yy_0_yyy_xxx[k] * ab_z + g_yy_0_yyy_xxxz[k];

                g_yy_0_yyyz_xxy[k] = -g_yy_0_yyy_xxy[k] * ab_z + g_yy_0_yyy_xxyz[k];

                g_yy_0_yyyz_xxz[k] = -g_yy_0_yyy_xxz[k] * ab_z + g_yy_0_yyy_xxzz[k];

                g_yy_0_yyyz_xyy[k] = -g_yy_0_yyy_xyy[k] * ab_z + g_yy_0_yyy_xyyz[k];

                g_yy_0_yyyz_xyz[k] = -g_yy_0_yyy_xyz[k] * ab_z + g_yy_0_yyy_xyzz[k];

                g_yy_0_yyyz_xzz[k] = -g_yy_0_yyy_xzz[k] * ab_z + g_yy_0_yyy_xzzz[k];

                g_yy_0_yyyz_yyy[k] = -g_yy_0_yyy_yyy[k] * ab_z + g_yy_0_yyy_yyyz[k];

                g_yy_0_yyyz_yyz[k] = -g_yy_0_yyy_yyz[k] * ab_z + g_yy_0_yyy_yyzz[k];

                g_yy_0_yyyz_yzz[k] = -g_yy_0_yyy_yzz[k] * ab_z + g_yy_0_yyy_yzzz[k];

                g_yy_0_yyyz_zzz[k] = -g_yy_0_yyy_zzz[k] * ab_z + g_yy_0_yyy_zzzz[k];
            }

            /// Set up 570-580 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 570 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 571 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 572 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 573 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 574 * ccomps * dcomps);

            auto g_yy_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 575 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 576 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 577 * ccomps * dcomps);

            auto g_yy_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 578 * ccomps * dcomps);

            auto g_yy_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 579 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyz_xxx, g_yy_0_yyz_xxxz, g_yy_0_yyz_xxy, g_yy_0_yyz_xxyz, g_yy_0_yyz_xxz, g_yy_0_yyz_xxzz, g_yy_0_yyz_xyy, g_yy_0_yyz_xyyz, g_yy_0_yyz_xyz, g_yy_0_yyz_xyzz, g_yy_0_yyz_xzz, g_yy_0_yyz_xzzz, g_yy_0_yyz_yyy, g_yy_0_yyz_yyyz, g_yy_0_yyz_yyz, g_yy_0_yyz_yyzz, g_yy_0_yyz_yzz, g_yy_0_yyz_yzzz, g_yy_0_yyz_zzz, g_yy_0_yyz_zzzz, g_yy_0_yyzz_xxx, g_yy_0_yyzz_xxy, g_yy_0_yyzz_xxz, g_yy_0_yyzz_xyy, g_yy_0_yyzz_xyz, g_yy_0_yyzz_xzz, g_yy_0_yyzz_yyy, g_yy_0_yyzz_yyz, g_yy_0_yyzz_yzz, g_yy_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyzz_xxx[k] = -g_yy_0_yyz_xxx[k] * ab_z + g_yy_0_yyz_xxxz[k];

                g_yy_0_yyzz_xxy[k] = -g_yy_0_yyz_xxy[k] * ab_z + g_yy_0_yyz_xxyz[k];

                g_yy_0_yyzz_xxz[k] = -g_yy_0_yyz_xxz[k] * ab_z + g_yy_0_yyz_xxzz[k];

                g_yy_0_yyzz_xyy[k] = -g_yy_0_yyz_xyy[k] * ab_z + g_yy_0_yyz_xyyz[k];

                g_yy_0_yyzz_xyz[k] = -g_yy_0_yyz_xyz[k] * ab_z + g_yy_0_yyz_xyzz[k];

                g_yy_0_yyzz_xzz[k] = -g_yy_0_yyz_xzz[k] * ab_z + g_yy_0_yyz_xzzz[k];

                g_yy_0_yyzz_yyy[k] = -g_yy_0_yyz_yyy[k] * ab_z + g_yy_0_yyz_yyyz[k];

                g_yy_0_yyzz_yyz[k] = -g_yy_0_yyz_yyz[k] * ab_z + g_yy_0_yyz_yyzz[k];

                g_yy_0_yyzz_yzz[k] = -g_yy_0_yyz_yzz[k] * ab_z + g_yy_0_yyz_yzzz[k];

                g_yy_0_yyzz_zzz[k] = -g_yy_0_yyz_zzz[k] * ab_z + g_yy_0_yyz_zzzz[k];
            }

            /// Set up 580-590 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 580 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 581 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 582 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 583 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 584 * ccomps * dcomps);

            auto g_yy_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 585 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 586 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 587 * ccomps * dcomps);

            auto g_yy_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 588 * ccomps * dcomps);

            auto g_yy_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 589 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yzz_xxx, g_yy_0_yzz_xxxz, g_yy_0_yzz_xxy, g_yy_0_yzz_xxyz, g_yy_0_yzz_xxz, g_yy_0_yzz_xxzz, g_yy_0_yzz_xyy, g_yy_0_yzz_xyyz, g_yy_0_yzz_xyz, g_yy_0_yzz_xyzz, g_yy_0_yzz_xzz, g_yy_0_yzz_xzzz, g_yy_0_yzz_yyy, g_yy_0_yzz_yyyz, g_yy_0_yzz_yyz, g_yy_0_yzz_yyzz, g_yy_0_yzz_yzz, g_yy_0_yzz_yzzz, g_yy_0_yzz_zzz, g_yy_0_yzz_zzzz, g_yy_0_yzzz_xxx, g_yy_0_yzzz_xxy, g_yy_0_yzzz_xxz, g_yy_0_yzzz_xyy, g_yy_0_yzzz_xyz, g_yy_0_yzzz_xzz, g_yy_0_yzzz_yyy, g_yy_0_yzzz_yyz, g_yy_0_yzzz_yzz, g_yy_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzzz_xxx[k] = -g_yy_0_yzz_xxx[k] * ab_z + g_yy_0_yzz_xxxz[k];

                g_yy_0_yzzz_xxy[k] = -g_yy_0_yzz_xxy[k] * ab_z + g_yy_0_yzz_xxyz[k];

                g_yy_0_yzzz_xxz[k] = -g_yy_0_yzz_xxz[k] * ab_z + g_yy_0_yzz_xxzz[k];

                g_yy_0_yzzz_xyy[k] = -g_yy_0_yzz_xyy[k] * ab_z + g_yy_0_yzz_xyyz[k];

                g_yy_0_yzzz_xyz[k] = -g_yy_0_yzz_xyz[k] * ab_z + g_yy_0_yzz_xyzz[k];

                g_yy_0_yzzz_xzz[k] = -g_yy_0_yzz_xzz[k] * ab_z + g_yy_0_yzz_xzzz[k];

                g_yy_0_yzzz_yyy[k] = -g_yy_0_yzz_yyy[k] * ab_z + g_yy_0_yzz_yyyz[k];

                g_yy_0_yzzz_yyz[k] = -g_yy_0_yzz_yyz[k] * ab_z + g_yy_0_yzz_yyzz[k];

                g_yy_0_yzzz_yzz[k] = -g_yy_0_yzz_yzz[k] * ab_z + g_yy_0_yzz_yzzz[k];

                g_yy_0_yzzz_zzz[k] = -g_yy_0_yzz_zzz[k] * ab_z + g_yy_0_yzz_zzzz[k];
            }

            /// Set up 590-600 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 590 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 591 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 592 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 593 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 594 * ccomps * dcomps);

            auto g_yy_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 595 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 596 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 597 * ccomps * dcomps);

            auto g_yy_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 598 * ccomps * dcomps);

            auto g_yy_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zzz_xxx, g_yy_0_zzz_xxxz, g_yy_0_zzz_xxy, g_yy_0_zzz_xxyz, g_yy_0_zzz_xxz, g_yy_0_zzz_xxzz, g_yy_0_zzz_xyy, g_yy_0_zzz_xyyz, g_yy_0_zzz_xyz, g_yy_0_zzz_xyzz, g_yy_0_zzz_xzz, g_yy_0_zzz_xzzz, g_yy_0_zzz_yyy, g_yy_0_zzz_yyyz, g_yy_0_zzz_yyz, g_yy_0_zzz_yyzz, g_yy_0_zzz_yzz, g_yy_0_zzz_yzzz, g_yy_0_zzz_zzz, g_yy_0_zzz_zzzz, g_yy_0_zzzz_xxx, g_yy_0_zzzz_xxy, g_yy_0_zzzz_xxz, g_yy_0_zzzz_xyy, g_yy_0_zzzz_xyz, g_yy_0_zzzz_xzz, g_yy_0_zzzz_yyy, g_yy_0_zzzz_yyz, g_yy_0_zzzz_yzz, g_yy_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzzz_xxx[k] = -g_yy_0_zzz_xxx[k] * ab_z + g_yy_0_zzz_xxxz[k];

                g_yy_0_zzzz_xxy[k] = -g_yy_0_zzz_xxy[k] * ab_z + g_yy_0_zzz_xxyz[k];

                g_yy_0_zzzz_xxz[k] = -g_yy_0_zzz_xxz[k] * ab_z + g_yy_0_zzz_xxzz[k];

                g_yy_0_zzzz_xyy[k] = -g_yy_0_zzz_xyy[k] * ab_z + g_yy_0_zzz_xyyz[k];

                g_yy_0_zzzz_xyz[k] = -g_yy_0_zzz_xyz[k] * ab_z + g_yy_0_zzz_xyzz[k];

                g_yy_0_zzzz_xzz[k] = -g_yy_0_zzz_xzz[k] * ab_z + g_yy_0_zzz_xzzz[k];

                g_yy_0_zzzz_yyy[k] = -g_yy_0_zzz_yyy[k] * ab_z + g_yy_0_zzz_yyyz[k];

                g_yy_0_zzzz_yyz[k] = -g_yy_0_zzz_yyz[k] * ab_z + g_yy_0_zzz_yyzz[k];

                g_yy_0_zzzz_yzz[k] = -g_yy_0_zzz_yzz[k] * ab_z + g_yy_0_zzz_yzzz[k];

                g_yy_0_zzzz_zzz[k] = -g_yy_0_zzz_zzz[k] * ab_z + g_yy_0_zzz_zzzz[k];
            }

            /// Set up 600-610 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 600 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 601 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 602 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 603 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 604 * ccomps * dcomps);

            auto g_yz_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 605 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 606 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 607 * ccomps * dcomps);

            auto g_yz_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 608 * ccomps * dcomps);

            auto g_yz_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 609 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxx_xxx, g_yz_0_xxx_xxxx, g_yz_0_xxx_xxxy, g_yz_0_xxx_xxxz, g_yz_0_xxx_xxy, g_yz_0_xxx_xxyy, g_yz_0_xxx_xxyz, g_yz_0_xxx_xxz, g_yz_0_xxx_xxzz, g_yz_0_xxx_xyy, g_yz_0_xxx_xyyy, g_yz_0_xxx_xyyz, g_yz_0_xxx_xyz, g_yz_0_xxx_xyzz, g_yz_0_xxx_xzz, g_yz_0_xxx_xzzz, g_yz_0_xxx_yyy, g_yz_0_xxx_yyz, g_yz_0_xxx_yzz, g_yz_0_xxx_zzz, g_yz_0_xxxx_xxx, g_yz_0_xxxx_xxy, g_yz_0_xxxx_xxz, g_yz_0_xxxx_xyy, g_yz_0_xxxx_xyz, g_yz_0_xxxx_xzz, g_yz_0_xxxx_yyy, g_yz_0_xxxx_yyz, g_yz_0_xxxx_yzz, g_yz_0_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxx_xxx[k] = -g_yz_0_xxx_xxx[k] * ab_x + g_yz_0_xxx_xxxx[k];

                g_yz_0_xxxx_xxy[k] = -g_yz_0_xxx_xxy[k] * ab_x + g_yz_0_xxx_xxxy[k];

                g_yz_0_xxxx_xxz[k] = -g_yz_0_xxx_xxz[k] * ab_x + g_yz_0_xxx_xxxz[k];

                g_yz_0_xxxx_xyy[k] = -g_yz_0_xxx_xyy[k] * ab_x + g_yz_0_xxx_xxyy[k];

                g_yz_0_xxxx_xyz[k] = -g_yz_0_xxx_xyz[k] * ab_x + g_yz_0_xxx_xxyz[k];

                g_yz_0_xxxx_xzz[k] = -g_yz_0_xxx_xzz[k] * ab_x + g_yz_0_xxx_xxzz[k];

                g_yz_0_xxxx_yyy[k] = -g_yz_0_xxx_yyy[k] * ab_x + g_yz_0_xxx_xyyy[k];

                g_yz_0_xxxx_yyz[k] = -g_yz_0_xxx_yyz[k] * ab_x + g_yz_0_xxx_xyyz[k];

                g_yz_0_xxxx_yzz[k] = -g_yz_0_xxx_yzz[k] * ab_x + g_yz_0_xxx_xyzz[k];

                g_yz_0_xxxx_zzz[k] = -g_yz_0_xxx_zzz[k] * ab_x + g_yz_0_xxx_xzzz[k];
            }

            /// Set up 610-620 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 610 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 611 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 612 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 613 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 614 * ccomps * dcomps);

            auto g_yz_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 615 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 616 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 617 * ccomps * dcomps);

            auto g_yz_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 618 * ccomps * dcomps);

            auto g_yz_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 619 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxy_xxx, g_yz_0_xxxy_xxy, g_yz_0_xxxy_xxz, g_yz_0_xxxy_xyy, g_yz_0_xxxy_xyz, g_yz_0_xxxy_xzz, g_yz_0_xxxy_yyy, g_yz_0_xxxy_yyz, g_yz_0_xxxy_yzz, g_yz_0_xxxy_zzz, g_yz_0_xxy_xxx, g_yz_0_xxy_xxxx, g_yz_0_xxy_xxxy, g_yz_0_xxy_xxxz, g_yz_0_xxy_xxy, g_yz_0_xxy_xxyy, g_yz_0_xxy_xxyz, g_yz_0_xxy_xxz, g_yz_0_xxy_xxzz, g_yz_0_xxy_xyy, g_yz_0_xxy_xyyy, g_yz_0_xxy_xyyz, g_yz_0_xxy_xyz, g_yz_0_xxy_xyzz, g_yz_0_xxy_xzz, g_yz_0_xxy_xzzz, g_yz_0_xxy_yyy, g_yz_0_xxy_yyz, g_yz_0_xxy_yzz, g_yz_0_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxy_xxx[k] = -g_yz_0_xxy_xxx[k] * ab_x + g_yz_0_xxy_xxxx[k];

                g_yz_0_xxxy_xxy[k] = -g_yz_0_xxy_xxy[k] * ab_x + g_yz_0_xxy_xxxy[k];

                g_yz_0_xxxy_xxz[k] = -g_yz_0_xxy_xxz[k] * ab_x + g_yz_0_xxy_xxxz[k];

                g_yz_0_xxxy_xyy[k] = -g_yz_0_xxy_xyy[k] * ab_x + g_yz_0_xxy_xxyy[k];

                g_yz_0_xxxy_xyz[k] = -g_yz_0_xxy_xyz[k] * ab_x + g_yz_0_xxy_xxyz[k];

                g_yz_0_xxxy_xzz[k] = -g_yz_0_xxy_xzz[k] * ab_x + g_yz_0_xxy_xxzz[k];

                g_yz_0_xxxy_yyy[k] = -g_yz_0_xxy_yyy[k] * ab_x + g_yz_0_xxy_xyyy[k];

                g_yz_0_xxxy_yyz[k] = -g_yz_0_xxy_yyz[k] * ab_x + g_yz_0_xxy_xyyz[k];

                g_yz_0_xxxy_yzz[k] = -g_yz_0_xxy_yzz[k] * ab_x + g_yz_0_xxy_xyzz[k];

                g_yz_0_xxxy_zzz[k] = -g_yz_0_xxy_zzz[k] * ab_x + g_yz_0_xxy_xzzz[k];
            }

            /// Set up 620-630 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 620 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 621 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 622 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 623 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 624 * ccomps * dcomps);

            auto g_yz_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 625 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 626 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 627 * ccomps * dcomps);

            auto g_yz_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 628 * ccomps * dcomps);

            auto g_yz_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxz_xxx, g_yz_0_xxxz_xxy, g_yz_0_xxxz_xxz, g_yz_0_xxxz_xyy, g_yz_0_xxxz_xyz, g_yz_0_xxxz_xzz, g_yz_0_xxxz_yyy, g_yz_0_xxxz_yyz, g_yz_0_xxxz_yzz, g_yz_0_xxxz_zzz, g_yz_0_xxz_xxx, g_yz_0_xxz_xxxx, g_yz_0_xxz_xxxy, g_yz_0_xxz_xxxz, g_yz_0_xxz_xxy, g_yz_0_xxz_xxyy, g_yz_0_xxz_xxyz, g_yz_0_xxz_xxz, g_yz_0_xxz_xxzz, g_yz_0_xxz_xyy, g_yz_0_xxz_xyyy, g_yz_0_xxz_xyyz, g_yz_0_xxz_xyz, g_yz_0_xxz_xyzz, g_yz_0_xxz_xzz, g_yz_0_xxz_xzzz, g_yz_0_xxz_yyy, g_yz_0_xxz_yyz, g_yz_0_xxz_yzz, g_yz_0_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxz_xxx[k] = -g_yz_0_xxz_xxx[k] * ab_x + g_yz_0_xxz_xxxx[k];

                g_yz_0_xxxz_xxy[k] = -g_yz_0_xxz_xxy[k] * ab_x + g_yz_0_xxz_xxxy[k];

                g_yz_0_xxxz_xxz[k] = -g_yz_0_xxz_xxz[k] * ab_x + g_yz_0_xxz_xxxz[k];

                g_yz_0_xxxz_xyy[k] = -g_yz_0_xxz_xyy[k] * ab_x + g_yz_0_xxz_xxyy[k];

                g_yz_0_xxxz_xyz[k] = -g_yz_0_xxz_xyz[k] * ab_x + g_yz_0_xxz_xxyz[k];

                g_yz_0_xxxz_xzz[k] = -g_yz_0_xxz_xzz[k] * ab_x + g_yz_0_xxz_xxzz[k];

                g_yz_0_xxxz_yyy[k] = -g_yz_0_xxz_yyy[k] * ab_x + g_yz_0_xxz_xyyy[k];

                g_yz_0_xxxz_yyz[k] = -g_yz_0_xxz_yyz[k] * ab_x + g_yz_0_xxz_xyyz[k];

                g_yz_0_xxxz_yzz[k] = -g_yz_0_xxz_yzz[k] * ab_x + g_yz_0_xxz_xyzz[k];

                g_yz_0_xxxz_zzz[k] = -g_yz_0_xxz_zzz[k] * ab_x + g_yz_0_xxz_xzzz[k];
            }

            /// Set up 630-640 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 630 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 631 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 632 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 633 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 634 * ccomps * dcomps);

            auto g_yz_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 635 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 636 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 637 * ccomps * dcomps);

            auto g_yz_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 638 * ccomps * dcomps);

            auto g_yz_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 639 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyy_xxx, g_yz_0_xxyy_xxy, g_yz_0_xxyy_xxz, g_yz_0_xxyy_xyy, g_yz_0_xxyy_xyz, g_yz_0_xxyy_xzz, g_yz_0_xxyy_yyy, g_yz_0_xxyy_yyz, g_yz_0_xxyy_yzz, g_yz_0_xxyy_zzz, g_yz_0_xyy_xxx, g_yz_0_xyy_xxxx, g_yz_0_xyy_xxxy, g_yz_0_xyy_xxxz, g_yz_0_xyy_xxy, g_yz_0_xyy_xxyy, g_yz_0_xyy_xxyz, g_yz_0_xyy_xxz, g_yz_0_xyy_xxzz, g_yz_0_xyy_xyy, g_yz_0_xyy_xyyy, g_yz_0_xyy_xyyz, g_yz_0_xyy_xyz, g_yz_0_xyy_xyzz, g_yz_0_xyy_xzz, g_yz_0_xyy_xzzz, g_yz_0_xyy_yyy, g_yz_0_xyy_yyz, g_yz_0_xyy_yzz, g_yz_0_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyy_xxx[k] = -g_yz_0_xyy_xxx[k] * ab_x + g_yz_0_xyy_xxxx[k];

                g_yz_0_xxyy_xxy[k] = -g_yz_0_xyy_xxy[k] * ab_x + g_yz_0_xyy_xxxy[k];

                g_yz_0_xxyy_xxz[k] = -g_yz_0_xyy_xxz[k] * ab_x + g_yz_0_xyy_xxxz[k];

                g_yz_0_xxyy_xyy[k] = -g_yz_0_xyy_xyy[k] * ab_x + g_yz_0_xyy_xxyy[k];

                g_yz_0_xxyy_xyz[k] = -g_yz_0_xyy_xyz[k] * ab_x + g_yz_0_xyy_xxyz[k];

                g_yz_0_xxyy_xzz[k] = -g_yz_0_xyy_xzz[k] * ab_x + g_yz_0_xyy_xxzz[k];

                g_yz_0_xxyy_yyy[k] = -g_yz_0_xyy_yyy[k] * ab_x + g_yz_0_xyy_xyyy[k];

                g_yz_0_xxyy_yyz[k] = -g_yz_0_xyy_yyz[k] * ab_x + g_yz_0_xyy_xyyz[k];

                g_yz_0_xxyy_yzz[k] = -g_yz_0_xyy_yzz[k] * ab_x + g_yz_0_xyy_xyzz[k];

                g_yz_0_xxyy_zzz[k] = -g_yz_0_xyy_zzz[k] * ab_x + g_yz_0_xyy_xzzz[k];
            }

            /// Set up 640-650 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 640 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 641 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 642 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 643 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 644 * ccomps * dcomps);

            auto g_yz_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 645 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 646 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 647 * ccomps * dcomps);

            auto g_yz_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 648 * ccomps * dcomps);

            auto g_yz_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 649 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyz_xxx, g_yz_0_xxyz_xxy, g_yz_0_xxyz_xxz, g_yz_0_xxyz_xyy, g_yz_0_xxyz_xyz, g_yz_0_xxyz_xzz, g_yz_0_xxyz_yyy, g_yz_0_xxyz_yyz, g_yz_0_xxyz_yzz, g_yz_0_xxyz_zzz, g_yz_0_xyz_xxx, g_yz_0_xyz_xxxx, g_yz_0_xyz_xxxy, g_yz_0_xyz_xxxz, g_yz_0_xyz_xxy, g_yz_0_xyz_xxyy, g_yz_0_xyz_xxyz, g_yz_0_xyz_xxz, g_yz_0_xyz_xxzz, g_yz_0_xyz_xyy, g_yz_0_xyz_xyyy, g_yz_0_xyz_xyyz, g_yz_0_xyz_xyz, g_yz_0_xyz_xyzz, g_yz_0_xyz_xzz, g_yz_0_xyz_xzzz, g_yz_0_xyz_yyy, g_yz_0_xyz_yyz, g_yz_0_xyz_yzz, g_yz_0_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyz_xxx[k] = -g_yz_0_xyz_xxx[k] * ab_x + g_yz_0_xyz_xxxx[k];

                g_yz_0_xxyz_xxy[k] = -g_yz_0_xyz_xxy[k] * ab_x + g_yz_0_xyz_xxxy[k];

                g_yz_0_xxyz_xxz[k] = -g_yz_0_xyz_xxz[k] * ab_x + g_yz_0_xyz_xxxz[k];

                g_yz_0_xxyz_xyy[k] = -g_yz_0_xyz_xyy[k] * ab_x + g_yz_0_xyz_xxyy[k];

                g_yz_0_xxyz_xyz[k] = -g_yz_0_xyz_xyz[k] * ab_x + g_yz_0_xyz_xxyz[k];

                g_yz_0_xxyz_xzz[k] = -g_yz_0_xyz_xzz[k] * ab_x + g_yz_0_xyz_xxzz[k];

                g_yz_0_xxyz_yyy[k] = -g_yz_0_xyz_yyy[k] * ab_x + g_yz_0_xyz_xyyy[k];

                g_yz_0_xxyz_yyz[k] = -g_yz_0_xyz_yyz[k] * ab_x + g_yz_0_xyz_xyyz[k];

                g_yz_0_xxyz_yzz[k] = -g_yz_0_xyz_yzz[k] * ab_x + g_yz_0_xyz_xyzz[k];

                g_yz_0_xxyz_zzz[k] = -g_yz_0_xyz_zzz[k] * ab_x + g_yz_0_xyz_xzzz[k];
            }

            /// Set up 650-660 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 650 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 651 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 652 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 653 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 654 * ccomps * dcomps);

            auto g_yz_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 655 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 656 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 657 * ccomps * dcomps);

            auto g_yz_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 658 * ccomps * dcomps);

            auto g_yz_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxzz_xxx, g_yz_0_xxzz_xxy, g_yz_0_xxzz_xxz, g_yz_0_xxzz_xyy, g_yz_0_xxzz_xyz, g_yz_0_xxzz_xzz, g_yz_0_xxzz_yyy, g_yz_0_xxzz_yyz, g_yz_0_xxzz_yzz, g_yz_0_xxzz_zzz, g_yz_0_xzz_xxx, g_yz_0_xzz_xxxx, g_yz_0_xzz_xxxy, g_yz_0_xzz_xxxz, g_yz_0_xzz_xxy, g_yz_0_xzz_xxyy, g_yz_0_xzz_xxyz, g_yz_0_xzz_xxz, g_yz_0_xzz_xxzz, g_yz_0_xzz_xyy, g_yz_0_xzz_xyyy, g_yz_0_xzz_xyyz, g_yz_0_xzz_xyz, g_yz_0_xzz_xyzz, g_yz_0_xzz_xzz, g_yz_0_xzz_xzzz, g_yz_0_xzz_yyy, g_yz_0_xzz_yyz, g_yz_0_xzz_yzz, g_yz_0_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxzz_xxx[k] = -g_yz_0_xzz_xxx[k] * ab_x + g_yz_0_xzz_xxxx[k];

                g_yz_0_xxzz_xxy[k] = -g_yz_0_xzz_xxy[k] * ab_x + g_yz_0_xzz_xxxy[k];

                g_yz_0_xxzz_xxz[k] = -g_yz_0_xzz_xxz[k] * ab_x + g_yz_0_xzz_xxxz[k];

                g_yz_0_xxzz_xyy[k] = -g_yz_0_xzz_xyy[k] * ab_x + g_yz_0_xzz_xxyy[k];

                g_yz_0_xxzz_xyz[k] = -g_yz_0_xzz_xyz[k] * ab_x + g_yz_0_xzz_xxyz[k];

                g_yz_0_xxzz_xzz[k] = -g_yz_0_xzz_xzz[k] * ab_x + g_yz_0_xzz_xxzz[k];

                g_yz_0_xxzz_yyy[k] = -g_yz_0_xzz_yyy[k] * ab_x + g_yz_0_xzz_xyyy[k];

                g_yz_0_xxzz_yyz[k] = -g_yz_0_xzz_yyz[k] * ab_x + g_yz_0_xzz_xyyz[k];

                g_yz_0_xxzz_yzz[k] = -g_yz_0_xzz_yzz[k] * ab_x + g_yz_0_xzz_xyzz[k];

                g_yz_0_xxzz_zzz[k] = -g_yz_0_xzz_zzz[k] * ab_x + g_yz_0_xzz_xzzz[k];
            }

            /// Set up 660-670 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 660 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 661 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 662 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 663 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 664 * ccomps * dcomps);

            auto g_yz_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 665 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 666 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 667 * ccomps * dcomps);

            auto g_yz_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 668 * ccomps * dcomps);

            auto g_yz_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 669 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyy_xxx, g_yz_0_xyyy_xxy, g_yz_0_xyyy_xxz, g_yz_0_xyyy_xyy, g_yz_0_xyyy_xyz, g_yz_0_xyyy_xzz, g_yz_0_xyyy_yyy, g_yz_0_xyyy_yyz, g_yz_0_xyyy_yzz, g_yz_0_xyyy_zzz, g_yz_0_yyy_xxx, g_yz_0_yyy_xxxx, g_yz_0_yyy_xxxy, g_yz_0_yyy_xxxz, g_yz_0_yyy_xxy, g_yz_0_yyy_xxyy, g_yz_0_yyy_xxyz, g_yz_0_yyy_xxz, g_yz_0_yyy_xxzz, g_yz_0_yyy_xyy, g_yz_0_yyy_xyyy, g_yz_0_yyy_xyyz, g_yz_0_yyy_xyz, g_yz_0_yyy_xyzz, g_yz_0_yyy_xzz, g_yz_0_yyy_xzzz, g_yz_0_yyy_yyy, g_yz_0_yyy_yyz, g_yz_0_yyy_yzz, g_yz_0_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyy_xxx[k] = -g_yz_0_yyy_xxx[k] * ab_x + g_yz_0_yyy_xxxx[k];

                g_yz_0_xyyy_xxy[k] = -g_yz_0_yyy_xxy[k] * ab_x + g_yz_0_yyy_xxxy[k];

                g_yz_0_xyyy_xxz[k] = -g_yz_0_yyy_xxz[k] * ab_x + g_yz_0_yyy_xxxz[k];

                g_yz_0_xyyy_xyy[k] = -g_yz_0_yyy_xyy[k] * ab_x + g_yz_0_yyy_xxyy[k];

                g_yz_0_xyyy_xyz[k] = -g_yz_0_yyy_xyz[k] * ab_x + g_yz_0_yyy_xxyz[k];

                g_yz_0_xyyy_xzz[k] = -g_yz_0_yyy_xzz[k] * ab_x + g_yz_0_yyy_xxzz[k];

                g_yz_0_xyyy_yyy[k] = -g_yz_0_yyy_yyy[k] * ab_x + g_yz_0_yyy_xyyy[k];

                g_yz_0_xyyy_yyz[k] = -g_yz_0_yyy_yyz[k] * ab_x + g_yz_0_yyy_xyyz[k];

                g_yz_0_xyyy_yzz[k] = -g_yz_0_yyy_yzz[k] * ab_x + g_yz_0_yyy_xyzz[k];

                g_yz_0_xyyy_zzz[k] = -g_yz_0_yyy_zzz[k] * ab_x + g_yz_0_yyy_xzzz[k];
            }

            /// Set up 670-680 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 670 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 671 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 672 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 673 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 674 * ccomps * dcomps);

            auto g_yz_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 675 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 676 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 677 * ccomps * dcomps);

            auto g_yz_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 678 * ccomps * dcomps);

            auto g_yz_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 679 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyz_xxx, g_yz_0_xyyz_xxy, g_yz_0_xyyz_xxz, g_yz_0_xyyz_xyy, g_yz_0_xyyz_xyz, g_yz_0_xyyz_xzz, g_yz_0_xyyz_yyy, g_yz_0_xyyz_yyz, g_yz_0_xyyz_yzz, g_yz_0_xyyz_zzz, g_yz_0_yyz_xxx, g_yz_0_yyz_xxxx, g_yz_0_yyz_xxxy, g_yz_0_yyz_xxxz, g_yz_0_yyz_xxy, g_yz_0_yyz_xxyy, g_yz_0_yyz_xxyz, g_yz_0_yyz_xxz, g_yz_0_yyz_xxzz, g_yz_0_yyz_xyy, g_yz_0_yyz_xyyy, g_yz_0_yyz_xyyz, g_yz_0_yyz_xyz, g_yz_0_yyz_xyzz, g_yz_0_yyz_xzz, g_yz_0_yyz_xzzz, g_yz_0_yyz_yyy, g_yz_0_yyz_yyz, g_yz_0_yyz_yzz, g_yz_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyz_xxx[k] = -g_yz_0_yyz_xxx[k] * ab_x + g_yz_0_yyz_xxxx[k];

                g_yz_0_xyyz_xxy[k] = -g_yz_0_yyz_xxy[k] * ab_x + g_yz_0_yyz_xxxy[k];

                g_yz_0_xyyz_xxz[k] = -g_yz_0_yyz_xxz[k] * ab_x + g_yz_0_yyz_xxxz[k];

                g_yz_0_xyyz_xyy[k] = -g_yz_0_yyz_xyy[k] * ab_x + g_yz_0_yyz_xxyy[k];

                g_yz_0_xyyz_xyz[k] = -g_yz_0_yyz_xyz[k] * ab_x + g_yz_0_yyz_xxyz[k];

                g_yz_0_xyyz_xzz[k] = -g_yz_0_yyz_xzz[k] * ab_x + g_yz_0_yyz_xxzz[k];

                g_yz_0_xyyz_yyy[k] = -g_yz_0_yyz_yyy[k] * ab_x + g_yz_0_yyz_xyyy[k];

                g_yz_0_xyyz_yyz[k] = -g_yz_0_yyz_yyz[k] * ab_x + g_yz_0_yyz_xyyz[k];

                g_yz_0_xyyz_yzz[k] = -g_yz_0_yyz_yzz[k] * ab_x + g_yz_0_yyz_xyzz[k];

                g_yz_0_xyyz_zzz[k] = -g_yz_0_yyz_zzz[k] * ab_x + g_yz_0_yyz_xzzz[k];
            }

            /// Set up 680-690 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 680 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 681 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 682 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 683 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 684 * ccomps * dcomps);

            auto g_yz_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 685 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 686 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 687 * ccomps * dcomps);

            auto g_yz_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 688 * ccomps * dcomps);

            auto g_yz_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyzz_xxx, g_yz_0_xyzz_xxy, g_yz_0_xyzz_xxz, g_yz_0_xyzz_xyy, g_yz_0_xyzz_xyz, g_yz_0_xyzz_xzz, g_yz_0_xyzz_yyy, g_yz_0_xyzz_yyz, g_yz_0_xyzz_yzz, g_yz_0_xyzz_zzz, g_yz_0_yzz_xxx, g_yz_0_yzz_xxxx, g_yz_0_yzz_xxxy, g_yz_0_yzz_xxxz, g_yz_0_yzz_xxy, g_yz_0_yzz_xxyy, g_yz_0_yzz_xxyz, g_yz_0_yzz_xxz, g_yz_0_yzz_xxzz, g_yz_0_yzz_xyy, g_yz_0_yzz_xyyy, g_yz_0_yzz_xyyz, g_yz_0_yzz_xyz, g_yz_0_yzz_xyzz, g_yz_0_yzz_xzz, g_yz_0_yzz_xzzz, g_yz_0_yzz_yyy, g_yz_0_yzz_yyz, g_yz_0_yzz_yzz, g_yz_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyzz_xxx[k] = -g_yz_0_yzz_xxx[k] * ab_x + g_yz_0_yzz_xxxx[k];

                g_yz_0_xyzz_xxy[k] = -g_yz_0_yzz_xxy[k] * ab_x + g_yz_0_yzz_xxxy[k];

                g_yz_0_xyzz_xxz[k] = -g_yz_0_yzz_xxz[k] * ab_x + g_yz_0_yzz_xxxz[k];

                g_yz_0_xyzz_xyy[k] = -g_yz_0_yzz_xyy[k] * ab_x + g_yz_0_yzz_xxyy[k];

                g_yz_0_xyzz_xyz[k] = -g_yz_0_yzz_xyz[k] * ab_x + g_yz_0_yzz_xxyz[k];

                g_yz_0_xyzz_xzz[k] = -g_yz_0_yzz_xzz[k] * ab_x + g_yz_0_yzz_xxzz[k];

                g_yz_0_xyzz_yyy[k] = -g_yz_0_yzz_yyy[k] * ab_x + g_yz_0_yzz_xyyy[k];

                g_yz_0_xyzz_yyz[k] = -g_yz_0_yzz_yyz[k] * ab_x + g_yz_0_yzz_xyyz[k];

                g_yz_0_xyzz_yzz[k] = -g_yz_0_yzz_yzz[k] * ab_x + g_yz_0_yzz_xyzz[k];

                g_yz_0_xyzz_zzz[k] = -g_yz_0_yzz_zzz[k] * ab_x + g_yz_0_yzz_xzzz[k];
            }

            /// Set up 690-700 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 690 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 691 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 692 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 693 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 694 * ccomps * dcomps);

            auto g_yz_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 695 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 696 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 697 * ccomps * dcomps);

            auto g_yz_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 698 * ccomps * dcomps);

            auto g_yz_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzzz_xxx, g_yz_0_xzzz_xxy, g_yz_0_xzzz_xxz, g_yz_0_xzzz_xyy, g_yz_0_xzzz_xyz, g_yz_0_xzzz_xzz, g_yz_0_xzzz_yyy, g_yz_0_xzzz_yyz, g_yz_0_xzzz_yzz, g_yz_0_xzzz_zzz, g_yz_0_zzz_xxx, g_yz_0_zzz_xxxx, g_yz_0_zzz_xxxy, g_yz_0_zzz_xxxz, g_yz_0_zzz_xxy, g_yz_0_zzz_xxyy, g_yz_0_zzz_xxyz, g_yz_0_zzz_xxz, g_yz_0_zzz_xxzz, g_yz_0_zzz_xyy, g_yz_0_zzz_xyyy, g_yz_0_zzz_xyyz, g_yz_0_zzz_xyz, g_yz_0_zzz_xyzz, g_yz_0_zzz_xzz, g_yz_0_zzz_xzzz, g_yz_0_zzz_yyy, g_yz_0_zzz_yyz, g_yz_0_zzz_yzz, g_yz_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzzz_xxx[k] = -g_yz_0_zzz_xxx[k] * ab_x + g_yz_0_zzz_xxxx[k];

                g_yz_0_xzzz_xxy[k] = -g_yz_0_zzz_xxy[k] * ab_x + g_yz_0_zzz_xxxy[k];

                g_yz_0_xzzz_xxz[k] = -g_yz_0_zzz_xxz[k] * ab_x + g_yz_0_zzz_xxxz[k];

                g_yz_0_xzzz_xyy[k] = -g_yz_0_zzz_xyy[k] * ab_x + g_yz_0_zzz_xxyy[k];

                g_yz_0_xzzz_xyz[k] = -g_yz_0_zzz_xyz[k] * ab_x + g_yz_0_zzz_xxyz[k];

                g_yz_0_xzzz_xzz[k] = -g_yz_0_zzz_xzz[k] * ab_x + g_yz_0_zzz_xxzz[k];

                g_yz_0_xzzz_yyy[k] = -g_yz_0_zzz_yyy[k] * ab_x + g_yz_0_zzz_xyyy[k];

                g_yz_0_xzzz_yyz[k] = -g_yz_0_zzz_yyz[k] * ab_x + g_yz_0_zzz_xyyz[k];

                g_yz_0_xzzz_yzz[k] = -g_yz_0_zzz_yzz[k] * ab_x + g_yz_0_zzz_xyzz[k];

                g_yz_0_xzzz_zzz[k] = -g_yz_0_zzz_zzz[k] * ab_x + g_yz_0_zzz_xzzz[k];
            }

            /// Set up 700-710 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 700 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 701 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 702 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 703 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 704 * ccomps * dcomps);

            auto g_yz_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 705 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 706 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 707 * ccomps * dcomps);

            auto g_yz_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 708 * ccomps * dcomps);

            auto g_yz_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 709 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyy_xxx, g_yz_0_yyy_xxxy, g_yz_0_yyy_xxy, g_yz_0_yyy_xxyy, g_yz_0_yyy_xxyz, g_yz_0_yyy_xxz, g_yz_0_yyy_xyy, g_yz_0_yyy_xyyy, g_yz_0_yyy_xyyz, g_yz_0_yyy_xyz, g_yz_0_yyy_xyzz, g_yz_0_yyy_xzz, g_yz_0_yyy_yyy, g_yz_0_yyy_yyyy, g_yz_0_yyy_yyyz, g_yz_0_yyy_yyz, g_yz_0_yyy_yyzz, g_yz_0_yyy_yzz, g_yz_0_yyy_yzzz, g_yz_0_yyy_zzz, g_yz_0_yyyy_xxx, g_yz_0_yyyy_xxy, g_yz_0_yyyy_xxz, g_yz_0_yyyy_xyy, g_yz_0_yyyy_xyz, g_yz_0_yyyy_xzz, g_yz_0_yyyy_yyy, g_yz_0_yyyy_yyz, g_yz_0_yyyy_yzz, g_yz_0_yyyy_zzz, g_z_0_yyy_xxx, g_z_0_yyy_xxy, g_z_0_yyy_xxz, g_z_0_yyy_xyy, g_z_0_yyy_xyz, g_z_0_yyy_xzz, g_z_0_yyy_yyy, g_z_0_yyy_yyz, g_z_0_yyy_yzz, g_z_0_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyy_xxx[k] = -g_z_0_yyy_xxx[k] - g_yz_0_yyy_xxx[k] * ab_y + g_yz_0_yyy_xxxy[k];

                g_yz_0_yyyy_xxy[k] = -g_z_0_yyy_xxy[k] - g_yz_0_yyy_xxy[k] * ab_y + g_yz_0_yyy_xxyy[k];

                g_yz_0_yyyy_xxz[k] = -g_z_0_yyy_xxz[k] - g_yz_0_yyy_xxz[k] * ab_y + g_yz_0_yyy_xxyz[k];

                g_yz_0_yyyy_xyy[k] = -g_z_0_yyy_xyy[k] - g_yz_0_yyy_xyy[k] * ab_y + g_yz_0_yyy_xyyy[k];

                g_yz_0_yyyy_xyz[k] = -g_z_0_yyy_xyz[k] - g_yz_0_yyy_xyz[k] * ab_y + g_yz_0_yyy_xyyz[k];

                g_yz_0_yyyy_xzz[k] = -g_z_0_yyy_xzz[k] - g_yz_0_yyy_xzz[k] * ab_y + g_yz_0_yyy_xyzz[k];

                g_yz_0_yyyy_yyy[k] = -g_z_0_yyy_yyy[k] - g_yz_0_yyy_yyy[k] * ab_y + g_yz_0_yyy_yyyy[k];

                g_yz_0_yyyy_yyz[k] = -g_z_0_yyy_yyz[k] - g_yz_0_yyy_yyz[k] * ab_y + g_yz_0_yyy_yyyz[k];

                g_yz_0_yyyy_yzz[k] = -g_z_0_yyy_yzz[k] - g_yz_0_yyy_yzz[k] * ab_y + g_yz_0_yyy_yyzz[k];

                g_yz_0_yyyy_zzz[k] = -g_z_0_yyy_zzz[k] - g_yz_0_yyy_zzz[k] * ab_y + g_yz_0_yyy_yzzz[k];
            }

            /// Set up 710-720 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 710 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 711 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 712 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 713 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 714 * ccomps * dcomps);

            auto g_yz_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 715 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 716 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 717 * ccomps * dcomps);

            auto g_yz_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 718 * ccomps * dcomps);

            auto g_yz_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyz_xxx, g_yz_0_yyyz_xxy, g_yz_0_yyyz_xxz, g_yz_0_yyyz_xyy, g_yz_0_yyyz_xyz, g_yz_0_yyyz_xzz, g_yz_0_yyyz_yyy, g_yz_0_yyyz_yyz, g_yz_0_yyyz_yzz, g_yz_0_yyyz_zzz, g_yz_0_yyz_xxx, g_yz_0_yyz_xxxy, g_yz_0_yyz_xxy, g_yz_0_yyz_xxyy, g_yz_0_yyz_xxyz, g_yz_0_yyz_xxz, g_yz_0_yyz_xyy, g_yz_0_yyz_xyyy, g_yz_0_yyz_xyyz, g_yz_0_yyz_xyz, g_yz_0_yyz_xyzz, g_yz_0_yyz_xzz, g_yz_0_yyz_yyy, g_yz_0_yyz_yyyy, g_yz_0_yyz_yyyz, g_yz_0_yyz_yyz, g_yz_0_yyz_yyzz, g_yz_0_yyz_yzz, g_yz_0_yyz_yzzz, g_yz_0_yyz_zzz, g_z_0_yyz_xxx, g_z_0_yyz_xxy, g_z_0_yyz_xxz, g_z_0_yyz_xyy, g_z_0_yyz_xyz, g_z_0_yyz_xzz, g_z_0_yyz_yyy, g_z_0_yyz_yyz, g_z_0_yyz_yzz, g_z_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyz_xxx[k] = -g_z_0_yyz_xxx[k] - g_yz_0_yyz_xxx[k] * ab_y + g_yz_0_yyz_xxxy[k];

                g_yz_0_yyyz_xxy[k] = -g_z_0_yyz_xxy[k] - g_yz_0_yyz_xxy[k] * ab_y + g_yz_0_yyz_xxyy[k];

                g_yz_0_yyyz_xxz[k] = -g_z_0_yyz_xxz[k] - g_yz_0_yyz_xxz[k] * ab_y + g_yz_0_yyz_xxyz[k];

                g_yz_0_yyyz_xyy[k] = -g_z_0_yyz_xyy[k] - g_yz_0_yyz_xyy[k] * ab_y + g_yz_0_yyz_xyyy[k];

                g_yz_0_yyyz_xyz[k] = -g_z_0_yyz_xyz[k] - g_yz_0_yyz_xyz[k] * ab_y + g_yz_0_yyz_xyyz[k];

                g_yz_0_yyyz_xzz[k] = -g_z_0_yyz_xzz[k] - g_yz_0_yyz_xzz[k] * ab_y + g_yz_0_yyz_xyzz[k];

                g_yz_0_yyyz_yyy[k] = -g_z_0_yyz_yyy[k] - g_yz_0_yyz_yyy[k] * ab_y + g_yz_0_yyz_yyyy[k];

                g_yz_0_yyyz_yyz[k] = -g_z_0_yyz_yyz[k] - g_yz_0_yyz_yyz[k] * ab_y + g_yz_0_yyz_yyyz[k];

                g_yz_0_yyyz_yzz[k] = -g_z_0_yyz_yzz[k] - g_yz_0_yyz_yzz[k] * ab_y + g_yz_0_yyz_yyzz[k];

                g_yz_0_yyyz_zzz[k] = -g_z_0_yyz_zzz[k] - g_yz_0_yyz_zzz[k] * ab_y + g_yz_0_yyz_yzzz[k];
            }

            /// Set up 720-730 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 720 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 721 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 722 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 723 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 724 * ccomps * dcomps);

            auto g_yz_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 725 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 726 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 727 * ccomps * dcomps);

            auto g_yz_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 728 * ccomps * dcomps);

            auto g_yz_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 729 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyzz_xxx, g_yz_0_yyzz_xxy, g_yz_0_yyzz_xxz, g_yz_0_yyzz_xyy, g_yz_0_yyzz_xyz, g_yz_0_yyzz_xzz, g_yz_0_yyzz_yyy, g_yz_0_yyzz_yyz, g_yz_0_yyzz_yzz, g_yz_0_yyzz_zzz, g_yz_0_yzz_xxx, g_yz_0_yzz_xxxy, g_yz_0_yzz_xxy, g_yz_0_yzz_xxyy, g_yz_0_yzz_xxyz, g_yz_0_yzz_xxz, g_yz_0_yzz_xyy, g_yz_0_yzz_xyyy, g_yz_0_yzz_xyyz, g_yz_0_yzz_xyz, g_yz_0_yzz_xyzz, g_yz_0_yzz_xzz, g_yz_0_yzz_yyy, g_yz_0_yzz_yyyy, g_yz_0_yzz_yyyz, g_yz_0_yzz_yyz, g_yz_0_yzz_yyzz, g_yz_0_yzz_yzz, g_yz_0_yzz_yzzz, g_yz_0_yzz_zzz, g_z_0_yzz_xxx, g_z_0_yzz_xxy, g_z_0_yzz_xxz, g_z_0_yzz_xyy, g_z_0_yzz_xyz, g_z_0_yzz_xzz, g_z_0_yzz_yyy, g_z_0_yzz_yyz, g_z_0_yzz_yzz, g_z_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyzz_xxx[k] = -g_z_0_yzz_xxx[k] - g_yz_0_yzz_xxx[k] * ab_y + g_yz_0_yzz_xxxy[k];

                g_yz_0_yyzz_xxy[k] = -g_z_0_yzz_xxy[k] - g_yz_0_yzz_xxy[k] * ab_y + g_yz_0_yzz_xxyy[k];

                g_yz_0_yyzz_xxz[k] = -g_z_0_yzz_xxz[k] - g_yz_0_yzz_xxz[k] * ab_y + g_yz_0_yzz_xxyz[k];

                g_yz_0_yyzz_xyy[k] = -g_z_0_yzz_xyy[k] - g_yz_0_yzz_xyy[k] * ab_y + g_yz_0_yzz_xyyy[k];

                g_yz_0_yyzz_xyz[k] = -g_z_0_yzz_xyz[k] - g_yz_0_yzz_xyz[k] * ab_y + g_yz_0_yzz_xyyz[k];

                g_yz_0_yyzz_xzz[k] = -g_z_0_yzz_xzz[k] - g_yz_0_yzz_xzz[k] * ab_y + g_yz_0_yzz_xyzz[k];

                g_yz_0_yyzz_yyy[k] = -g_z_0_yzz_yyy[k] - g_yz_0_yzz_yyy[k] * ab_y + g_yz_0_yzz_yyyy[k];

                g_yz_0_yyzz_yyz[k] = -g_z_0_yzz_yyz[k] - g_yz_0_yzz_yyz[k] * ab_y + g_yz_0_yzz_yyyz[k];

                g_yz_0_yyzz_yzz[k] = -g_z_0_yzz_yzz[k] - g_yz_0_yzz_yzz[k] * ab_y + g_yz_0_yzz_yyzz[k];

                g_yz_0_yyzz_zzz[k] = -g_z_0_yzz_zzz[k] - g_yz_0_yzz_zzz[k] * ab_y + g_yz_0_yzz_yzzz[k];
            }

            /// Set up 730-740 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 730 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 731 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 732 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 733 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 734 * ccomps * dcomps);

            auto g_yz_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 735 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 736 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 737 * ccomps * dcomps);

            auto g_yz_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 738 * ccomps * dcomps);

            auto g_yz_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 739 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzzz_xxx, g_yz_0_yzzz_xxy, g_yz_0_yzzz_xxz, g_yz_0_yzzz_xyy, g_yz_0_yzzz_xyz, g_yz_0_yzzz_xzz, g_yz_0_yzzz_yyy, g_yz_0_yzzz_yyz, g_yz_0_yzzz_yzz, g_yz_0_yzzz_zzz, g_yz_0_zzz_xxx, g_yz_0_zzz_xxxy, g_yz_0_zzz_xxy, g_yz_0_zzz_xxyy, g_yz_0_zzz_xxyz, g_yz_0_zzz_xxz, g_yz_0_zzz_xyy, g_yz_0_zzz_xyyy, g_yz_0_zzz_xyyz, g_yz_0_zzz_xyz, g_yz_0_zzz_xyzz, g_yz_0_zzz_xzz, g_yz_0_zzz_yyy, g_yz_0_zzz_yyyy, g_yz_0_zzz_yyyz, g_yz_0_zzz_yyz, g_yz_0_zzz_yyzz, g_yz_0_zzz_yzz, g_yz_0_zzz_yzzz, g_yz_0_zzz_zzz, g_z_0_zzz_xxx, g_z_0_zzz_xxy, g_z_0_zzz_xxz, g_z_0_zzz_xyy, g_z_0_zzz_xyz, g_z_0_zzz_xzz, g_z_0_zzz_yyy, g_z_0_zzz_yyz, g_z_0_zzz_yzz, g_z_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzzz_xxx[k] = -g_z_0_zzz_xxx[k] - g_yz_0_zzz_xxx[k] * ab_y + g_yz_0_zzz_xxxy[k];

                g_yz_0_yzzz_xxy[k] = -g_z_0_zzz_xxy[k] - g_yz_0_zzz_xxy[k] * ab_y + g_yz_0_zzz_xxyy[k];

                g_yz_0_yzzz_xxz[k] = -g_z_0_zzz_xxz[k] - g_yz_0_zzz_xxz[k] * ab_y + g_yz_0_zzz_xxyz[k];

                g_yz_0_yzzz_xyy[k] = -g_z_0_zzz_xyy[k] - g_yz_0_zzz_xyy[k] * ab_y + g_yz_0_zzz_xyyy[k];

                g_yz_0_yzzz_xyz[k] = -g_z_0_zzz_xyz[k] - g_yz_0_zzz_xyz[k] * ab_y + g_yz_0_zzz_xyyz[k];

                g_yz_0_yzzz_xzz[k] = -g_z_0_zzz_xzz[k] - g_yz_0_zzz_xzz[k] * ab_y + g_yz_0_zzz_xyzz[k];

                g_yz_0_yzzz_yyy[k] = -g_z_0_zzz_yyy[k] - g_yz_0_zzz_yyy[k] * ab_y + g_yz_0_zzz_yyyy[k];

                g_yz_0_yzzz_yyz[k] = -g_z_0_zzz_yyz[k] - g_yz_0_zzz_yyz[k] * ab_y + g_yz_0_zzz_yyyz[k];

                g_yz_0_yzzz_yzz[k] = -g_z_0_zzz_yzz[k] - g_yz_0_zzz_yzz[k] * ab_y + g_yz_0_zzz_yyzz[k];

                g_yz_0_yzzz_zzz[k] = -g_z_0_zzz_zzz[k] - g_yz_0_zzz_zzz[k] * ab_y + g_yz_0_zzz_yzzz[k];
            }

            /// Set up 740-750 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 740 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 741 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 742 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 743 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 744 * ccomps * dcomps);

            auto g_yz_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 745 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 746 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 747 * ccomps * dcomps);

            auto g_yz_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 748 * ccomps * dcomps);

            auto g_yz_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzz_xxx, g_y_0_zzz_xxy, g_y_0_zzz_xxz, g_y_0_zzz_xyy, g_y_0_zzz_xyz, g_y_0_zzz_xzz, g_y_0_zzz_yyy, g_y_0_zzz_yyz, g_y_0_zzz_yzz, g_y_0_zzz_zzz, g_yz_0_zzz_xxx, g_yz_0_zzz_xxxz, g_yz_0_zzz_xxy, g_yz_0_zzz_xxyz, g_yz_0_zzz_xxz, g_yz_0_zzz_xxzz, g_yz_0_zzz_xyy, g_yz_0_zzz_xyyz, g_yz_0_zzz_xyz, g_yz_0_zzz_xyzz, g_yz_0_zzz_xzz, g_yz_0_zzz_xzzz, g_yz_0_zzz_yyy, g_yz_0_zzz_yyyz, g_yz_0_zzz_yyz, g_yz_0_zzz_yyzz, g_yz_0_zzz_yzz, g_yz_0_zzz_yzzz, g_yz_0_zzz_zzz, g_yz_0_zzz_zzzz, g_yz_0_zzzz_xxx, g_yz_0_zzzz_xxy, g_yz_0_zzzz_xxz, g_yz_0_zzzz_xyy, g_yz_0_zzzz_xyz, g_yz_0_zzzz_xzz, g_yz_0_zzzz_yyy, g_yz_0_zzzz_yyz, g_yz_0_zzzz_yzz, g_yz_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzzz_xxx[k] = -g_y_0_zzz_xxx[k] - g_yz_0_zzz_xxx[k] * ab_z + g_yz_0_zzz_xxxz[k];

                g_yz_0_zzzz_xxy[k] = -g_y_0_zzz_xxy[k] - g_yz_0_zzz_xxy[k] * ab_z + g_yz_0_zzz_xxyz[k];

                g_yz_0_zzzz_xxz[k] = -g_y_0_zzz_xxz[k] - g_yz_0_zzz_xxz[k] * ab_z + g_yz_0_zzz_xxzz[k];

                g_yz_0_zzzz_xyy[k] = -g_y_0_zzz_xyy[k] - g_yz_0_zzz_xyy[k] * ab_z + g_yz_0_zzz_xyyz[k];

                g_yz_0_zzzz_xyz[k] = -g_y_0_zzz_xyz[k] - g_yz_0_zzz_xyz[k] * ab_z + g_yz_0_zzz_xyzz[k];

                g_yz_0_zzzz_xzz[k] = -g_y_0_zzz_xzz[k] - g_yz_0_zzz_xzz[k] * ab_z + g_yz_0_zzz_xzzz[k];

                g_yz_0_zzzz_yyy[k] = -g_y_0_zzz_yyy[k] - g_yz_0_zzz_yyy[k] * ab_z + g_yz_0_zzz_yyyz[k];

                g_yz_0_zzzz_yyz[k] = -g_y_0_zzz_yyz[k] - g_yz_0_zzz_yyz[k] * ab_z + g_yz_0_zzz_yyzz[k];

                g_yz_0_zzzz_yzz[k] = -g_y_0_zzz_yzz[k] - g_yz_0_zzz_yzz[k] * ab_z + g_yz_0_zzz_yzzz[k];

                g_yz_0_zzzz_zzz[k] = -g_y_0_zzz_zzz[k] - g_yz_0_zzz_zzz[k] * ab_z + g_yz_0_zzz_zzzz[k];
            }

            /// Set up 750-760 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 750 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 751 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 752 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 753 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 754 * ccomps * dcomps);

            auto g_zz_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 755 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 756 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 757 * ccomps * dcomps);

            auto g_zz_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 758 * ccomps * dcomps);

            auto g_zz_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 759 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxx_xxx, g_zz_0_xxx_xxxx, g_zz_0_xxx_xxxy, g_zz_0_xxx_xxxz, g_zz_0_xxx_xxy, g_zz_0_xxx_xxyy, g_zz_0_xxx_xxyz, g_zz_0_xxx_xxz, g_zz_0_xxx_xxzz, g_zz_0_xxx_xyy, g_zz_0_xxx_xyyy, g_zz_0_xxx_xyyz, g_zz_0_xxx_xyz, g_zz_0_xxx_xyzz, g_zz_0_xxx_xzz, g_zz_0_xxx_xzzz, g_zz_0_xxx_yyy, g_zz_0_xxx_yyz, g_zz_0_xxx_yzz, g_zz_0_xxx_zzz, g_zz_0_xxxx_xxx, g_zz_0_xxxx_xxy, g_zz_0_xxxx_xxz, g_zz_0_xxxx_xyy, g_zz_0_xxxx_xyz, g_zz_0_xxxx_xzz, g_zz_0_xxxx_yyy, g_zz_0_xxxx_yyz, g_zz_0_xxxx_yzz, g_zz_0_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxx_xxx[k] = -g_zz_0_xxx_xxx[k] * ab_x + g_zz_0_xxx_xxxx[k];

                g_zz_0_xxxx_xxy[k] = -g_zz_0_xxx_xxy[k] * ab_x + g_zz_0_xxx_xxxy[k];

                g_zz_0_xxxx_xxz[k] = -g_zz_0_xxx_xxz[k] * ab_x + g_zz_0_xxx_xxxz[k];

                g_zz_0_xxxx_xyy[k] = -g_zz_0_xxx_xyy[k] * ab_x + g_zz_0_xxx_xxyy[k];

                g_zz_0_xxxx_xyz[k] = -g_zz_0_xxx_xyz[k] * ab_x + g_zz_0_xxx_xxyz[k];

                g_zz_0_xxxx_xzz[k] = -g_zz_0_xxx_xzz[k] * ab_x + g_zz_0_xxx_xxzz[k];

                g_zz_0_xxxx_yyy[k] = -g_zz_0_xxx_yyy[k] * ab_x + g_zz_0_xxx_xyyy[k];

                g_zz_0_xxxx_yyz[k] = -g_zz_0_xxx_yyz[k] * ab_x + g_zz_0_xxx_xyyz[k];

                g_zz_0_xxxx_yzz[k] = -g_zz_0_xxx_yzz[k] * ab_x + g_zz_0_xxx_xyzz[k];

                g_zz_0_xxxx_zzz[k] = -g_zz_0_xxx_zzz[k] * ab_x + g_zz_0_xxx_xzzz[k];
            }

            /// Set up 760-770 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 760 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 761 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 762 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 763 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 764 * ccomps * dcomps);

            auto g_zz_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 765 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 766 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 767 * ccomps * dcomps);

            auto g_zz_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 768 * ccomps * dcomps);

            auto g_zz_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 769 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxy_xxx, g_zz_0_xxxy_xxy, g_zz_0_xxxy_xxz, g_zz_0_xxxy_xyy, g_zz_0_xxxy_xyz, g_zz_0_xxxy_xzz, g_zz_0_xxxy_yyy, g_zz_0_xxxy_yyz, g_zz_0_xxxy_yzz, g_zz_0_xxxy_zzz, g_zz_0_xxy_xxx, g_zz_0_xxy_xxxx, g_zz_0_xxy_xxxy, g_zz_0_xxy_xxxz, g_zz_0_xxy_xxy, g_zz_0_xxy_xxyy, g_zz_0_xxy_xxyz, g_zz_0_xxy_xxz, g_zz_0_xxy_xxzz, g_zz_0_xxy_xyy, g_zz_0_xxy_xyyy, g_zz_0_xxy_xyyz, g_zz_0_xxy_xyz, g_zz_0_xxy_xyzz, g_zz_0_xxy_xzz, g_zz_0_xxy_xzzz, g_zz_0_xxy_yyy, g_zz_0_xxy_yyz, g_zz_0_xxy_yzz, g_zz_0_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxy_xxx[k] = -g_zz_0_xxy_xxx[k] * ab_x + g_zz_0_xxy_xxxx[k];

                g_zz_0_xxxy_xxy[k] = -g_zz_0_xxy_xxy[k] * ab_x + g_zz_0_xxy_xxxy[k];

                g_zz_0_xxxy_xxz[k] = -g_zz_0_xxy_xxz[k] * ab_x + g_zz_0_xxy_xxxz[k];

                g_zz_0_xxxy_xyy[k] = -g_zz_0_xxy_xyy[k] * ab_x + g_zz_0_xxy_xxyy[k];

                g_zz_0_xxxy_xyz[k] = -g_zz_0_xxy_xyz[k] * ab_x + g_zz_0_xxy_xxyz[k];

                g_zz_0_xxxy_xzz[k] = -g_zz_0_xxy_xzz[k] * ab_x + g_zz_0_xxy_xxzz[k];

                g_zz_0_xxxy_yyy[k] = -g_zz_0_xxy_yyy[k] * ab_x + g_zz_0_xxy_xyyy[k];

                g_zz_0_xxxy_yyz[k] = -g_zz_0_xxy_yyz[k] * ab_x + g_zz_0_xxy_xyyz[k];

                g_zz_0_xxxy_yzz[k] = -g_zz_0_xxy_yzz[k] * ab_x + g_zz_0_xxy_xyzz[k];

                g_zz_0_xxxy_zzz[k] = -g_zz_0_xxy_zzz[k] * ab_x + g_zz_0_xxy_xzzz[k];
            }

            /// Set up 770-780 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 770 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 771 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 772 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 773 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 774 * ccomps * dcomps);

            auto g_zz_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 775 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 776 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 777 * ccomps * dcomps);

            auto g_zz_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 778 * ccomps * dcomps);

            auto g_zz_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxz_xxx, g_zz_0_xxxz_xxy, g_zz_0_xxxz_xxz, g_zz_0_xxxz_xyy, g_zz_0_xxxz_xyz, g_zz_0_xxxz_xzz, g_zz_0_xxxz_yyy, g_zz_0_xxxz_yyz, g_zz_0_xxxz_yzz, g_zz_0_xxxz_zzz, g_zz_0_xxz_xxx, g_zz_0_xxz_xxxx, g_zz_0_xxz_xxxy, g_zz_0_xxz_xxxz, g_zz_0_xxz_xxy, g_zz_0_xxz_xxyy, g_zz_0_xxz_xxyz, g_zz_0_xxz_xxz, g_zz_0_xxz_xxzz, g_zz_0_xxz_xyy, g_zz_0_xxz_xyyy, g_zz_0_xxz_xyyz, g_zz_0_xxz_xyz, g_zz_0_xxz_xyzz, g_zz_0_xxz_xzz, g_zz_0_xxz_xzzz, g_zz_0_xxz_yyy, g_zz_0_xxz_yyz, g_zz_0_xxz_yzz, g_zz_0_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxz_xxx[k] = -g_zz_0_xxz_xxx[k] * ab_x + g_zz_0_xxz_xxxx[k];

                g_zz_0_xxxz_xxy[k] = -g_zz_0_xxz_xxy[k] * ab_x + g_zz_0_xxz_xxxy[k];

                g_zz_0_xxxz_xxz[k] = -g_zz_0_xxz_xxz[k] * ab_x + g_zz_0_xxz_xxxz[k];

                g_zz_0_xxxz_xyy[k] = -g_zz_0_xxz_xyy[k] * ab_x + g_zz_0_xxz_xxyy[k];

                g_zz_0_xxxz_xyz[k] = -g_zz_0_xxz_xyz[k] * ab_x + g_zz_0_xxz_xxyz[k];

                g_zz_0_xxxz_xzz[k] = -g_zz_0_xxz_xzz[k] * ab_x + g_zz_0_xxz_xxzz[k];

                g_zz_0_xxxz_yyy[k] = -g_zz_0_xxz_yyy[k] * ab_x + g_zz_0_xxz_xyyy[k];

                g_zz_0_xxxz_yyz[k] = -g_zz_0_xxz_yyz[k] * ab_x + g_zz_0_xxz_xyyz[k];

                g_zz_0_xxxz_yzz[k] = -g_zz_0_xxz_yzz[k] * ab_x + g_zz_0_xxz_xyzz[k];

                g_zz_0_xxxz_zzz[k] = -g_zz_0_xxz_zzz[k] * ab_x + g_zz_0_xxz_xzzz[k];
            }

            /// Set up 780-790 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 780 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 781 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 782 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 783 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 784 * ccomps * dcomps);

            auto g_zz_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 785 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 786 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 787 * ccomps * dcomps);

            auto g_zz_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 788 * ccomps * dcomps);

            auto g_zz_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 789 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyy_xxx, g_zz_0_xxyy_xxy, g_zz_0_xxyy_xxz, g_zz_0_xxyy_xyy, g_zz_0_xxyy_xyz, g_zz_0_xxyy_xzz, g_zz_0_xxyy_yyy, g_zz_0_xxyy_yyz, g_zz_0_xxyy_yzz, g_zz_0_xxyy_zzz, g_zz_0_xyy_xxx, g_zz_0_xyy_xxxx, g_zz_0_xyy_xxxy, g_zz_0_xyy_xxxz, g_zz_0_xyy_xxy, g_zz_0_xyy_xxyy, g_zz_0_xyy_xxyz, g_zz_0_xyy_xxz, g_zz_0_xyy_xxzz, g_zz_0_xyy_xyy, g_zz_0_xyy_xyyy, g_zz_0_xyy_xyyz, g_zz_0_xyy_xyz, g_zz_0_xyy_xyzz, g_zz_0_xyy_xzz, g_zz_0_xyy_xzzz, g_zz_0_xyy_yyy, g_zz_0_xyy_yyz, g_zz_0_xyy_yzz, g_zz_0_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyy_xxx[k] = -g_zz_0_xyy_xxx[k] * ab_x + g_zz_0_xyy_xxxx[k];

                g_zz_0_xxyy_xxy[k] = -g_zz_0_xyy_xxy[k] * ab_x + g_zz_0_xyy_xxxy[k];

                g_zz_0_xxyy_xxz[k] = -g_zz_0_xyy_xxz[k] * ab_x + g_zz_0_xyy_xxxz[k];

                g_zz_0_xxyy_xyy[k] = -g_zz_0_xyy_xyy[k] * ab_x + g_zz_0_xyy_xxyy[k];

                g_zz_0_xxyy_xyz[k] = -g_zz_0_xyy_xyz[k] * ab_x + g_zz_0_xyy_xxyz[k];

                g_zz_0_xxyy_xzz[k] = -g_zz_0_xyy_xzz[k] * ab_x + g_zz_0_xyy_xxzz[k];

                g_zz_0_xxyy_yyy[k] = -g_zz_0_xyy_yyy[k] * ab_x + g_zz_0_xyy_xyyy[k];

                g_zz_0_xxyy_yyz[k] = -g_zz_0_xyy_yyz[k] * ab_x + g_zz_0_xyy_xyyz[k];

                g_zz_0_xxyy_yzz[k] = -g_zz_0_xyy_yzz[k] * ab_x + g_zz_0_xyy_xyzz[k];

                g_zz_0_xxyy_zzz[k] = -g_zz_0_xyy_zzz[k] * ab_x + g_zz_0_xyy_xzzz[k];
            }

            /// Set up 790-800 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 790 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 791 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 792 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 793 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 794 * ccomps * dcomps);

            auto g_zz_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 795 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 796 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 797 * ccomps * dcomps);

            auto g_zz_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 798 * ccomps * dcomps);

            auto g_zz_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 799 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyz_xxx, g_zz_0_xxyz_xxy, g_zz_0_xxyz_xxz, g_zz_0_xxyz_xyy, g_zz_0_xxyz_xyz, g_zz_0_xxyz_xzz, g_zz_0_xxyz_yyy, g_zz_0_xxyz_yyz, g_zz_0_xxyz_yzz, g_zz_0_xxyz_zzz, g_zz_0_xyz_xxx, g_zz_0_xyz_xxxx, g_zz_0_xyz_xxxy, g_zz_0_xyz_xxxz, g_zz_0_xyz_xxy, g_zz_0_xyz_xxyy, g_zz_0_xyz_xxyz, g_zz_0_xyz_xxz, g_zz_0_xyz_xxzz, g_zz_0_xyz_xyy, g_zz_0_xyz_xyyy, g_zz_0_xyz_xyyz, g_zz_0_xyz_xyz, g_zz_0_xyz_xyzz, g_zz_0_xyz_xzz, g_zz_0_xyz_xzzz, g_zz_0_xyz_yyy, g_zz_0_xyz_yyz, g_zz_0_xyz_yzz, g_zz_0_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyz_xxx[k] = -g_zz_0_xyz_xxx[k] * ab_x + g_zz_0_xyz_xxxx[k];

                g_zz_0_xxyz_xxy[k] = -g_zz_0_xyz_xxy[k] * ab_x + g_zz_0_xyz_xxxy[k];

                g_zz_0_xxyz_xxz[k] = -g_zz_0_xyz_xxz[k] * ab_x + g_zz_0_xyz_xxxz[k];

                g_zz_0_xxyz_xyy[k] = -g_zz_0_xyz_xyy[k] * ab_x + g_zz_0_xyz_xxyy[k];

                g_zz_0_xxyz_xyz[k] = -g_zz_0_xyz_xyz[k] * ab_x + g_zz_0_xyz_xxyz[k];

                g_zz_0_xxyz_xzz[k] = -g_zz_0_xyz_xzz[k] * ab_x + g_zz_0_xyz_xxzz[k];

                g_zz_0_xxyz_yyy[k] = -g_zz_0_xyz_yyy[k] * ab_x + g_zz_0_xyz_xyyy[k];

                g_zz_0_xxyz_yyz[k] = -g_zz_0_xyz_yyz[k] * ab_x + g_zz_0_xyz_xyyz[k];

                g_zz_0_xxyz_yzz[k] = -g_zz_0_xyz_yzz[k] * ab_x + g_zz_0_xyz_xyzz[k];

                g_zz_0_xxyz_zzz[k] = -g_zz_0_xyz_zzz[k] * ab_x + g_zz_0_xyz_xzzz[k];
            }

            /// Set up 800-810 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 800 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 801 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 802 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 803 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 804 * ccomps * dcomps);

            auto g_zz_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 805 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 806 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 807 * ccomps * dcomps);

            auto g_zz_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 808 * ccomps * dcomps);

            auto g_zz_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxzz_xxx, g_zz_0_xxzz_xxy, g_zz_0_xxzz_xxz, g_zz_0_xxzz_xyy, g_zz_0_xxzz_xyz, g_zz_0_xxzz_xzz, g_zz_0_xxzz_yyy, g_zz_0_xxzz_yyz, g_zz_0_xxzz_yzz, g_zz_0_xxzz_zzz, g_zz_0_xzz_xxx, g_zz_0_xzz_xxxx, g_zz_0_xzz_xxxy, g_zz_0_xzz_xxxz, g_zz_0_xzz_xxy, g_zz_0_xzz_xxyy, g_zz_0_xzz_xxyz, g_zz_0_xzz_xxz, g_zz_0_xzz_xxzz, g_zz_0_xzz_xyy, g_zz_0_xzz_xyyy, g_zz_0_xzz_xyyz, g_zz_0_xzz_xyz, g_zz_0_xzz_xyzz, g_zz_0_xzz_xzz, g_zz_0_xzz_xzzz, g_zz_0_xzz_yyy, g_zz_0_xzz_yyz, g_zz_0_xzz_yzz, g_zz_0_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxzz_xxx[k] = -g_zz_0_xzz_xxx[k] * ab_x + g_zz_0_xzz_xxxx[k];

                g_zz_0_xxzz_xxy[k] = -g_zz_0_xzz_xxy[k] * ab_x + g_zz_0_xzz_xxxy[k];

                g_zz_0_xxzz_xxz[k] = -g_zz_0_xzz_xxz[k] * ab_x + g_zz_0_xzz_xxxz[k];

                g_zz_0_xxzz_xyy[k] = -g_zz_0_xzz_xyy[k] * ab_x + g_zz_0_xzz_xxyy[k];

                g_zz_0_xxzz_xyz[k] = -g_zz_0_xzz_xyz[k] * ab_x + g_zz_0_xzz_xxyz[k];

                g_zz_0_xxzz_xzz[k] = -g_zz_0_xzz_xzz[k] * ab_x + g_zz_0_xzz_xxzz[k];

                g_zz_0_xxzz_yyy[k] = -g_zz_0_xzz_yyy[k] * ab_x + g_zz_0_xzz_xyyy[k];

                g_zz_0_xxzz_yyz[k] = -g_zz_0_xzz_yyz[k] * ab_x + g_zz_0_xzz_xyyz[k];

                g_zz_0_xxzz_yzz[k] = -g_zz_0_xzz_yzz[k] * ab_x + g_zz_0_xzz_xyzz[k];

                g_zz_0_xxzz_zzz[k] = -g_zz_0_xzz_zzz[k] * ab_x + g_zz_0_xzz_xzzz[k];
            }

            /// Set up 810-820 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 810 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 811 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 812 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 813 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 814 * ccomps * dcomps);

            auto g_zz_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 815 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 816 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 817 * ccomps * dcomps);

            auto g_zz_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 818 * ccomps * dcomps);

            auto g_zz_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 819 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyy_xxx, g_zz_0_xyyy_xxy, g_zz_0_xyyy_xxz, g_zz_0_xyyy_xyy, g_zz_0_xyyy_xyz, g_zz_0_xyyy_xzz, g_zz_0_xyyy_yyy, g_zz_0_xyyy_yyz, g_zz_0_xyyy_yzz, g_zz_0_xyyy_zzz, g_zz_0_yyy_xxx, g_zz_0_yyy_xxxx, g_zz_0_yyy_xxxy, g_zz_0_yyy_xxxz, g_zz_0_yyy_xxy, g_zz_0_yyy_xxyy, g_zz_0_yyy_xxyz, g_zz_0_yyy_xxz, g_zz_0_yyy_xxzz, g_zz_0_yyy_xyy, g_zz_0_yyy_xyyy, g_zz_0_yyy_xyyz, g_zz_0_yyy_xyz, g_zz_0_yyy_xyzz, g_zz_0_yyy_xzz, g_zz_0_yyy_xzzz, g_zz_0_yyy_yyy, g_zz_0_yyy_yyz, g_zz_0_yyy_yzz, g_zz_0_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyy_xxx[k] = -g_zz_0_yyy_xxx[k] * ab_x + g_zz_0_yyy_xxxx[k];

                g_zz_0_xyyy_xxy[k] = -g_zz_0_yyy_xxy[k] * ab_x + g_zz_0_yyy_xxxy[k];

                g_zz_0_xyyy_xxz[k] = -g_zz_0_yyy_xxz[k] * ab_x + g_zz_0_yyy_xxxz[k];

                g_zz_0_xyyy_xyy[k] = -g_zz_0_yyy_xyy[k] * ab_x + g_zz_0_yyy_xxyy[k];

                g_zz_0_xyyy_xyz[k] = -g_zz_0_yyy_xyz[k] * ab_x + g_zz_0_yyy_xxyz[k];

                g_zz_0_xyyy_xzz[k] = -g_zz_0_yyy_xzz[k] * ab_x + g_zz_0_yyy_xxzz[k];

                g_zz_0_xyyy_yyy[k] = -g_zz_0_yyy_yyy[k] * ab_x + g_zz_0_yyy_xyyy[k];

                g_zz_0_xyyy_yyz[k] = -g_zz_0_yyy_yyz[k] * ab_x + g_zz_0_yyy_xyyz[k];

                g_zz_0_xyyy_yzz[k] = -g_zz_0_yyy_yzz[k] * ab_x + g_zz_0_yyy_xyzz[k];

                g_zz_0_xyyy_zzz[k] = -g_zz_0_yyy_zzz[k] * ab_x + g_zz_0_yyy_xzzz[k];
            }

            /// Set up 820-830 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 820 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 821 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 822 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 823 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 824 * ccomps * dcomps);

            auto g_zz_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 825 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 826 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 827 * ccomps * dcomps);

            auto g_zz_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 828 * ccomps * dcomps);

            auto g_zz_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 829 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyz_xxx, g_zz_0_xyyz_xxy, g_zz_0_xyyz_xxz, g_zz_0_xyyz_xyy, g_zz_0_xyyz_xyz, g_zz_0_xyyz_xzz, g_zz_0_xyyz_yyy, g_zz_0_xyyz_yyz, g_zz_0_xyyz_yzz, g_zz_0_xyyz_zzz, g_zz_0_yyz_xxx, g_zz_0_yyz_xxxx, g_zz_0_yyz_xxxy, g_zz_0_yyz_xxxz, g_zz_0_yyz_xxy, g_zz_0_yyz_xxyy, g_zz_0_yyz_xxyz, g_zz_0_yyz_xxz, g_zz_0_yyz_xxzz, g_zz_0_yyz_xyy, g_zz_0_yyz_xyyy, g_zz_0_yyz_xyyz, g_zz_0_yyz_xyz, g_zz_0_yyz_xyzz, g_zz_0_yyz_xzz, g_zz_0_yyz_xzzz, g_zz_0_yyz_yyy, g_zz_0_yyz_yyz, g_zz_0_yyz_yzz, g_zz_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyz_xxx[k] = -g_zz_0_yyz_xxx[k] * ab_x + g_zz_0_yyz_xxxx[k];

                g_zz_0_xyyz_xxy[k] = -g_zz_0_yyz_xxy[k] * ab_x + g_zz_0_yyz_xxxy[k];

                g_zz_0_xyyz_xxz[k] = -g_zz_0_yyz_xxz[k] * ab_x + g_zz_0_yyz_xxxz[k];

                g_zz_0_xyyz_xyy[k] = -g_zz_0_yyz_xyy[k] * ab_x + g_zz_0_yyz_xxyy[k];

                g_zz_0_xyyz_xyz[k] = -g_zz_0_yyz_xyz[k] * ab_x + g_zz_0_yyz_xxyz[k];

                g_zz_0_xyyz_xzz[k] = -g_zz_0_yyz_xzz[k] * ab_x + g_zz_0_yyz_xxzz[k];

                g_zz_0_xyyz_yyy[k] = -g_zz_0_yyz_yyy[k] * ab_x + g_zz_0_yyz_xyyy[k];

                g_zz_0_xyyz_yyz[k] = -g_zz_0_yyz_yyz[k] * ab_x + g_zz_0_yyz_xyyz[k];

                g_zz_0_xyyz_yzz[k] = -g_zz_0_yyz_yzz[k] * ab_x + g_zz_0_yyz_xyzz[k];

                g_zz_0_xyyz_zzz[k] = -g_zz_0_yyz_zzz[k] * ab_x + g_zz_0_yyz_xzzz[k];
            }

            /// Set up 830-840 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 830 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 831 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 832 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 833 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 834 * ccomps * dcomps);

            auto g_zz_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 835 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 836 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 837 * ccomps * dcomps);

            auto g_zz_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 838 * ccomps * dcomps);

            auto g_zz_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyzz_xxx, g_zz_0_xyzz_xxy, g_zz_0_xyzz_xxz, g_zz_0_xyzz_xyy, g_zz_0_xyzz_xyz, g_zz_0_xyzz_xzz, g_zz_0_xyzz_yyy, g_zz_0_xyzz_yyz, g_zz_0_xyzz_yzz, g_zz_0_xyzz_zzz, g_zz_0_yzz_xxx, g_zz_0_yzz_xxxx, g_zz_0_yzz_xxxy, g_zz_0_yzz_xxxz, g_zz_0_yzz_xxy, g_zz_0_yzz_xxyy, g_zz_0_yzz_xxyz, g_zz_0_yzz_xxz, g_zz_0_yzz_xxzz, g_zz_0_yzz_xyy, g_zz_0_yzz_xyyy, g_zz_0_yzz_xyyz, g_zz_0_yzz_xyz, g_zz_0_yzz_xyzz, g_zz_0_yzz_xzz, g_zz_0_yzz_xzzz, g_zz_0_yzz_yyy, g_zz_0_yzz_yyz, g_zz_0_yzz_yzz, g_zz_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyzz_xxx[k] = -g_zz_0_yzz_xxx[k] * ab_x + g_zz_0_yzz_xxxx[k];

                g_zz_0_xyzz_xxy[k] = -g_zz_0_yzz_xxy[k] * ab_x + g_zz_0_yzz_xxxy[k];

                g_zz_0_xyzz_xxz[k] = -g_zz_0_yzz_xxz[k] * ab_x + g_zz_0_yzz_xxxz[k];

                g_zz_0_xyzz_xyy[k] = -g_zz_0_yzz_xyy[k] * ab_x + g_zz_0_yzz_xxyy[k];

                g_zz_0_xyzz_xyz[k] = -g_zz_0_yzz_xyz[k] * ab_x + g_zz_0_yzz_xxyz[k];

                g_zz_0_xyzz_xzz[k] = -g_zz_0_yzz_xzz[k] * ab_x + g_zz_0_yzz_xxzz[k];

                g_zz_0_xyzz_yyy[k] = -g_zz_0_yzz_yyy[k] * ab_x + g_zz_0_yzz_xyyy[k];

                g_zz_0_xyzz_yyz[k] = -g_zz_0_yzz_yyz[k] * ab_x + g_zz_0_yzz_xyyz[k];

                g_zz_0_xyzz_yzz[k] = -g_zz_0_yzz_yzz[k] * ab_x + g_zz_0_yzz_xyzz[k];

                g_zz_0_xyzz_zzz[k] = -g_zz_0_yzz_zzz[k] * ab_x + g_zz_0_yzz_xzzz[k];
            }

            /// Set up 840-850 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 840 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 841 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 842 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 843 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 844 * ccomps * dcomps);

            auto g_zz_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 845 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 846 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 847 * ccomps * dcomps);

            auto g_zz_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 848 * ccomps * dcomps);

            auto g_zz_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 849 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzzz_xxx, g_zz_0_xzzz_xxy, g_zz_0_xzzz_xxz, g_zz_0_xzzz_xyy, g_zz_0_xzzz_xyz, g_zz_0_xzzz_xzz, g_zz_0_xzzz_yyy, g_zz_0_xzzz_yyz, g_zz_0_xzzz_yzz, g_zz_0_xzzz_zzz, g_zz_0_zzz_xxx, g_zz_0_zzz_xxxx, g_zz_0_zzz_xxxy, g_zz_0_zzz_xxxz, g_zz_0_zzz_xxy, g_zz_0_zzz_xxyy, g_zz_0_zzz_xxyz, g_zz_0_zzz_xxz, g_zz_0_zzz_xxzz, g_zz_0_zzz_xyy, g_zz_0_zzz_xyyy, g_zz_0_zzz_xyyz, g_zz_0_zzz_xyz, g_zz_0_zzz_xyzz, g_zz_0_zzz_xzz, g_zz_0_zzz_xzzz, g_zz_0_zzz_yyy, g_zz_0_zzz_yyz, g_zz_0_zzz_yzz, g_zz_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzzz_xxx[k] = -g_zz_0_zzz_xxx[k] * ab_x + g_zz_0_zzz_xxxx[k];

                g_zz_0_xzzz_xxy[k] = -g_zz_0_zzz_xxy[k] * ab_x + g_zz_0_zzz_xxxy[k];

                g_zz_0_xzzz_xxz[k] = -g_zz_0_zzz_xxz[k] * ab_x + g_zz_0_zzz_xxxz[k];

                g_zz_0_xzzz_xyy[k] = -g_zz_0_zzz_xyy[k] * ab_x + g_zz_0_zzz_xxyy[k];

                g_zz_0_xzzz_xyz[k] = -g_zz_0_zzz_xyz[k] * ab_x + g_zz_0_zzz_xxyz[k];

                g_zz_0_xzzz_xzz[k] = -g_zz_0_zzz_xzz[k] * ab_x + g_zz_0_zzz_xxzz[k];

                g_zz_0_xzzz_yyy[k] = -g_zz_0_zzz_yyy[k] * ab_x + g_zz_0_zzz_xyyy[k];

                g_zz_0_xzzz_yyz[k] = -g_zz_0_zzz_yyz[k] * ab_x + g_zz_0_zzz_xyyz[k];

                g_zz_0_xzzz_yzz[k] = -g_zz_0_zzz_yzz[k] * ab_x + g_zz_0_zzz_xyzz[k];

                g_zz_0_xzzz_zzz[k] = -g_zz_0_zzz_zzz[k] * ab_x + g_zz_0_zzz_xzzz[k];
            }

            /// Set up 850-860 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 850 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 851 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 852 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 853 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 854 * ccomps * dcomps);

            auto g_zz_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 855 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 856 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 857 * ccomps * dcomps);

            auto g_zz_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 858 * ccomps * dcomps);

            auto g_zz_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 859 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyy_xxx, g_zz_0_yyy_xxxy, g_zz_0_yyy_xxy, g_zz_0_yyy_xxyy, g_zz_0_yyy_xxyz, g_zz_0_yyy_xxz, g_zz_0_yyy_xyy, g_zz_0_yyy_xyyy, g_zz_0_yyy_xyyz, g_zz_0_yyy_xyz, g_zz_0_yyy_xyzz, g_zz_0_yyy_xzz, g_zz_0_yyy_yyy, g_zz_0_yyy_yyyy, g_zz_0_yyy_yyyz, g_zz_0_yyy_yyz, g_zz_0_yyy_yyzz, g_zz_0_yyy_yzz, g_zz_0_yyy_yzzz, g_zz_0_yyy_zzz, g_zz_0_yyyy_xxx, g_zz_0_yyyy_xxy, g_zz_0_yyyy_xxz, g_zz_0_yyyy_xyy, g_zz_0_yyyy_xyz, g_zz_0_yyyy_xzz, g_zz_0_yyyy_yyy, g_zz_0_yyyy_yyz, g_zz_0_yyyy_yzz, g_zz_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyy_xxx[k] = -g_zz_0_yyy_xxx[k] * ab_y + g_zz_0_yyy_xxxy[k];

                g_zz_0_yyyy_xxy[k] = -g_zz_0_yyy_xxy[k] * ab_y + g_zz_0_yyy_xxyy[k];

                g_zz_0_yyyy_xxz[k] = -g_zz_0_yyy_xxz[k] * ab_y + g_zz_0_yyy_xxyz[k];

                g_zz_0_yyyy_xyy[k] = -g_zz_0_yyy_xyy[k] * ab_y + g_zz_0_yyy_xyyy[k];

                g_zz_0_yyyy_xyz[k] = -g_zz_0_yyy_xyz[k] * ab_y + g_zz_0_yyy_xyyz[k];

                g_zz_0_yyyy_xzz[k] = -g_zz_0_yyy_xzz[k] * ab_y + g_zz_0_yyy_xyzz[k];

                g_zz_0_yyyy_yyy[k] = -g_zz_0_yyy_yyy[k] * ab_y + g_zz_0_yyy_yyyy[k];

                g_zz_0_yyyy_yyz[k] = -g_zz_0_yyy_yyz[k] * ab_y + g_zz_0_yyy_yyyz[k];

                g_zz_0_yyyy_yzz[k] = -g_zz_0_yyy_yzz[k] * ab_y + g_zz_0_yyy_yyzz[k];

                g_zz_0_yyyy_zzz[k] = -g_zz_0_yyy_zzz[k] * ab_y + g_zz_0_yyy_yzzz[k];
            }

            /// Set up 860-870 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 860 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 861 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 862 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 863 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 864 * ccomps * dcomps);

            auto g_zz_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 865 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 866 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 867 * ccomps * dcomps);

            auto g_zz_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 868 * ccomps * dcomps);

            auto g_zz_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyz_xxx, g_zz_0_yyyz_xxy, g_zz_0_yyyz_xxz, g_zz_0_yyyz_xyy, g_zz_0_yyyz_xyz, g_zz_0_yyyz_xzz, g_zz_0_yyyz_yyy, g_zz_0_yyyz_yyz, g_zz_0_yyyz_yzz, g_zz_0_yyyz_zzz, g_zz_0_yyz_xxx, g_zz_0_yyz_xxxy, g_zz_0_yyz_xxy, g_zz_0_yyz_xxyy, g_zz_0_yyz_xxyz, g_zz_0_yyz_xxz, g_zz_0_yyz_xyy, g_zz_0_yyz_xyyy, g_zz_0_yyz_xyyz, g_zz_0_yyz_xyz, g_zz_0_yyz_xyzz, g_zz_0_yyz_xzz, g_zz_0_yyz_yyy, g_zz_0_yyz_yyyy, g_zz_0_yyz_yyyz, g_zz_0_yyz_yyz, g_zz_0_yyz_yyzz, g_zz_0_yyz_yzz, g_zz_0_yyz_yzzz, g_zz_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyz_xxx[k] = -g_zz_0_yyz_xxx[k] * ab_y + g_zz_0_yyz_xxxy[k];

                g_zz_0_yyyz_xxy[k] = -g_zz_0_yyz_xxy[k] * ab_y + g_zz_0_yyz_xxyy[k];

                g_zz_0_yyyz_xxz[k] = -g_zz_0_yyz_xxz[k] * ab_y + g_zz_0_yyz_xxyz[k];

                g_zz_0_yyyz_xyy[k] = -g_zz_0_yyz_xyy[k] * ab_y + g_zz_0_yyz_xyyy[k];

                g_zz_0_yyyz_xyz[k] = -g_zz_0_yyz_xyz[k] * ab_y + g_zz_0_yyz_xyyz[k];

                g_zz_0_yyyz_xzz[k] = -g_zz_0_yyz_xzz[k] * ab_y + g_zz_0_yyz_xyzz[k];

                g_zz_0_yyyz_yyy[k] = -g_zz_0_yyz_yyy[k] * ab_y + g_zz_0_yyz_yyyy[k];

                g_zz_0_yyyz_yyz[k] = -g_zz_0_yyz_yyz[k] * ab_y + g_zz_0_yyz_yyyz[k];

                g_zz_0_yyyz_yzz[k] = -g_zz_0_yyz_yzz[k] * ab_y + g_zz_0_yyz_yyzz[k];

                g_zz_0_yyyz_zzz[k] = -g_zz_0_yyz_zzz[k] * ab_y + g_zz_0_yyz_yzzz[k];
            }

            /// Set up 870-880 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 870 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 871 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 872 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 873 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 874 * ccomps * dcomps);

            auto g_zz_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 875 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 876 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 877 * ccomps * dcomps);

            auto g_zz_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 878 * ccomps * dcomps);

            auto g_zz_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 879 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyzz_xxx, g_zz_0_yyzz_xxy, g_zz_0_yyzz_xxz, g_zz_0_yyzz_xyy, g_zz_0_yyzz_xyz, g_zz_0_yyzz_xzz, g_zz_0_yyzz_yyy, g_zz_0_yyzz_yyz, g_zz_0_yyzz_yzz, g_zz_0_yyzz_zzz, g_zz_0_yzz_xxx, g_zz_0_yzz_xxxy, g_zz_0_yzz_xxy, g_zz_0_yzz_xxyy, g_zz_0_yzz_xxyz, g_zz_0_yzz_xxz, g_zz_0_yzz_xyy, g_zz_0_yzz_xyyy, g_zz_0_yzz_xyyz, g_zz_0_yzz_xyz, g_zz_0_yzz_xyzz, g_zz_0_yzz_xzz, g_zz_0_yzz_yyy, g_zz_0_yzz_yyyy, g_zz_0_yzz_yyyz, g_zz_0_yzz_yyz, g_zz_0_yzz_yyzz, g_zz_0_yzz_yzz, g_zz_0_yzz_yzzz, g_zz_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyzz_xxx[k] = -g_zz_0_yzz_xxx[k] * ab_y + g_zz_0_yzz_xxxy[k];

                g_zz_0_yyzz_xxy[k] = -g_zz_0_yzz_xxy[k] * ab_y + g_zz_0_yzz_xxyy[k];

                g_zz_0_yyzz_xxz[k] = -g_zz_0_yzz_xxz[k] * ab_y + g_zz_0_yzz_xxyz[k];

                g_zz_0_yyzz_xyy[k] = -g_zz_0_yzz_xyy[k] * ab_y + g_zz_0_yzz_xyyy[k];

                g_zz_0_yyzz_xyz[k] = -g_zz_0_yzz_xyz[k] * ab_y + g_zz_0_yzz_xyyz[k];

                g_zz_0_yyzz_xzz[k] = -g_zz_0_yzz_xzz[k] * ab_y + g_zz_0_yzz_xyzz[k];

                g_zz_0_yyzz_yyy[k] = -g_zz_0_yzz_yyy[k] * ab_y + g_zz_0_yzz_yyyy[k];

                g_zz_0_yyzz_yyz[k] = -g_zz_0_yzz_yyz[k] * ab_y + g_zz_0_yzz_yyyz[k];

                g_zz_0_yyzz_yzz[k] = -g_zz_0_yzz_yzz[k] * ab_y + g_zz_0_yzz_yyzz[k];

                g_zz_0_yyzz_zzz[k] = -g_zz_0_yzz_zzz[k] * ab_y + g_zz_0_yzz_yzzz[k];
            }

            /// Set up 880-890 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 880 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 881 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 882 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 883 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 884 * ccomps * dcomps);

            auto g_zz_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 885 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 886 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 887 * ccomps * dcomps);

            auto g_zz_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 888 * ccomps * dcomps);

            auto g_zz_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 889 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzzz_xxx, g_zz_0_yzzz_xxy, g_zz_0_yzzz_xxz, g_zz_0_yzzz_xyy, g_zz_0_yzzz_xyz, g_zz_0_yzzz_xzz, g_zz_0_yzzz_yyy, g_zz_0_yzzz_yyz, g_zz_0_yzzz_yzz, g_zz_0_yzzz_zzz, g_zz_0_zzz_xxx, g_zz_0_zzz_xxxy, g_zz_0_zzz_xxy, g_zz_0_zzz_xxyy, g_zz_0_zzz_xxyz, g_zz_0_zzz_xxz, g_zz_0_zzz_xyy, g_zz_0_zzz_xyyy, g_zz_0_zzz_xyyz, g_zz_0_zzz_xyz, g_zz_0_zzz_xyzz, g_zz_0_zzz_xzz, g_zz_0_zzz_yyy, g_zz_0_zzz_yyyy, g_zz_0_zzz_yyyz, g_zz_0_zzz_yyz, g_zz_0_zzz_yyzz, g_zz_0_zzz_yzz, g_zz_0_zzz_yzzz, g_zz_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzzz_xxx[k] = -g_zz_0_zzz_xxx[k] * ab_y + g_zz_0_zzz_xxxy[k];

                g_zz_0_yzzz_xxy[k] = -g_zz_0_zzz_xxy[k] * ab_y + g_zz_0_zzz_xxyy[k];

                g_zz_0_yzzz_xxz[k] = -g_zz_0_zzz_xxz[k] * ab_y + g_zz_0_zzz_xxyz[k];

                g_zz_0_yzzz_xyy[k] = -g_zz_0_zzz_xyy[k] * ab_y + g_zz_0_zzz_xyyy[k];

                g_zz_0_yzzz_xyz[k] = -g_zz_0_zzz_xyz[k] * ab_y + g_zz_0_zzz_xyyz[k];

                g_zz_0_yzzz_xzz[k] = -g_zz_0_zzz_xzz[k] * ab_y + g_zz_0_zzz_xyzz[k];

                g_zz_0_yzzz_yyy[k] = -g_zz_0_zzz_yyy[k] * ab_y + g_zz_0_zzz_yyyy[k];

                g_zz_0_yzzz_yyz[k] = -g_zz_0_zzz_yyz[k] * ab_y + g_zz_0_zzz_yyyz[k];

                g_zz_0_yzzz_yzz[k] = -g_zz_0_zzz_yzz[k] * ab_y + g_zz_0_zzz_yyzz[k];

                g_zz_0_yzzz_zzz[k] = -g_zz_0_zzz_zzz[k] * ab_y + g_zz_0_zzz_yzzz[k];
            }

            /// Set up 890-900 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 890 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 891 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 892 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 893 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 894 * ccomps * dcomps);

            auto g_zz_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 895 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 896 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 897 * ccomps * dcomps);

            auto g_zz_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 898 * ccomps * dcomps);

            auto g_zz_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_xxx, g_z_0_zzz_xxy, g_z_0_zzz_xxz, g_z_0_zzz_xyy, g_z_0_zzz_xyz, g_z_0_zzz_xzz, g_z_0_zzz_yyy, g_z_0_zzz_yyz, g_z_0_zzz_yzz, g_z_0_zzz_zzz, g_zz_0_zzz_xxx, g_zz_0_zzz_xxxz, g_zz_0_zzz_xxy, g_zz_0_zzz_xxyz, g_zz_0_zzz_xxz, g_zz_0_zzz_xxzz, g_zz_0_zzz_xyy, g_zz_0_zzz_xyyz, g_zz_0_zzz_xyz, g_zz_0_zzz_xyzz, g_zz_0_zzz_xzz, g_zz_0_zzz_xzzz, g_zz_0_zzz_yyy, g_zz_0_zzz_yyyz, g_zz_0_zzz_yyz, g_zz_0_zzz_yyzz, g_zz_0_zzz_yzz, g_zz_0_zzz_yzzz, g_zz_0_zzz_zzz, g_zz_0_zzz_zzzz, g_zz_0_zzzz_xxx, g_zz_0_zzzz_xxy, g_zz_0_zzzz_xxz, g_zz_0_zzzz_xyy, g_zz_0_zzzz_xyz, g_zz_0_zzzz_xzz, g_zz_0_zzzz_yyy, g_zz_0_zzzz_yyz, g_zz_0_zzzz_yzz, g_zz_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzzz_xxx[k] = -2.0 * g_z_0_zzz_xxx[k] - g_zz_0_zzz_xxx[k] * ab_z + g_zz_0_zzz_xxxz[k];

                g_zz_0_zzzz_xxy[k] = -2.0 * g_z_0_zzz_xxy[k] - g_zz_0_zzz_xxy[k] * ab_z + g_zz_0_zzz_xxyz[k];

                g_zz_0_zzzz_xxz[k] = -2.0 * g_z_0_zzz_xxz[k] - g_zz_0_zzz_xxz[k] * ab_z + g_zz_0_zzz_xxzz[k];

                g_zz_0_zzzz_xyy[k] = -2.0 * g_z_0_zzz_xyy[k] - g_zz_0_zzz_xyy[k] * ab_z + g_zz_0_zzz_xyyz[k];

                g_zz_0_zzzz_xyz[k] = -2.0 * g_z_0_zzz_xyz[k] - g_zz_0_zzz_xyz[k] * ab_z + g_zz_0_zzz_xyzz[k];

                g_zz_0_zzzz_xzz[k] = -2.0 * g_z_0_zzz_xzz[k] - g_zz_0_zzz_xzz[k] * ab_z + g_zz_0_zzz_xzzz[k];

                g_zz_0_zzzz_yyy[k] = -2.0 * g_z_0_zzz_yyy[k] - g_zz_0_zzz_yyy[k] * ab_z + g_zz_0_zzz_yyyz[k];

                g_zz_0_zzzz_yyz[k] = -2.0 * g_z_0_zzz_yyz[k] - g_zz_0_zzz_yyz[k] * ab_z + g_zz_0_zzz_yyzz[k];

                g_zz_0_zzzz_yzz[k] = -2.0 * g_z_0_zzz_yzz[k] - g_zz_0_zzz_yzz[k] * ab_z + g_zz_0_zzz_yzzz[k];

                g_zz_0_zzzz_zzz[k] = -2.0 * g_z_0_zzz_zzz[k] - g_zz_0_zzz_zzz[k] * ab_z + g_zz_0_zzz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

