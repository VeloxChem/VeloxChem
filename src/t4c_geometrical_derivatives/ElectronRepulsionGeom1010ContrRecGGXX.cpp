#include "ElectronRepulsionGeom1010ContrRecGGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_ggxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_ggxx,
                                              const size_t idx_geom_0010_fgxx,
                                              const size_t idx_geom_1010_fgxx,
                                              const size_t idx_geom_1010_fhxx,
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
            /// Set up components of auxilary buffer : FGSS

            const auto fg_geom_0010_off = idx_geom_0010_fgxx + i * dcomps + j;

            auto g_0_0_x_0_xxx_xxxx = cbuffer.data(fg_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_xxxy = cbuffer.data(fg_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_xxxz = cbuffer.data(fg_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_xxyy = cbuffer.data(fg_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_xxyz = cbuffer.data(fg_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_xxzz = cbuffer.data(fg_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_xyyy = cbuffer.data(fg_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_xyyz = cbuffer.data(fg_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_xyzz = cbuffer.data(fg_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_xzzz = cbuffer.data(fg_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_yyyy = cbuffer.data(fg_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_yyyz = cbuffer.data(fg_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_yyzz = cbuffer.data(fg_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_yzzz = cbuffer.data(fg_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_zzzz = cbuffer.data(fg_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xxxx = cbuffer.data(fg_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xxxy = cbuffer.data(fg_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xxxz = cbuffer.data(fg_geom_0010_off + 17 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xxyy = cbuffer.data(fg_geom_0010_off + 18 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xxyz = cbuffer.data(fg_geom_0010_off + 19 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xxzz = cbuffer.data(fg_geom_0010_off + 20 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xyyy = cbuffer.data(fg_geom_0010_off + 21 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xyyz = cbuffer.data(fg_geom_0010_off + 22 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xyzz = cbuffer.data(fg_geom_0010_off + 23 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_xzzz = cbuffer.data(fg_geom_0010_off + 24 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_yyyy = cbuffer.data(fg_geom_0010_off + 25 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_yyyz = cbuffer.data(fg_geom_0010_off + 26 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_yyzz = cbuffer.data(fg_geom_0010_off + 27 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_yzzz = cbuffer.data(fg_geom_0010_off + 28 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_zzzz = cbuffer.data(fg_geom_0010_off + 29 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xxxx = cbuffer.data(fg_geom_0010_off + 30 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xxxy = cbuffer.data(fg_geom_0010_off + 31 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xxxz = cbuffer.data(fg_geom_0010_off + 32 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xxyy = cbuffer.data(fg_geom_0010_off + 33 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xxyz = cbuffer.data(fg_geom_0010_off + 34 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xxzz = cbuffer.data(fg_geom_0010_off + 35 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xyyy = cbuffer.data(fg_geom_0010_off + 36 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xyyz = cbuffer.data(fg_geom_0010_off + 37 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xyzz = cbuffer.data(fg_geom_0010_off + 38 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_xzzz = cbuffer.data(fg_geom_0010_off + 39 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_yyyy = cbuffer.data(fg_geom_0010_off + 40 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_yyyz = cbuffer.data(fg_geom_0010_off + 41 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_yyzz = cbuffer.data(fg_geom_0010_off + 42 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_yzzz = cbuffer.data(fg_geom_0010_off + 43 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_zzzz = cbuffer.data(fg_geom_0010_off + 44 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xxxx = cbuffer.data(fg_geom_0010_off + 45 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xxxy = cbuffer.data(fg_geom_0010_off + 46 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xxxz = cbuffer.data(fg_geom_0010_off + 47 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xxyy = cbuffer.data(fg_geom_0010_off + 48 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xxyz = cbuffer.data(fg_geom_0010_off + 49 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xxzz = cbuffer.data(fg_geom_0010_off + 50 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xyyy = cbuffer.data(fg_geom_0010_off + 51 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xyyz = cbuffer.data(fg_geom_0010_off + 52 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xyzz = cbuffer.data(fg_geom_0010_off + 53 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_xzzz = cbuffer.data(fg_geom_0010_off + 54 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_yyyy = cbuffer.data(fg_geom_0010_off + 55 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_yyyz = cbuffer.data(fg_geom_0010_off + 56 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_yyzz = cbuffer.data(fg_geom_0010_off + 57 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_yzzz = cbuffer.data(fg_geom_0010_off + 58 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_zzzz = cbuffer.data(fg_geom_0010_off + 59 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xxxx = cbuffer.data(fg_geom_0010_off + 60 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xxxy = cbuffer.data(fg_geom_0010_off + 61 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xxxz = cbuffer.data(fg_geom_0010_off + 62 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xxyy = cbuffer.data(fg_geom_0010_off + 63 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xxyz = cbuffer.data(fg_geom_0010_off + 64 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xxzz = cbuffer.data(fg_geom_0010_off + 65 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xyyy = cbuffer.data(fg_geom_0010_off + 66 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xyyz = cbuffer.data(fg_geom_0010_off + 67 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xyzz = cbuffer.data(fg_geom_0010_off + 68 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_xzzz = cbuffer.data(fg_geom_0010_off + 69 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_yyyy = cbuffer.data(fg_geom_0010_off + 70 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_yyyz = cbuffer.data(fg_geom_0010_off + 71 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_yyzz = cbuffer.data(fg_geom_0010_off + 72 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_yzzz = cbuffer.data(fg_geom_0010_off + 73 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_zzzz = cbuffer.data(fg_geom_0010_off + 74 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xxxx = cbuffer.data(fg_geom_0010_off + 75 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xxxy = cbuffer.data(fg_geom_0010_off + 76 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xxxz = cbuffer.data(fg_geom_0010_off + 77 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xxyy = cbuffer.data(fg_geom_0010_off + 78 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xxyz = cbuffer.data(fg_geom_0010_off + 79 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xxzz = cbuffer.data(fg_geom_0010_off + 80 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xyyy = cbuffer.data(fg_geom_0010_off + 81 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xyyz = cbuffer.data(fg_geom_0010_off + 82 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xyzz = cbuffer.data(fg_geom_0010_off + 83 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_xzzz = cbuffer.data(fg_geom_0010_off + 84 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_yyyy = cbuffer.data(fg_geom_0010_off + 85 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_yyyz = cbuffer.data(fg_geom_0010_off + 86 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_yyzz = cbuffer.data(fg_geom_0010_off + 87 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_yzzz = cbuffer.data(fg_geom_0010_off + 88 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_zzzz = cbuffer.data(fg_geom_0010_off + 89 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xxxx = cbuffer.data(fg_geom_0010_off + 90 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xxxy = cbuffer.data(fg_geom_0010_off + 91 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xxxz = cbuffer.data(fg_geom_0010_off + 92 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xxyy = cbuffer.data(fg_geom_0010_off + 93 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xxyz = cbuffer.data(fg_geom_0010_off + 94 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xxzz = cbuffer.data(fg_geom_0010_off + 95 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xyyy = cbuffer.data(fg_geom_0010_off + 96 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xyyz = cbuffer.data(fg_geom_0010_off + 97 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xyzz = cbuffer.data(fg_geom_0010_off + 98 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_xzzz = cbuffer.data(fg_geom_0010_off + 99 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_yyyy = cbuffer.data(fg_geom_0010_off + 100 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_yyyz = cbuffer.data(fg_geom_0010_off + 101 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_yyzz = cbuffer.data(fg_geom_0010_off + 102 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_yzzz = cbuffer.data(fg_geom_0010_off + 103 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_zzzz = cbuffer.data(fg_geom_0010_off + 104 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xxxx = cbuffer.data(fg_geom_0010_off + 105 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xxxy = cbuffer.data(fg_geom_0010_off + 106 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xxxz = cbuffer.data(fg_geom_0010_off + 107 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xxyy = cbuffer.data(fg_geom_0010_off + 108 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xxyz = cbuffer.data(fg_geom_0010_off + 109 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xxzz = cbuffer.data(fg_geom_0010_off + 110 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xyyy = cbuffer.data(fg_geom_0010_off + 111 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xyyz = cbuffer.data(fg_geom_0010_off + 112 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xyzz = cbuffer.data(fg_geom_0010_off + 113 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_xzzz = cbuffer.data(fg_geom_0010_off + 114 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_yyyy = cbuffer.data(fg_geom_0010_off + 115 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_yyyz = cbuffer.data(fg_geom_0010_off + 116 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_yyzz = cbuffer.data(fg_geom_0010_off + 117 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_yzzz = cbuffer.data(fg_geom_0010_off + 118 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_zzzz = cbuffer.data(fg_geom_0010_off + 119 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xxxx = cbuffer.data(fg_geom_0010_off + 120 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xxxy = cbuffer.data(fg_geom_0010_off + 121 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xxxz = cbuffer.data(fg_geom_0010_off + 122 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xxyy = cbuffer.data(fg_geom_0010_off + 123 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xxyz = cbuffer.data(fg_geom_0010_off + 124 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xxzz = cbuffer.data(fg_geom_0010_off + 125 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xyyy = cbuffer.data(fg_geom_0010_off + 126 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xyyz = cbuffer.data(fg_geom_0010_off + 127 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xyzz = cbuffer.data(fg_geom_0010_off + 128 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_xzzz = cbuffer.data(fg_geom_0010_off + 129 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_yyyy = cbuffer.data(fg_geom_0010_off + 130 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_yyyz = cbuffer.data(fg_geom_0010_off + 131 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_yyzz = cbuffer.data(fg_geom_0010_off + 132 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_yzzz = cbuffer.data(fg_geom_0010_off + 133 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_zzzz = cbuffer.data(fg_geom_0010_off + 134 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xxxx = cbuffer.data(fg_geom_0010_off + 135 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xxxy = cbuffer.data(fg_geom_0010_off + 136 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xxxz = cbuffer.data(fg_geom_0010_off + 137 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xxyy = cbuffer.data(fg_geom_0010_off + 138 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xxyz = cbuffer.data(fg_geom_0010_off + 139 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xxzz = cbuffer.data(fg_geom_0010_off + 140 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xyyy = cbuffer.data(fg_geom_0010_off + 141 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xyyz = cbuffer.data(fg_geom_0010_off + 142 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xyzz = cbuffer.data(fg_geom_0010_off + 143 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_xzzz = cbuffer.data(fg_geom_0010_off + 144 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_yyyy = cbuffer.data(fg_geom_0010_off + 145 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_yyyz = cbuffer.data(fg_geom_0010_off + 146 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_yyzz = cbuffer.data(fg_geom_0010_off + 147 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_yzzz = cbuffer.data(fg_geom_0010_off + 148 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_zzzz = cbuffer.data(fg_geom_0010_off + 149 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xxxx = cbuffer.data(fg_geom_0010_off + 150 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xxxy = cbuffer.data(fg_geom_0010_off + 151 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xxxz = cbuffer.data(fg_geom_0010_off + 152 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xxyy = cbuffer.data(fg_geom_0010_off + 153 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xxyz = cbuffer.data(fg_geom_0010_off + 154 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xxzz = cbuffer.data(fg_geom_0010_off + 155 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xyyy = cbuffer.data(fg_geom_0010_off + 156 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xyyz = cbuffer.data(fg_geom_0010_off + 157 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xyzz = cbuffer.data(fg_geom_0010_off + 158 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_xzzz = cbuffer.data(fg_geom_0010_off + 159 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_yyyy = cbuffer.data(fg_geom_0010_off + 160 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_yyyz = cbuffer.data(fg_geom_0010_off + 161 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_yyzz = cbuffer.data(fg_geom_0010_off + 162 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_yzzz = cbuffer.data(fg_geom_0010_off + 163 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_zzzz = cbuffer.data(fg_geom_0010_off + 164 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xxxx = cbuffer.data(fg_geom_0010_off + 165 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xxxy = cbuffer.data(fg_geom_0010_off + 166 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xxxz = cbuffer.data(fg_geom_0010_off + 167 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xxyy = cbuffer.data(fg_geom_0010_off + 168 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xxyz = cbuffer.data(fg_geom_0010_off + 169 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xxzz = cbuffer.data(fg_geom_0010_off + 170 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xyyy = cbuffer.data(fg_geom_0010_off + 171 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xyyz = cbuffer.data(fg_geom_0010_off + 172 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xyzz = cbuffer.data(fg_geom_0010_off + 173 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_xzzz = cbuffer.data(fg_geom_0010_off + 174 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_yyyy = cbuffer.data(fg_geom_0010_off + 175 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_yyyz = cbuffer.data(fg_geom_0010_off + 176 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_yyzz = cbuffer.data(fg_geom_0010_off + 177 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_yzzz = cbuffer.data(fg_geom_0010_off + 178 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_zzzz = cbuffer.data(fg_geom_0010_off + 179 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xxxx = cbuffer.data(fg_geom_0010_off + 180 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xxxy = cbuffer.data(fg_geom_0010_off + 181 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xxxz = cbuffer.data(fg_geom_0010_off + 182 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xxyy = cbuffer.data(fg_geom_0010_off + 183 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xxyz = cbuffer.data(fg_geom_0010_off + 184 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xxzz = cbuffer.data(fg_geom_0010_off + 185 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xyyy = cbuffer.data(fg_geom_0010_off + 186 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xyyz = cbuffer.data(fg_geom_0010_off + 187 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xyzz = cbuffer.data(fg_geom_0010_off + 188 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_xzzz = cbuffer.data(fg_geom_0010_off + 189 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_yyyy = cbuffer.data(fg_geom_0010_off + 190 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_yyyz = cbuffer.data(fg_geom_0010_off + 191 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_yyzz = cbuffer.data(fg_geom_0010_off + 192 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_yzzz = cbuffer.data(fg_geom_0010_off + 193 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_zzzz = cbuffer.data(fg_geom_0010_off + 194 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xxxx = cbuffer.data(fg_geom_0010_off + 195 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xxxy = cbuffer.data(fg_geom_0010_off + 196 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xxxz = cbuffer.data(fg_geom_0010_off + 197 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xxyy = cbuffer.data(fg_geom_0010_off + 198 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xxyz = cbuffer.data(fg_geom_0010_off + 199 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xxzz = cbuffer.data(fg_geom_0010_off + 200 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xyyy = cbuffer.data(fg_geom_0010_off + 201 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xyyz = cbuffer.data(fg_geom_0010_off + 202 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xyzz = cbuffer.data(fg_geom_0010_off + 203 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_xzzz = cbuffer.data(fg_geom_0010_off + 204 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_yyyy = cbuffer.data(fg_geom_0010_off + 205 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_yyyz = cbuffer.data(fg_geom_0010_off + 206 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_yyzz = cbuffer.data(fg_geom_0010_off + 207 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_yzzz = cbuffer.data(fg_geom_0010_off + 208 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_zzzz = cbuffer.data(fg_geom_0010_off + 209 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xxxx = cbuffer.data(fg_geom_0010_off + 210 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xxxy = cbuffer.data(fg_geom_0010_off + 211 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xxxz = cbuffer.data(fg_geom_0010_off + 212 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xxyy = cbuffer.data(fg_geom_0010_off + 213 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xxyz = cbuffer.data(fg_geom_0010_off + 214 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xxzz = cbuffer.data(fg_geom_0010_off + 215 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xyyy = cbuffer.data(fg_geom_0010_off + 216 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xyyz = cbuffer.data(fg_geom_0010_off + 217 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xyzz = cbuffer.data(fg_geom_0010_off + 218 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_xzzz = cbuffer.data(fg_geom_0010_off + 219 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_yyyy = cbuffer.data(fg_geom_0010_off + 220 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_yyyz = cbuffer.data(fg_geom_0010_off + 221 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_yyzz = cbuffer.data(fg_geom_0010_off + 222 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_yzzz = cbuffer.data(fg_geom_0010_off + 223 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_zzzz = cbuffer.data(fg_geom_0010_off + 224 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xxxx = cbuffer.data(fg_geom_0010_off + 225 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xxxy = cbuffer.data(fg_geom_0010_off + 226 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xxxz = cbuffer.data(fg_geom_0010_off + 227 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xxyy = cbuffer.data(fg_geom_0010_off + 228 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xxyz = cbuffer.data(fg_geom_0010_off + 229 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xxzz = cbuffer.data(fg_geom_0010_off + 230 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xyyy = cbuffer.data(fg_geom_0010_off + 231 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xyyz = cbuffer.data(fg_geom_0010_off + 232 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xyzz = cbuffer.data(fg_geom_0010_off + 233 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_xzzz = cbuffer.data(fg_geom_0010_off + 234 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_yyyy = cbuffer.data(fg_geom_0010_off + 235 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_yyyz = cbuffer.data(fg_geom_0010_off + 236 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_yyzz = cbuffer.data(fg_geom_0010_off + 237 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_yzzz = cbuffer.data(fg_geom_0010_off + 238 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_zzzz = cbuffer.data(fg_geom_0010_off + 239 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xxxx = cbuffer.data(fg_geom_0010_off + 240 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xxxy = cbuffer.data(fg_geom_0010_off + 241 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xxxz = cbuffer.data(fg_geom_0010_off + 242 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xxyy = cbuffer.data(fg_geom_0010_off + 243 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xxyz = cbuffer.data(fg_geom_0010_off + 244 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xxzz = cbuffer.data(fg_geom_0010_off + 245 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xyyy = cbuffer.data(fg_geom_0010_off + 246 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xyyz = cbuffer.data(fg_geom_0010_off + 247 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xyzz = cbuffer.data(fg_geom_0010_off + 248 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_xzzz = cbuffer.data(fg_geom_0010_off + 249 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_yyyy = cbuffer.data(fg_geom_0010_off + 250 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_yyyz = cbuffer.data(fg_geom_0010_off + 251 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_yyzz = cbuffer.data(fg_geom_0010_off + 252 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_yzzz = cbuffer.data(fg_geom_0010_off + 253 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_zzzz = cbuffer.data(fg_geom_0010_off + 254 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xxxx = cbuffer.data(fg_geom_0010_off + 255 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xxxy = cbuffer.data(fg_geom_0010_off + 256 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xxxz = cbuffer.data(fg_geom_0010_off + 257 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xxyy = cbuffer.data(fg_geom_0010_off + 258 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xxyz = cbuffer.data(fg_geom_0010_off + 259 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xxzz = cbuffer.data(fg_geom_0010_off + 260 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xyyy = cbuffer.data(fg_geom_0010_off + 261 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xyyz = cbuffer.data(fg_geom_0010_off + 262 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xyzz = cbuffer.data(fg_geom_0010_off + 263 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_xzzz = cbuffer.data(fg_geom_0010_off + 264 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_yyyy = cbuffer.data(fg_geom_0010_off + 265 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_yyyz = cbuffer.data(fg_geom_0010_off + 266 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_yyzz = cbuffer.data(fg_geom_0010_off + 267 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_yzzz = cbuffer.data(fg_geom_0010_off + 268 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_zzzz = cbuffer.data(fg_geom_0010_off + 269 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xxxx = cbuffer.data(fg_geom_0010_off + 270 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xxxy = cbuffer.data(fg_geom_0010_off + 271 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xxxz = cbuffer.data(fg_geom_0010_off + 272 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xxyy = cbuffer.data(fg_geom_0010_off + 273 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xxyz = cbuffer.data(fg_geom_0010_off + 274 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xxzz = cbuffer.data(fg_geom_0010_off + 275 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xyyy = cbuffer.data(fg_geom_0010_off + 276 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xyyz = cbuffer.data(fg_geom_0010_off + 277 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xyzz = cbuffer.data(fg_geom_0010_off + 278 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_xzzz = cbuffer.data(fg_geom_0010_off + 279 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_yyyy = cbuffer.data(fg_geom_0010_off + 280 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_yyyz = cbuffer.data(fg_geom_0010_off + 281 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_yyzz = cbuffer.data(fg_geom_0010_off + 282 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_yzzz = cbuffer.data(fg_geom_0010_off + 283 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_zzzz = cbuffer.data(fg_geom_0010_off + 284 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xxxx = cbuffer.data(fg_geom_0010_off + 285 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xxxy = cbuffer.data(fg_geom_0010_off + 286 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xxxz = cbuffer.data(fg_geom_0010_off + 287 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xxyy = cbuffer.data(fg_geom_0010_off + 288 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xxyz = cbuffer.data(fg_geom_0010_off + 289 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xxzz = cbuffer.data(fg_geom_0010_off + 290 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xyyy = cbuffer.data(fg_geom_0010_off + 291 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xyyz = cbuffer.data(fg_geom_0010_off + 292 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xyzz = cbuffer.data(fg_geom_0010_off + 293 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_xzzz = cbuffer.data(fg_geom_0010_off + 294 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_yyyy = cbuffer.data(fg_geom_0010_off + 295 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_yyyz = cbuffer.data(fg_geom_0010_off + 296 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_yyzz = cbuffer.data(fg_geom_0010_off + 297 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_yzzz = cbuffer.data(fg_geom_0010_off + 298 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_zzzz = cbuffer.data(fg_geom_0010_off + 299 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xxxx = cbuffer.data(fg_geom_0010_off + 300 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xxxy = cbuffer.data(fg_geom_0010_off + 301 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xxxz = cbuffer.data(fg_geom_0010_off + 302 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xxyy = cbuffer.data(fg_geom_0010_off + 303 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xxyz = cbuffer.data(fg_geom_0010_off + 304 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xxzz = cbuffer.data(fg_geom_0010_off + 305 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xyyy = cbuffer.data(fg_geom_0010_off + 306 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xyyz = cbuffer.data(fg_geom_0010_off + 307 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xyzz = cbuffer.data(fg_geom_0010_off + 308 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_xzzz = cbuffer.data(fg_geom_0010_off + 309 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_yyyy = cbuffer.data(fg_geom_0010_off + 310 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_yyyz = cbuffer.data(fg_geom_0010_off + 311 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_yyzz = cbuffer.data(fg_geom_0010_off + 312 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_yzzz = cbuffer.data(fg_geom_0010_off + 313 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_zzzz = cbuffer.data(fg_geom_0010_off + 314 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xxxx = cbuffer.data(fg_geom_0010_off + 315 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xxxy = cbuffer.data(fg_geom_0010_off + 316 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xxxz = cbuffer.data(fg_geom_0010_off + 317 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xxyy = cbuffer.data(fg_geom_0010_off + 318 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xxyz = cbuffer.data(fg_geom_0010_off + 319 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xxzz = cbuffer.data(fg_geom_0010_off + 320 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xyyy = cbuffer.data(fg_geom_0010_off + 321 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xyyz = cbuffer.data(fg_geom_0010_off + 322 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xyzz = cbuffer.data(fg_geom_0010_off + 323 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_xzzz = cbuffer.data(fg_geom_0010_off + 324 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_yyyy = cbuffer.data(fg_geom_0010_off + 325 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_yyyz = cbuffer.data(fg_geom_0010_off + 326 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_yyzz = cbuffer.data(fg_geom_0010_off + 327 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_yzzz = cbuffer.data(fg_geom_0010_off + 328 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_zzzz = cbuffer.data(fg_geom_0010_off + 329 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xxxx = cbuffer.data(fg_geom_0010_off + 330 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xxxy = cbuffer.data(fg_geom_0010_off + 331 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xxxz = cbuffer.data(fg_geom_0010_off + 332 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xxyy = cbuffer.data(fg_geom_0010_off + 333 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xxyz = cbuffer.data(fg_geom_0010_off + 334 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xxzz = cbuffer.data(fg_geom_0010_off + 335 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xyyy = cbuffer.data(fg_geom_0010_off + 336 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xyyz = cbuffer.data(fg_geom_0010_off + 337 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xyzz = cbuffer.data(fg_geom_0010_off + 338 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_xzzz = cbuffer.data(fg_geom_0010_off + 339 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_yyyy = cbuffer.data(fg_geom_0010_off + 340 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_yyyz = cbuffer.data(fg_geom_0010_off + 341 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_yyzz = cbuffer.data(fg_geom_0010_off + 342 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_yzzz = cbuffer.data(fg_geom_0010_off + 343 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_zzzz = cbuffer.data(fg_geom_0010_off + 344 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xxxx = cbuffer.data(fg_geom_0010_off + 345 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xxxy = cbuffer.data(fg_geom_0010_off + 346 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xxxz = cbuffer.data(fg_geom_0010_off + 347 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xxyy = cbuffer.data(fg_geom_0010_off + 348 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xxyz = cbuffer.data(fg_geom_0010_off + 349 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xxzz = cbuffer.data(fg_geom_0010_off + 350 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xyyy = cbuffer.data(fg_geom_0010_off + 351 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xyyz = cbuffer.data(fg_geom_0010_off + 352 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xyzz = cbuffer.data(fg_geom_0010_off + 353 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_xzzz = cbuffer.data(fg_geom_0010_off + 354 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_yyyy = cbuffer.data(fg_geom_0010_off + 355 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_yyyz = cbuffer.data(fg_geom_0010_off + 356 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_yyzz = cbuffer.data(fg_geom_0010_off + 357 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_yzzz = cbuffer.data(fg_geom_0010_off + 358 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_zzzz = cbuffer.data(fg_geom_0010_off + 359 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xxxx = cbuffer.data(fg_geom_0010_off + 360 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xxxy = cbuffer.data(fg_geom_0010_off + 361 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xxxz = cbuffer.data(fg_geom_0010_off + 362 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xxyy = cbuffer.data(fg_geom_0010_off + 363 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xxyz = cbuffer.data(fg_geom_0010_off + 364 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xxzz = cbuffer.data(fg_geom_0010_off + 365 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xyyy = cbuffer.data(fg_geom_0010_off + 366 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xyyz = cbuffer.data(fg_geom_0010_off + 367 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xyzz = cbuffer.data(fg_geom_0010_off + 368 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_xzzz = cbuffer.data(fg_geom_0010_off + 369 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_yyyy = cbuffer.data(fg_geom_0010_off + 370 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_yyyz = cbuffer.data(fg_geom_0010_off + 371 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_yyzz = cbuffer.data(fg_geom_0010_off + 372 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_yzzz = cbuffer.data(fg_geom_0010_off + 373 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_zzzz = cbuffer.data(fg_geom_0010_off + 374 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xxxx = cbuffer.data(fg_geom_0010_off + 375 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xxxy = cbuffer.data(fg_geom_0010_off + 376 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xxxz = cbuffer.data(fg_geom_0010_off + 377 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xxyy = cbuffer.data(fg_geom_0010_off + 378 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xxyz = cbuffer.data(fg_geom_0010_off + 379 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xxzz = cbuffer.data(fg_geom_0010_off + 380 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xyyy = cbuffer.data(fg_geom_0010_off + 381 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xyyz = cbuffer.data(fg_geom_0010_off + 382 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xyzz = cbuffer.data(fg_geom_0010_off + 383 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_xzzz = cbuffer.data(fg_geom_0010_off + 384 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_yyyy = cbuffer.data(fg_geom_0010_off + 385 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_yyyz = cbuffer.data(fg_geom_0010_off + 386 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_yyzz = cbuffer.data(fg_geom_0010_off + 387 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_yzzz = cbuffer.data(fg_geom_0010_off + 388 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_zzzz = cbuffer.data(fg_geom_0010_off + 389 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xxxx = cbuffer.data(fg_geom_0010_off + 390 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xxxy = cbuffer.data(fg_geom_0010_off + 391 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xxxz = cbuffer.data(fg_geom_0010_off + 392 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xxyy = cbuffer.data(fg_geom_0010_off + 393 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xxyz = cbuffer.data(fg_geom_0010_off + 394 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xxzz = cbuffer.data(fg_geom_0010_off + 395 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xyyy = cbuffer.data(fg_geom_0010_off + 396 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xyyz = cbuffer.data(fg_geom_0010_off + 397 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xyzz = cbuffer.data(fg_geom_0010_off + 398 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_xzzz = cbuffer.data(fg_geom_0010_off + 399 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_yyyy = cbuffer.data(fg_geom_0010_off + 400 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_yyyz = cbuffer.data(fg_geom_0010_off + 401 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_yyzz = cbuffer.data(fg_geom_0010_off + 402 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_yzzz = cbuffer.data(fg_geom_0010_off + 403 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_zzzz = cbuffer.data(fg_geom_0010_off + 404 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xxxx = cbuffer.data(fg_geom_0010_off + 405 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xxxy = cbuffer.data(fg_geom_0010_off + 406 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xxxz = cbuffer.data(fg_geom_0010_off + 407 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xxyy = cbuffer.data(fg_geom_0010_off + 408 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xxyz = cbuffer.data(fg_geom_0010_off + 409 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xxzz = cbuffer.data(fg_geom_0010_off + 410 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xyyy = cbuffer.data(fg_geom_0010_off + 411 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xyyz = cbuffer.data(fg_geom_0010_off + 412 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xyzz = cbuffer.data(fg_geom_0010_off + 413 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_xzzz = cbuffer.data(fg_geom_0010_off + 414 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_yyyy = cbuffer.data(fg_geom_0010_off + 415 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_yyyz = cbuffer.data(fg_geom_0010_off + 416 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_yyzz = cbuffer.data(fg_geom_0010_off + 417 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_yzzz = cbuffer.data(fg_geom_0010_off + 418 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_zzzz = cbuffer.data(fg_geom_0010_off + 419 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xxxx = cbuffer.data(fg_geom_0010_off + 420 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xxxy = cbuffer.data(fg_geom_0010_off + 421 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xxxz = cbuffer.data(fg_geom_0010_off + 422 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xxyy = cbuffer.data(fg_geom_0010_off + 423 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xxyz = cbuffer.data(fg_geom_0010_off + 424 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xxzz = cbuffer.data(fg_geom_0010_off + 425 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xyyy = cbuffer.data(fg_geom_0010_off + 426 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xyyz = cbuffer.data(fg_geom_0010_off + 427 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xyzz = cbuffer.data(fg_geom_0010_off + 428 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_xzzz = cbuffer.data(fg_geom_0010_off + 429 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_yyyy = cbuffer.data(fg_geom_0010_off + 430 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_yyyz = cbuffer.data(fg_geom_0010_off + 431 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_yyzz = cbuffer.data(fg_geom_0010_off + 432 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_yzzz = cbuffer.data(fg_geom_0010_off + 433 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_zzzz = cbuffer.data(fg_geom_0010_off + 434 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xxxx = cbuffer.data(fg_geom_0010_off + 435 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xxxy = cbuffer.data(fg_geom_0010_off + 436 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xxxz = cbuffer.data(fg_geom_0010_off + 437 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xxyy = cbuffer.data(fg_geom_0010_off + 438 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xxyz = cbuffer.data(fg_geom_0010_off + 439 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xxzz = cbuffer.data(fg_geom_0010_off + 440 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xyyy = cbuffer.data(fg_geom_0010_off + 441 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xyyz = cbuffer.data(fg_geom_0010_off + 442 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xyzz = cbuffer.data(fg_geom_0010_off + 443 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_xzzz = cbuffer.data(fg_geom_0010_off + 444 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_yyyy = cbuffer.data(fg_geom_0010_off + 445 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_yyyz = cbuffer.data(fg_geom_0010_off + 446 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_yyzz = cbuffer.data(fg_geom_0010_off + 447 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_yzzz = cbuffer.data(fg_geom_0010_off + 448 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_zzzz = cbuffer.data(fg_geom_0010_off + 449 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FGSS

            const auto fg_geom_1010_off = idx_geom_1010_fgxx + i * dcomps + j;

            auto g_x_0_x_0_xxx_xxxx = cbuffer.data(fg_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxxy = cbuffer.data(fg_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxxz = cbuffer.data(fg_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxyy = cbuffer.data(fg_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxyz = cbuffer.data(fg_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxzz = cbuffer.data(fg_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xyyy = cbuffer.data(fg_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xyyz = cbuffer.data(fg_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xyzz = cbuffer.data(fg_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xzzz = cbuffer.data(fg_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yyyy = cbuffer.data(fg_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yyyz = cbuffer.data(fg_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yyzz = cbuffer.data(fg_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yzzz = cbuffer.data(fg_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_zzzz = cbuffer.data(fg_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxxx = cbuffer.data(fg_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxxy = cbuffer.data(fg_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxxz = cbuffer.data(fg_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxyy = cbuffer.data(fg_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxyz = cbuffer.data(fg_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxzz = cbuffer.data(fg_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xyyy = cbuffer.data(fg_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xyyz = cbuffer.data(fg_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xyzz = cbuffer.data(fg_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xzzz = cbuffer.data(fg_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yyyy = cbuffer.data(fg_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yyyz = cbuffer.data(fg_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yyzz = cbuffer.data(fg_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yzzz = cbuffer.data(fg_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_zzzz = cbuffer.data(fg_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxxx = cbuffer.data(fg_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxxy = cbuffer.data(fg_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxxz = cbuffer.data(fg_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxyy = cbuffer.data(fg_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxyz = cbuffer.data(fg_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxzz = cbuffer.data(fg_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xyyy = cbuffer.data(fg_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xyyz = cbuffer.data(fg_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xyzz = cbuffer.data(fg_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xzzz = cbuffer.data(fg_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yyyy = cbuffer.data(fg_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yyyz = cbuffer.data(fg_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yyzz = cbuffer.data(fg_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yzzz = cbuffer.data(fg_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_zzzz = cbuffer.data(fg_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxxx = cbuffer.data(fg_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxxy = cbuffer.data(fg_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxxz = cbuffer.data(fg_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxyy = cbuffer.data(fg_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxyz = cbuffer.data(fg_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxzz = cbuffer.data(fg_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xyyy = cbuffer.data(fg_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xyyz = cbuffer.data(fg_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xyzz = cbuffer.data(fg_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xzzz = cbuffer.data(fg_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yyyy = cbuffer.data(fg_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yyyz = cbuffer.data(fg_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yyzz = cbuffer.data(fg_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yzzz = cbuffer.data(fg_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_zzzz = cbuffer.data(fg_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxxx = cbuffer.data(fg_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxxy = cbuffer.data(fg_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxxz = cbuffer.data(fg_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxyy = cbuffer.data(fg_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxyz = cbuffer.data(fg_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxzz = cbuffer.data(fg_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xyyy = cbuffer.data(fg_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xyyz = cbuffer.data(fg_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xyzz = cbuffer.data(fg_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xzzz = cbuffer.data(fg_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yyyy = cbuffer.data(fg_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yyyz = cbuffer.data(fg_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yyzz = cbuffer.data(fg_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yzzz = cbuffer.data(fg_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_zzzz = cbuffer.data(fg_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxxx = cbuffer.data(fg_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxxy = cbuffer.data(fg_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxxz = cbuffer.data(fg_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxyy = cbuffer.data(fg_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxyz = cbuffer.data(fg_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxzz = cbuffer.data(fg_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xyyy = cbuffer.data(fg_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xyyz = cbuffer.data(fg_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xyzz = cbuffer.data(fg_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xzzz = cbuffer.data(fg_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yyyy = cbuffer.data(fg_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yyyz = cbuffer.data(fg_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yyzz = cbuffer.data(fg_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yzzz = cbuffer.data(fg_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_zzzz = cbuffer.data(fg_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxxx = cbuffer.data(fg_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxxy = cbuffer.data(fg_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxxz = cbuffer.data(fg_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxyy = cbuffer.data(fg_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxyz = cbuffer.data(fg_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxzz = cbuffer.data(fg_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xyyy = cbuffer.data(fg_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xyyz = cbuffer.data(fg_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xyzz = cbuffer.data(fg_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xzzz = cbuffer.data(fg_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yyyy = cbuffer.data(fg_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yyyz = cbuffer.data(fg_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yyzz = cbuffer.data(fg_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yzzz = cbuffer.data(fg_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_zzzz = cbuffer.data(fg_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxxx = cbuffer.data(fg_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxxy = cbuffer.data(fg_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxxz = cbuffer.data(fg_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxyy = cbuffer.data(fg_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxyz = cbuffer.data(fg_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxzz = cbuffer.data(fg_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xyyy = cbuffer.data(fg_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xyyz = cbuffer.data(fg_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xyzz = cbuffer.data(fg_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xzzz = cbuffer.data(fg_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yyyy = cbuffer.data(fg_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yyyz = cbuffer.data(fg_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yyzz = cbuffer.data(fg_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yzzz = cbuffer.data(fg_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_zzzz = cbuffer.data(fg_geom_1010_off + 119 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxxx = cbuffer.data(fg_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxxy = cbuffer.data(fg_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxxz = cbuffer.data(fg_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxyy = cbuffer.data(fg_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxyz = cbuffer.data(fg_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxzz = cbuffer.data(fg_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xyyy = cbuffer.data(fg_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xyyz = cbuffer.data(fg_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xyzz = cbuffer.data(fg_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xzzz = cbuffer.data(fg_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yyyy = cbuffer.data(fg_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yyyz = cbuffer.data(fg_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yyzz = cbuffer.data(fg_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yzzz = cbuffer.data(fg_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_zzzz = cbuffer.data(fg_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxxx = cbuffer.data(fg_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxxy = cbuffer.data(fg_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxxz = cbuffer.data(fg_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxyy = cbuffer.data(fg_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxyz = cbuffer.data(fg_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxzz = cbuffer.data(fg_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xyyy = cbuffer.data(fg_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xyyz = cbuffer.data(fg_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xyzz = cbuffer.data(fg_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xzzz = cbuffer.data(fg_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yyyy = cbuffer.data(fg_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yyyz = cbuffer.data(fg_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yyzz = cbuffer.data(fg_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yzzz = cbuffer.data(fg_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_zzzz = cbuffer.data(fg_geom_1010_off + 149 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxxx = cbuffer.data(fg_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxxy = cbuffer.data(fg_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxxz = cbuffer.data(fg_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxyy = cbuffer.data(fg_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxyz = cbuffer.data(fg_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxzz = cbuffer.data(fg_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xyyy = cbuffer.data(fg_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xyyz = cbuffer.data(fg_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xyzz = cbuffer.data(fg_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xzzz = cbuffer.data(fg_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yyyy = cbuffer.data(fg_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yyyz = cbuffer.data(fg_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yyzz = cbuffer.data(fg_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yzzz = cbuffer.data(fg_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_zzzz = cbuffer.data(fg_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxxx = cbuffer.data(fg_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxxy = cbuffer.data(fg_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxxz = cbuffer.data(fg_geom_1010_off + 167 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxyy = cbuffer.data(fg_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxyz = cbuffer.data(fg_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxzz = cbuffer.data(fg_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xyyy = cbuffer.data(fg_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xyyz = cbuffer.data(fg_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xyzz = cbuffer.data(fg_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xzzz = cbuffer.data(fg_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yyyy = cbuffer.data(fg_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yyyz = cbuffer.data(fg_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yyzz = cbuffer.data(fg_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yzzz = cbuffer.data(fg_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_zzzz = cbuffer.data(fg_geom_1010_off + 179 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxxx = cbuffer.data(fg_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxxy = cbuffer.data(fg_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxxz = cbuffer.data(fg_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxyy = cbuffer.data(fg_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxyz = cbuffer.data(fg_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxzz = cbuffer.data(fg_geom_1010_off + 185 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xyyy = cbuffer.data(fg_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xyyz = cbuffer.data(fg_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xyzz = cbuffer.data(fg_geom_1010_off + 188 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xzzz = cbuffer.data(fg_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yyyy = cbuffer.data(fg_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yyyz = cbuffer.data(fg_geom_1010_off + 191 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yyzz = cbuffer.data(fg_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yzzz = cbuffer.data(fg_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_zzzz = cbuffer.data(fg_geom_1010_off + 194 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxxx = cbuffer.data(fg_geom_1010_off + 195 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxxy = cbuffer.data(fg_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxxz = cbuffer.data(fg_geom_1010_off + 197 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxyy = cbuffer.data(fg_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxyz = cbuffer.data(fg_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxzz = cbuffer.data(fg_geom_1010_off + 200 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xyyy = cbuffer.data(fg_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xyyz = cbuffer.data(fg_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xyzz = cbuffer.data(fg_geom_1010_off + 203 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xzzz = cbuffer.data(fg_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yyyy = cbuffer.data(fg_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yyyz = cbuffer.data(fg_geom_1010_off + 206 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yyzz = cbuffer.data(fg_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yzzz = cbuffer.data(fg_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_zzzz = cbuffer.data(fg_geom_1010_off + 209 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxxx = cbuffer.data(fg_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxxy = cbuffer.data(fg_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxxz = cbuffer.data(fg_geom_1010_off + 212 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxyy = cbuffer.data(fg_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxyz = cbuffer.data(fg_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxzz = cbuffer.data(fg_geom_1010_off + 215 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xyyy = cbuffer.data(fg_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xyyz = cbuffer.data(fg_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xyzz = cbuffer.data(fg_geom_1010_off + 218 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xzzz = cbuffer.data(fg_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yyyy = cbuffer.data(fg_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yyyz = cbuffer.data(fg_geom_1010_off + 221 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yyzz = cbuffer.data(fg_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yzzz = cbuffer.data(fg_geom_1010_off + 223 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_zzzz = cbuffer.data(fg_geom_1010_off + 224 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxxx = cbuffer.data(fg_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxxy = cbuffer.data(fg_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxxz = cbuffer.data(fg_geom_1010_off + 227 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxyy = cbuffer.data(fg_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxyz = cbuffer.data(fg_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxzz = cbuffer.data(fg_geom_1010_off + 230 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xyyy = cbuffer.data(fg_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xyyz = cbuffer.data(fg_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xyzz = cbuffer.data(fg_geom_1010_off + 233 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xzzz = cbuffer.data(fg_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yyyy = cbuffer.data(fg_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yyyz = cbuffer.data(fg_geom_1010_off + 236 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yyzz = cbuffer.data(fg_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yzzz = cbuffer.data(fg_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_zzzz = cbuffer.data(fg_geom_1010_off + 239 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxxx = cbuffer.data(fg_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxxy = cbuffer.data(fg_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxxz = cbuffer.data(fg_geom_1010_off + 242 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxyy = cbuffer.data(fg_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxyz = cbuffer.data(fg_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxzz = cbuffer.data(fg_geom_1010_off + 245 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xyyy = cbuffer.data(fg_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xyyz = cbuffer.data(fg_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xyzz = cbuffer.data(fg_geom_1010_off + 248 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xzzz = cbuffer.data(fg_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yyyy = cbuffer.data(fg_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yyyz = cbuffer.data(fg_geom_1010_off + 251 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yyzz = cbuffer.data(fg_geom_1010_off + 252 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yzzz = cbuffer.data(fg_geom_1010_off + 253 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_zzzz = cbuffer.data(fg_geom_1010_off + 254 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxxx = cbuffer.data(fg_geom_1010_off + 255 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxxy = cbuffer.data(fg_geom_1010_off + 256 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxxz = cbuffer.data(fg_geom_1010_off + 257 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxyy = cbuffer.data(fg_geom_1010_off + 258 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxyz = cbuffer.data(fg_geom_1010_off + 259 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxzz = cbuffer.data(fg_geom_1010_off + 260 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xyyy = cbuffer.data(fg_geom_1010_off + 261 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xyyz = cbuffer.data(fg_geom_1010_off + 262 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xyzz = cbuffer.data(fg_geom_1010_off + 263 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xzzz = cbuffer.data(fg_geom_1010_off + 264 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yyyy = cbuffer.data(fg_geom_1010_off + 265 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yyyz = cbuffer.data(fg_geom_1010_off + 266 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yyzz = cbuffer.data(fg_geom_1010_off + 267 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yzzz = cbuffer.data(fg_geom_1010_off + 268 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_zzzz = cbuffer.data(fg_geom_1010_off + 269 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxxx = cbuffer.data(fg_geom_1010_off + 270 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxxy = cbuffer.data(fg_geom_1010_off + 271 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxxz = cbuffer.data(fg_geom_1010_off + 272 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxyy = cbuffer.data(fg_geom_1010_off + 273 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxyz = cbuffer.data(fg_geom_1010_off + 274 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxzz = cbuffer.data(fg_geom_1010_off + 275 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xyyy = cbuffer.data(fg_geom_1010_off + 276 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xyyz = cbuffer.data(fg_geom_1010_off + 277 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xyzz = cbuffer.data(fg_geom_1010_off + 278 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xzzz = cbuffer.data(fg_geom_1010_off + 279 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yyyy = cbuffer.data(fg_geom_1010_off + 280 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yyyz = cbuffer.data(fg_geom_1010_off + 281 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yyzz = cbuffer.data(fg_geom_1010_off + 282 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yzzz = cbuffer.data(fg_geom_1010_off + 283 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_zzzz = cbuffer.data(fg_geom_1010_off + 284 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxxx = cbuffer.data(fg_geom_1010_off + 285 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxxy = cbuffer.data(fg_geom_1010_off + 286 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxxz = cbuffer.data(fg_geom_1010_off + 287 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxyy = cbuffer.data(fg_geom_1010_off + 288 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxyz = cbuffer.data(fg_geom_1010_off + 289 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxzz = cbuffer.data(fg_geom_1010_off + 290 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xyyy = cbuffer.data(fg_geom_1010_off + 291 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xyyz = cbuffer.data(fg_geom_1010_off + 292 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xyzz = cbuffer.data(fg_geom_1010_off + 293 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xzzz = cbuffer.data(fg_geom_1010_off + 294 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yyyy = cbuffer.data(fg_geom_1010_off + 295 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yyyz = cbuffer.data(fg_geom_1010_off + 296 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yyzz = cbuffer.data(fg_geom_1010_off + 297 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yzzz = cbuffer.data(fg_geom_1010_off + 298 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_zzzz = cbuffer.data(fg_geom_1010_off + 299 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxxx = cbuffer.data(fg_geom_1010_off + 300 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxxy = cbuffer.data(fg_geom_1010_off + 301 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxxz = cbuffer.data(fg_geom_1010_off + 302 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxyy = cbuffer.data(fg_geom_1010_off + 303 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxyz = cbuffer.data(fg_geom_1010_off + 304 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxzz = cbuffer.data(fg_geom_1010_off + 305 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xyyy = cbuffer.data(fg_geom_1010_off + 306 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xyyz = cbuffer.data(fg_geom_1010_off + 307 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xyzz = cbuffer.data(fg_geom_1010_off + 308 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xzzz = cbuffer.data(fg_geom_1010_off + 309 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yyyy = cbuffer.data(fg_geom_1010_off + 310 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yyyz = cbuffer.data(fg_geom_1010_off + 311 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yyzz = cbuffer.data(fg_geom_1010_off + 312 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yzzz = cbuffer.data(fg_geom_1010_off + 313 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_zzzz = cbuffer.data(fg_geom_1010_off + 314 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxxx = cbuffer.data(fg_geom_1010_off + 315 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxxy = cbuffer.data(fg_geom_1010_off + 316 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxxz = cbuffer.data(fg_geom_1010_off + 317 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxyy = cbuffer.data(fg_geom_1010_off + 318 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxyz = cbuffer.data(fg_geom_1010_off + 319 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxzz = cbuffer.data(fg_geom_1010_off + 320 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xyyy = cbuffer.data(fg_geom_1010_off + 321 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xyyz = cbuffer.data(fg_geom_1010_off + 322 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xyzz = cbuffer.data(fg_geom_1010_off + 323 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xzzz = cbuffer.data(fg_geom_1010_off + 324 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yyyy = cbuffer.data(fg_geom_1010_off + 325 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yyyz = cbuffer.data(fg_geom_1010_off + 326 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yyzz = cbuffer.data(fg_geom_1010_off + 327 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yzzz = cbuffer.data(fg_geom_1010_off + 328 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_zzzz = cbuffer.data(fg_geom_1010_off + 329 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxxx = cbuffer.data(fg_geom_1010_off + 330 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxxy = cbuffer.data(fg_geom_1010_off + 331 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxxz = cbuffer.data(fg_geom_1010_off + 332 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxyy = cbuffer.data(fg_geom_1010_off + 333 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxyz = cbuffer.data(fg_geom_1010_off + 334 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxzz = cbuffer.data(fg_geom_1010_off + 335 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xyyy = cbuffer.data(fg_geom_1010_off + 336 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xyyz = cbuffer.data(fg_geom_1010_off + 337 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xyzz = cbuffer.data(fg_geom_1010_off + 338 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xzzz = cbuffer.data(fg_geom_1010_off + 339 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yyyy = cbuffer.data(fg_geom_1010_off + 340 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yyyz = cbuffer.data(fg_geom_1010_off + 341 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yyzz = cbuffer.data(fg_geom_1010_off + 342 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yzzz = cbuffer.data(fg_geom_1010_off + 343 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_zzzz = cbuffer.data(fg_geom_1010_off + 344 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxxx = cbuffer.data(fg_geom_1010_off + 345 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxxy = cbuffer.data(fg_geom_1010_off + 346 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxxz = cbuffer.data(fg_geom_1010_off + 347 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxyy = cbuffer.data(fg_geom_1010_off + 348 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxyz = cbuffer.data(fg_geom_1010_off + 349 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxzz = cbuffer.data(fg_geom_1010_off + 350 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xyyy = cbuffer.data(fg_geom_1010_off + 351 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xyyz = cbuffer.data(fg_geom_1010_off + 352 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xyzz = cbuffer.data(fg_geom_1010_off + 353 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xzzz = cbuffer.data(fg_geom_1010_off + 354 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yyyy = cbuffer.data(fg_geom_1010_off + 355 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yyyz = cbuffer.data(fg_geom_1010_off + 356 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yyzz = cbuffer.data(fg_geom_1010_off + 357 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yzzz = cbuffer.data(fg_geom_1010_off + 358 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_zzzz = cbuffer.data(fg_geom_1010_off + 359 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxxx = cbuffer.data(fg_geom_1010_off + 360 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxxy = cbuffer.data(fg_geom_1010_off + 361 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxxz = cbuffer.data(fg_geom_1010_off + 362 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxyy = cbuffer.data(fg_geom_1010_off + 363 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxyz = cbuffer.data(fg_geom_1010_off + 364 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxzz = cbuffer.data(fg_geom_1010_off + 365 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xyyy = cbuffer.data(fg_geom_1010_off + 366 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xyyz = cbuffer.data(fg_geom_1010_off + 367 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xyzz = cbuffer.data(fg_geom_1010_off + 368 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xzzz = cbuffer.data(fg_geom_1010_off + 369 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yyyy = cbuffer.data(fg_geom_1010_off + 370 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yyyz = cbuffer.data(fg_geom_1010_off + 371 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yyzz = cbuffer.data(fg_geom_1010_off + 372 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yzzz = cbuffer.data(fg_geom_1010_off + 373 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_zzzz = cbuffer.data(fg_geom_1010_off + 374 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxxx = cbuffer.data(fg_geom_1010_off + 375 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxxy = cbuffer.data(fg_geom_1010_off + 376 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxxz = cbuffer.data(fg_geom_1010_off + 377 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxyy = cbuffer.data(fg_geom_1010_off + 378 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxyz = cbuffer.data(fg_geom_1010_off + 379 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxzz = cbuffer.data(fg_geom_1010_off + 380 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xyyy = cbuffer.data(fg_geom_1010_off + 381 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xyyz = cbuffer.data(fg_geom_1010_off + 382 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xyzz = cbuffer.data(fg_geom_1010_off + 383 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xzzz = cbuffer.data(fg_geom_1010_off + 384 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yyyy = cbuffer.data(fg_geom_1010_off + 385 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yyyz = cbuffer.data(fg_geom_1010_off + 386 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yyzz = cbuffer.data(fg_geom_1010_off + 387 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yzzz = cbuffer.data(fg_geom_1010_off + 388 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_zzzz = cbuffer.data(fg_geom_1010_off + 389 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxxx = cbuffer.data(fg_geom_1010_off + 390 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxxy = cbuffer.data(fg_geom_1010_off + 391 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxxz = cbuffer.data(fg_geom_1010_off + 392 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxyy = cbuffer.data(fg_geom_1010_off + 393 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxyz = cbuffer.data(fg_geom_1010_off + 394 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxzz = cbuffer.data(fg_geom_1010_off + 395 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xyyy = cbuffer.data(fg_geom_1010_off + 396 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xyyz = cbuffer.data(fg_geom_1010_off + 397 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xyzz = cbuffer.data(fg_geom_1010_off + 398 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xzzz = cbuffer.data(fg_geom_1010_off + 399 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yyyy = cbuffer.data(fg_geom_1010_off + 400 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yyyz = cbuffer.data(fg_geom_1010_off + 401 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yyzz = cbuffer.data(fg_geom_1010_off + 402 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yzzz = cbuffer.data(fg_geom_1010_off + 403 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_zzzz = cbuffer.data(fg_geom_1010_off + 404 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxxx = cbuffer.data(fg_geom_1010_off + 405 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxxy = cbuffer.data(fg_geom_1010_off + 406 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxxz = cbuffer.data(fg_geom_1010_off + 407 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxyy = cbuffer.data(fg_geom_1010_off + 408 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxyz = cbuffer.data(fg_geom_1010_off + 409 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxzz = cbuffer.data(fg_geom_1010_off + 410 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xyyy = cbuffer.data(fg_geom_1010_off + 411 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xyyz = cbuffer.data(fg_geom_1010_off + 412 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xyzz = cbuffer.data(fg_geom_1010_off + 413 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xzzz = cbuffer.data(fg_geom_1010_off + 414 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yyyy = cbuffer.data(fg_geom_1010_off + 415 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yyyz = cbuffer.data(fg_geom_1010_off + 416 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yyzz = cbuffer.data(fg_geom_1010_off + 417 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yzzz = cbuffer.data(fg_geom_1010_off + 418 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_zzzz = cbuffer.data(fg_geom_1010_off + 419 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxxx = cbuffer.data(fg_geom_1010_off + 420 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxxy = cbuffer.data(fg_geom_1010_off + 421 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxxz = cbuffer.data(fg_geom_1010_off + 422 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxyy = cbuffer.data(fg_geom_1010_off + 423 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxyz = cbuffer.data(fg_geom_1010_off + 424 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxzz = cbuffer.data(fg_geom_1010_off + 425 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xyyy = cbuffer.data(fg_geom_1010_off + 426 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xyyz = cbuffer.data(fg_geom_1010_off + 427 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xyzz = cbuffer.data(fg_geom_1010_off + 428 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xzzz = cbuffer.data(fg_geom_1010_off + 429 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yyyy = cbuffer.data(fg_geom_1010_off + 430 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yyyz = cbuffer.data(fg_geom_1010_off + 431 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yyzz = cbuffer.data(fg_geom_1010_off + 432 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yzzz = cbuffer.data(fg_geom_1010_off + 433 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_zzzz = cbuffer.data(fg_geom_1010_off + 434 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxxx = cbuffer.data(fg_geom_1010_off + 435 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxxy = cbuffer.data(fg_geom_1010_off + 436 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxxz = cbuffer.data(fg_geom_1010_off + 437 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxyy = cbuffer.data(fg_geom_1010_off + 438 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxyz = cbuffer.data(fg_geom_1010_off + 439 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxzz = cbuffer.data(fg_geom_1010_off + 440 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xyyy = cbuffer.data(fg_geom_1010_off + 441 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xyyz = cbuffer.data(fg_geom_1010_off + 442 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xyzz = cbuffer.data(fg_geom_1010_off + 443 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xzzz = cbuffer.data(fg_geom_1010_off + 444 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yyyy = cbuffer.data(fg_geom_1010_off + 445 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yyyz = cbuffer.data(fg_geom_1010_off + 446 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yyzz = cbuffer.data(fg_geom_1010_off + 447 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yzzz = cbuffer.data(fg_geom_1010_off + 448 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_zzzz = cbuffer.data(fg_geom_1010_off + 449 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxxx = cbuffer.data(fg_geom_1010_off + 450 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxxy = cbuffer.data(fg_geom_1010_off + 451 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxxz = cbuffer.data(fg_geom_1010_off + 452 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxyy = cbuffer.data(fg_geom_1010_off + 453 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxyz = cbuffer.data(fg_geom_1010_off + 454 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxzz = cbuffer.data(fg_geom_1010_off + 455 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xyyy = cbuffer.data(fg_geom_1010_off + 456 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xyyz = cbuffer.data(fg_geom_1010_off + 457 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xyzz = cbuffer.data(fg_geom_1010_off + 458 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xzzz = cbuffer.data(fg_geom_1010_off + 459 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yyyy = cbuffer.data(fg_geom_1010_off + 460 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yyyz = cbuffer.data(fg_geom_1010_off + 461 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yyzz = cbuffer.data(fg_geom_1010_off + 462 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yzzz = cbuffer.data(fg_geom_1010_off + 463 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_zzzz = cbuffer.data(fg_geom_1010_off + 464 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxxx = cbuffer.data(fg_geom_1010_off + 465 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxxy = cbuffer.data(fg_geom_1010_off + 466 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxxz = cbuffer.data(fg_geom_1010_off + 467 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxyy = cbuffer.data(fg_geom_1010_off + 468 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxyz = cbuffer.data(fg_geom_1010_off + 469 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxzz = cbuffer.data(fg_geom_1010_off + 470 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xyyy = cbuffer.data(fg_geom_1010_off + 471 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xyyz = cbuffer.data(fg_geom_1010_off + 472 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xyzz = cbuffer.data(fg_geom_1010_off + 473 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xzzz = cbuffer.data(fg_geom_1010_off + 474 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yyyy = cbuffer.data(fg_geom_1010_off + 475 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yyyz = cbuffer.data(fg_geom_1010_off + 476 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yyzz = cbuffer.data(fg_geom_1010_off + 477 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yzzz = cbuffer.data(fg_geom_1010_off + 478 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_zzzz = cbuffer.data(fg_geom_1010_off + 479 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxxx = cbuffer.data(fg_geom_1010_off + 480 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxxy = cbuffer.data(fg_geom_1010_off + 481 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxxz = cbuffer.data(fg_geom_1010_off + 482 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxyy = cbuffer.data(fg_geom_1010_off + 483 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxyz = cbuffer.data(fg_geom_1010_off + 484 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxzz = cbuffer.data(fg_geom_1010_off + 485 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xyyy = cbuffer.data(fg_geom_1010_off + 486 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xyyz = cbuffer.data(fg_geom_1010_off + 487 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xyzz = cbuffer.data(fg_geom_1010_off + 488 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xzzz = cbuffer.data(fg_geom_1010_off + 489 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yyyy = cbuffer.data(fg_geom_1010_off + 490 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yyyz = cbuffer.data(fg_geom_1010_off + 491 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yyzz = cbuffer.data(fg_geom_1010_off + 492 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yzzz = cbuffer.data(fg_geom_1010_off + 493 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_zzzz = cbuffer.data(fg_geom_1010_off + 494 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxxx = cbuffer.data(fg_geom_1010_off + 495 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxxy = cbuffer.data(fg_geom_1010_off + 496 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxxz = cbuffer.data(fg_geom_1010_off + 497 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxyy = cbuffer.data(fg_geom_1010_off + 498 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxyz = cbuffer.data(fg_geom_1010_off + 499 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxzz = cbuffer.data(fg_geom_1010_off + 500 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xyyy = cbuffer.data(fg_geom_1010_off + 501 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xyyz = cbuffer.data(fg_geom_1010_off + 502 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xyzz = cbuffer.data(fg_geom_1010_off + 503 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xzzz = cbuffer.data(fg_geom_1010_off + 504 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yyyy = cbuffer.data(fg_geom_1010_off + 505 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yyyz = cbuffer.data(fg_geom_1010_off + 506 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yyzz = cbuffer.data(fg_geom_1010_off + 507 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yzzz = cbuffer.data(fg_geom_1010_off + 508 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_zzzz = cbuffer.data(fg_geom_1010_off + 509 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxxx = cbuffer.data(fg_geom_1010_off + 510 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxxy = cbuffer.data(fg_geom_1010_off + 511 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxxz = cbuffer.data(fg_geom_1010_off + 512 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxyy = cbuffer.data(fg_geom_1010_off + 513 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxyz = cbuffer.data(fg_geom_1010_off + 514 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxzz = cbuffer.data(fg_geom_1010_off + 515 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xyyy = cbuffer.data(fg_geom_1010_off + 516 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xyyz = cbuffer.data(fg_geom_1010_off + 517 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xyzz = cbuffer.data(fg_geom_1010_off + 518 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xzzz = cbuffer.data(fg_geom_1010_off + 519 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yyyy = cbuffer.data(fg_geom_1010_off + 520 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yyyz = cbuffer.data(fg_geom_1010_off + 521 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yyzz = cbuffer.data(fg_geom_1010_off + 522 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yzzz = cbuffer.data(fg_geom_1010_off + 523 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_zzzz = cbuffer.data(fg_geom_1010_off + 524 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxxx = cbuffer.data(fg_geom_1010_off + 525 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxxy = cbuffer.data(fg_geom_1010_off + 526 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxxz = cbuffer.data(fg_geom_1010_off + 527 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxyy = cbuffer.data(fg_geom_1010_off + 528 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxyz = cbuffer.data(fg_geom_1010_off + 529 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxzz = cbuffer.data(fg_geom_1010_off + 530 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xyyy = cbuffer.data(fg_geom_1010_off + 531 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xyyz = cbuffer.data(fg_geom_1010_off + 532 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xyzz = cbuffer.data(fg_geom_1010_off + 533 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xzzz = cbuffer.data(fg_geom_1010_off + 534 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yyyy = cbuffer.data(fg_geom_1010_off + 535 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yyyz = cbuffer.data(fg_geom_1010_off + 536 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yyzz = cbuffer.data(fg_geom_1010_off + 537 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yzzz = cbuffer.data(fg_geom_1010_off + 538 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_zzzz = cbuffer.data(fg_geom_1010_off + 539 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxxx = cbuffer.data(fg_geom_1010_off + 540 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxxy = cbuffer.data(fg_geom_1010_off + 541 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxxz = cbuffer.data(fg_geom_1010_off + 542 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxyy = cbuffer.data(fg_geom_1010_off + 543 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxyz = cbuffer.data(fg_geom_1010_off + 544 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxzz = cbuffer.data(fg_geom_1010_off + 545 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xyyy = cbuffer.data(fg_geom_1010_off + 546 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xyyz = cbuffer.data(fg_geom_1010_off + 547 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xyzz = cbuffer.data(fg_geom_1010_off + 548 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xzzz = cbuffer.data(fg_geom_1010_off + 549 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yyyy = cbuffer.data(fg_geom_1010_off + 550 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yyyz = cbuffer.data(fg_geom_1010_off + 551 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yyzz = cbuffer.data(fg_geom_1010_off + 552 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yzzz = cbuffer.data(fg_geom_1010_off + 553 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_zzzz = cbuffer.data(fg_geom_1010_off + 554 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxxx = cbuffer.data(fg_geom_1010_off + 555 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxxy = cbuffer.data(fg_geom_1010_off + 556 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxxz = cbuffer.data(fg_geom_1010_off + 557 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxyy = cbuffer.data(fg_geom_1010_off + 558 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxyz = cbuffer.data(fg_geom_1010_off + 559 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxzz = cbuffer.data(fg_geom_1010_off + 560 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xyyy = cbuffer.data(fg_geom_1010_off + 561 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xyyz = cbuffer.data(fg_geom_1010_off + 562 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xyzz = cbuffer.data(fg_geom_1010_off + 563 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xzzz = cbuffer.data(fg_geom_1010_off + 564 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yyyy = cbuffer.data(fg_geom_1010_off + 565 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yyyz = cbuffer.data(fg_geom_1010_off + 566 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yyzz = cbuffer.data(fg_geom_1010_off + 567 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yzzz = cbuffer.data(fg_geom_1010_off + 568 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_zzzz = cbuffer.data(fg_geom_1010_off + 569 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxxx = cbuffer.data(fg_geom_1010_off + 570 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxxy = cbuffer.data(fg_geom_1010_off + 571 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxxz = cbuffer.data(fg_geom_1010_off + 572 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxyy = cbuffer.data(fg_geom_1010_off + 573 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxyz = cbuffer.data(fg_geom_1010_off + 574 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxzz = cbuffer.data(fg_geom_1010_off + 575 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xyyy = cbuffer.data(fg_geom_1010_off + 576 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xyyz = cbuffer.data(fg_geom_1010_off + 577 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xyzz = cbuffer.data(fg_geom_1010_off + 578 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xzzz = cbuffer.data(fg_geom_1010_off + 579 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yyyy = cbuffer.data(fg_geom_1010_off + 580 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yyyz = cbuffer.data(fg_geom_1010_off + 581 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yyzz = cbuffer.data(fg_geom_1010_off + 582 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yzzz = cbuffer.data(fg_geom_1010_off + 583 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_zzzz = cbuffer.data(fg_geom_1010_off + 584 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxxx = cbuffer.data(fg_geom_1010_off + 585 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxxy = cbuffer.data(fg_geom_1010_off + 586 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxxz = cbuffer.data(fg_geom_1010_off + 587 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxyy = cbuffer.data(fg_geom_1010_off + 588 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxyz = cbuffer.data(fg_geom_1010_off + 589 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxzz = cbuffer.data(fg_geom_1010_off + 590 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xyyy = cbuffer.data(fg_geom_1010_off + 591 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xyyz = cbuffer.data(fg_geom_1010_off + 592 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xyzz = cbuffer.data(fg_geom_1010_off + 593 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xzzz = cbuffer.data(fg_geom_1010_off + 594 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yyyy = cbuffer.data(fg_geom_1010_off + 595 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yyyz = cbuffer.data(fg_geom_1010_off + 596 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yyzz = cbuffer.data(fg_geom_1010_off + 597 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yzzz = cbuffer.data(fg_geom_1010_off + 598 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_zzzz = cbuffer.data(fg_geom_1010_off + 599 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxxx = cbuffer.data(fg_geom_1010_off + 600 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxxy = cbuffer.data(fg_geom_1010_off + 601 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxxz = cbuffer.data(fg_geom_1010_off + 602 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxyy = cbuffer.data(fg_geom_1010_off + 603 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxyz = cbuffer.data(fg_geom_1010_off + 604 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxzz = cbuffer.data(fg_geom_1010_off + 605 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xyyy = cbuffer.data(fg_geom_1010_off + 606 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xyyz = cbuffer.data(fg_geom_1010_off + 607 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xyzz = cbuffer.data(fg_geom_1010_off + 608 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xzzz = cbuffer.data(fg_geom_1010_off + 609 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yyyy = cbuffer.data(fg_geom_1010_off + 610 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yyyz = cbuffer.data(fg_geom_1010_off + 611 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yyzz = cbuffer.data(fg_geom_1010_off + 612 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yzzz = cbuffer.data(fg_geom_1010_off + 613 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_zzzz = cbuffer.data(fg_geom_1010_off + 614 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxxx = cbuffer.data(fg_geom_1010_off + 615 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxxy = cbuffer.data(fg_geom_1010_off + 616 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxxz = cbuffer.data(fg_geom_1010_off + 617 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxyy = cbuffer.data(fg_geom_1010_off + 618 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxyz = cbuffer.data(fg_geom_1010_off + 619 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxzz = cbuffer.data(fg_geom_1010_off + 620 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xyyy = cbuffer.data(fg_geom_1010_off + 621 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xyyz = cbuffer.data(fg_geom_1010_off + 622 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xyzz = cbuffer.data(fg_geom_1010_off + 623 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xzzz = cbuffer.data(fg_geom_1010_off + 624 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yyyy = cbuffer.data(fg_geom_1010_off + 625 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yyyz = cbuffer.data(fg_geom_1010_off + 626 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yyzz = cbuffer.data(fg_geom_1010_off + 627 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yzzz = cbuffer.data(fg_geom_1010_off + 628 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_zzzz = cbuffer.data(fg_geom_1010_off + 629 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxxx = cbuffer.data(fg_geom_1010_off + 630 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxxy = cbuffer.data(fg_geom_1010_off + 631 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxxz = cbuffer.data(fg_geom_1010_off + 632 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxyy = cbuffer.data(fg_geom_1010_off + 633 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxyz = cbuffer.data(fg_geom_1010_off + 634 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxzz = cbuffer.data(fg_geom_1010_off + 635 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xyyy = cbuffer.data(fg_geom_1010_off + 636 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xyyz = cbuffer.data(fg_geom_1010_off + 637 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xyzz = cbuffer.data(fg_geom_1010_off + 638 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xzzz = cbuffer.data(fg_geom_1010_off + 639 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yyyy = cbuffer.data(fg_geom_1010_off + 640 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yyyz = cbuffer.data(fg_geom_1010_off + 641 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yyzz = cbuffer.data(fg_geom_1010_off + 642 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yzzz = cbuffer.data(fg_geom_1010_off + 643 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_zzzz = cbuffer.data(fg_geom_1010_off + 644 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxxx = cbuffer.data(fg_geom_1010_off + 645 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxxy = cbuffer.data(fg_geom_1010_off + 646 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxxz = cbuffer.data(fg_geom_1010_off + 647 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxyy = cbuffer.data(fg_geom_1010_off + 648 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxyz = cbuffer.data(fg_geom_1010_off + 649 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxzz = cbuffer.data(fg_geom_1010_off + 650 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xyyy = cbuffer.data(fg_geom_1010_off + 651 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xyyz = cbuffer.data(fg_geom_1010_off + 652 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xyzz = cbuffer.data(fg_geom_1010_off + 653 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xzzz = cbuffer.data(fg_geom_1010_off + 654 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yyyy = cbuffer.data(fg_geom_1010_off + 655 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yyyz = cbuffer.data(fg_geom_1010_off + 656 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yyzz = cbuffer.data(fg_geom_1010_off + 657 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yzzz = cbuffer.data(fg_geom_1010_off + 658 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_zzzz = cbuffer.data(fg_geom_1010_off + 659 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxxx = cbuffer.data(fg_geom_1010_off + 660 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxxy = cbuffer.data(fg_geom_1010_off + 661 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxxz = cbuffer.data(fg_geom_1010_off + 662 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxyy = cbuffer.data(fg_geom_1010_off + 663 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxyz = cbuffer.data(fg_geom_1010_off + 664 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxzz = cbuffer.data(fg_geom_1010_off + 665 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xyyy = cbuffer.data(fg_geom_1010_off + 666 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xyyz = cbuffer.data(fg_geom_1010_off + 667 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xyzz = cbuffer.data(fg_geom_1010_off + 668 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xzzz = cbuffer.data(fg_geom_1010_off + 669 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yyyy = cbuffer.data(fg_geom_1010_off + 670 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yyyz = cbuffer.data(fg_geom_1010_off + 671 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yyzz = cbuffer.data(fg_geom_1010_off + 672 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yzzz = cbuffer.data(fg_geom_1010_off + 673 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_zzzz = cbuffer.data(fg_geom_1010_off + 674 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxxx = cbuffer.data(fg_geom_1010_off + 675 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxxy = cbuffer.data(fg_geom_1010_off + 676 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxxz = cbuffer.data(fg_geom_1010_off + 677 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxyy = cbuffer.data(fg_geom_1010_off + 678 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxyz = cbuffer.data(fg_geom_1010_off + 679 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxzz = cbuffer.data(fg_geom_1010_off + 680 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xyyy = cbuffer.data(fg_geom_1010_off + 681 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xyyz = cbuffer.data(fg_geom_1010_off + 682 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xyzz = cbuffer.data(fg_geom_1010_off + 683 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xzzz = cbuffer.data(fg_geom_1010_off + 684 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yyyy = cbuffer.data(fg_geom_1010_off + 685 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yyyz = cbuffer.data(fg_geom_1010_off + 686 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yyzz = cbuffer.data(fg_geom_1010_off + 687 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yzzz = cbuffer.data(fg_geom_1010_off + 688 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_zzzz = cbuffer.data(fg_geom_1010_off + 689 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxxx = cbuffer.data(fg_geom_1010_off + 690 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxxy = cbuffer.data(fg_geom_1010_off + 691 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxxz = cbuffer.data(fg_geom_1010_off + 692 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxyy = cbuffer.data(fg_geom_1010_off + 693 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxyz = cbuffer.data(fg_geom_1010_off + 694 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxzz = cbuffer.data(fg_geom_1010_off + 695 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xyyy = cbuffer.data(fg_geom_1010_off + 696 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xyyz = cbuffer.data(fg_geom_1010_off + 697 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xyzz = cbuffer.data(fg_geom_1010_off + 698 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xzzz = cbuffer.data(fg_geom_1010_off + 699 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yyyy = cbuffer.data(fg_geom_1010_off + 700 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yyyz = cbuffer.data(fg_geom_1010_off + 701 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yyzz = cbuffer.data(fg_geom_1010_off + 702 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yzzz = cbuffer.data(fg_geom_1010_off + 703 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_zzzz = cbuffer.data(fg_geom_1010_off + 704 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxxx = cbuffer.data(fg_geom_1010_off + 705 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxxy = cbuffer.data(fg_geom_1010_off + 706 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxxz = cbuffer.data(fg_geom_1010_off + 707 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxyy = cbuffer.data(fg_geom_1010_off + 708 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxyz = cbuffer.data(fg_geom_1010_off + 709 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxzz = cbuffer.data(fg_geom_1010_off + 710 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xyyy = cbuffer.data(fg_geom_1010_off + 711 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xyyz = cbuffer.data(fg_geom_1010_off + 712 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xyzz = cbuffer.data(fg_geom_1010_off + 713 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xzzz = cbuffer.data(fg_geom_1010_off + 714 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yyyy = cbuffer.data(fg_geom_1010_off + 715 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yyyz = cbuffer.data(fg_geom_1010_off + 716 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yyzz = cbuffer.data(fg_geom_1010_off + 717 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yzzz = cbuffer.data(fg_geom_1010_off + 718 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_zzzz = cbuffer.data(fg_geom_1010_off + 719 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxxx = cbuffer.data(fg_geom_1010_off + 720 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxxy = cbuffer.data(fg_geom_1010_off + 721 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxxz = cbuffer.data(fg_geom_1010_off + 722 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxyy = cbuffer.data(fg_geom_1010_off + 723 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxyz = cbuffer.data(fg_geom_1010_off + 724 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxzz = cbuffer.data(fg_geom_1010_off + 725 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xyyy = cbuffer.data(fg_geom_1010_off + 726 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xyyz = cbuffer.data(fg_geom_1010_off + 727 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xyzz = cbuffer.data(fg_geom_1010_off + 728 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xzzz = cbuffer.data(fg_geom_1010_off + 729 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yyyy = cbuffer.data(fg_geom_1010_off + 730 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yyyz = cbuffer.data(fg_geom_1010_off + 731 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yyzz = cbuffer.data(fg_geom_1010_off + 732 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yzzz = cbuffer.data(fg_geom_1010_off + 733 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_zzzz = cbuffer.data(fg_geom_1010_off + 734 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxxx = cbuffer.data(fg_geom_1010_off + 735 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxxy = cbuffer.data(fg_geom_1010_off + 736 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxxz = cbuffer.data(fg_geom_1010_off + 737 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxyy = cbuffer.data(fg_geom_1010_off + 738 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxyz = cbuffer.data(fg_geom_1010_off + 739 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxzz = cbuffer.data(fg_geom_1010_off + 740 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xyyy = cbuffer.data(fg_geom_1010_off + 741 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xyyz = cbuffer.data(fg_geom_1010_off + 742 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xyzz = cbuffer.data(fg_geom_1010_off + 743 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xzzz = cbuffer.data(fg_geom_1010_off + 744 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yyyy = cbuffer.data(fg_geom_1010_off + 745 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yyyz = cbuffer.data(fg_geom_1010_off + 746 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yyzz = cbuffer.data(fg_geom_1010_off + 747 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yzzz = cbuffer.data(fg_geom_1010_off + 748 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_zzzz = cbuffer.data(fg_geom_1010_off + 749 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxxx = cbuffer.data(fg_geom_1010_off + 750 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxxy = cbuffer.data(fg_geom_1010_off + 751 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxxz = cbuffer.data(fg_geom_1010_off + 752 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxyy = cbuffer.data(fg_geom_1010_off + 753 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxyz = cbuffer.data(fg_geom_1010_off + 754 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxzz = cbuffer.data(fg_geom_1010_off + 755 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xyyy = cbuffer.data(fg_geom_1010_off + 756 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xyyz = cbuffer.data(fg_geom_1010_off + 757 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xyzz = cbuffer.data(fg_geom_1010_off + 758 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xzzz = cbuffer.data(fg_geom_1010_off + 759 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yyyy = cbuffer.data(fg_geom_1010_off + 760 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yyyz = cbuffer.data(fg_geom_1010_off + 761 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yyzz = cbuffer.data(fg_geom_1010_off + 762 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yzzz = cbuffer.data(fg_geom_1010_off + 763 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_zzzz = cbuffer.data(fg_geom_1010_off + 764 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxxx = cbuffer.data(fg_geom_1010_off + 765 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxxy = cbuffer.data(fg_geom_1010_off + 766 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxxz = cbuffer.data(fg_geom_1010_off + 767 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxyy = cbuffer.data(fg_geom_1010_off + 768 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxyz = cbuffer.data(fg_geom_1010_off + 769 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxzz = cbuffer.data(fg_geom_1010_off + 770 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xyyy = cbuffer.data(fg_geom_1010_off + 771 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xyyz = cbuffer.data(fg_geom_1010_off + 772 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xyzz = cbuffer.data(fg_geom_1010_off + 773 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xzzz = cbuffer.data(fg_geom_1010_off + 774 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yyyy = cbuffer.data(fg_geom_1010_off + 775 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yyyz = cbuffer.data(fg_geom_1010_off + 776 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yyzz = cbuffer.data(fg_geom_1010_off + 777 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yzzz = cbuffer.data(fg_geom_1010_off + 778 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_zzzz = cbuffer.data(fg_geom_1010_off + 779 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxxx = cbuffer.data(fg_geom_1010_off + 780 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxxy = cbuffer.data(fg_geom_1010_off + 781 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxxz = cbuffer.data(fg_geom_1010_off + 782 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxyy = cbuffer.data(fg_geom_1010_off + 783 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxyz = cbuffer.data(fg_geom_1010_off + 784 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxzz = cbuffer.data(fg_geom_1010_off + 785 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xyyy = cbuffer.data(fg_geom_1010_off + 786 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xyyz = cbuffer.data(fg_geom_1010_off + 787 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xyzz = cbuffer.data(fg_geom_1010_off + 788 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xzzz = cbuffer.data(fg_geom_1010_off + 789 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yyyy = cbuffer.data(fg_geom_1010_off + 790 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yyyz = cbuffer.data(fg_geom_1010_off + 791 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yyzz = cbuffer.data(fg_geom_1010_off + 792 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yzzz = cbuffer.data(fg_geom_1010_off + 793 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_zzzz = cbuffer.data(fg_geom_1010_off + 794 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxxx = cbuffer.data(fg_geom_1010_off + 795 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxxy = cbuffer.data(fg_geom_1010_off + 796 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxxz = cbuffer.data(fg_geom_1010_off + 797 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxyy = cbuffer.data(fg_geom_1010_off + 798 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxyz = cbuffer.data(fg_geom_1010_off + 799 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxzz = cbuffer.data(fg_geom_1010_off + 800 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xyyy = cbuffer.data(fg_geom_1010_off + 801 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xyyz = cbuffer.data(fg_geom_1010_off + 802 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xyzz = cbuffer.data(fg_geom_1010_off + 803 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xzzz = cbuffer.data(fg_geom_1010_off + 804 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yyyy = cbuffer.data(fg_geom_1010_off + 805 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yyyz = cbuffer.data(fg_geom_1010_off + 806 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yyzz = cbuffer.data(fg_geom_1010_off + 807 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yzzz = cbuffer.data(fg_geom_1010_off + 808 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_zzzz = cbuffer.data(fg_geom_1010_off + 809 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxxx = cbuffer.data(fg_geom_1010_off + 810 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxxy = cbuffer.data(fg_geom_1010_off + 811 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxxz = cbuffer.data(fg_geom_1010_off + 812 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxyy = cbuffer.data(fg_geom_1010_off + 813 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxyz = cbuffer.data(fg_geom_1010_off + 814 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxzz = cbuffer.data(fg_geom_1010_off + 815 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xyyy = cbuffer.data(fg_geom_1010_off + 816 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xyyz = cbuffer.data(fg_geom_1010_off + 817 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xyzz = cbuffer.data(fg_geom_1010_off + 818 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xzzz = cbuffer.data(fg_geom_1010_off + 819 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yyyy = cbuffer.data(fg_geom_1010_off + 820 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yyyz = cbuffer.data(fg_geom_1010_off + 821 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yyzz = cbuffer.data(fg_geom_1010_off + 822 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yzzz = cbuffer.data(fg_geom_1010_off + 823 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_zzzz = cbuffer.data(fg_geom_1010_off + 824 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxxx = cbuffer.data(fg_geom_1010_off + 825 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxxy = cbuffer.data(fg_geom_1010_off + 826 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxxz = cbuffer.data(fg_geom_1010_off + 827 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxyy = cbuffer.data(fg_geom_1010_off + 828 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxyz = cbuffer.data(fg_geom_1010_off + 829 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxzz = cbuffer.data(fg_geom_1010_off + 830 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xyyy = cbuffer.data(fg_geom_1010_off + 831 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xyyz = cbuffer.data(fg_geom_1010_off + 832 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xyzz = cbuffer.data(fg_geom_1010_off + 833 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xzzz = cbuffer.data(fg_geom_1010_off + 834 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yyyy = cbuffer.data(fg_geom_1010_off + 835 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yyyz = cbuffer.data(fg_geom_1010_off + 836 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yyzz = cbuffer.data(fg_geom_1010_off + 837 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yzzz = cbuffer.data(fg_geom_1010_off + 838 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_zzzz = cbuffer.data(fg_geom_1010_off + 839 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxxx = cbuffer.data(fg_geom_1010_off + 840 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxxy = cbuffer.data(fg_geom_1010_off + 841 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxxz = cbuffer.data(fg_geom_1010_off + 842 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxyy = cbuffer.data(fg_geom_1010_off + 843 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxyz = cbuffer.data(fg_geom_1010_off + 844 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxzz = cbuffer.data(fg_geom_1010_off + 845 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xyyy = cbuffer.data(fg_geom_1010_off + 846 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xyyz = cbuffer.data(fg_geom_1010_off + 847 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xyzz = cbuffer.data(fg_geom_1010_off + 848 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xzzz = cbuffer.data(fg_geom_1010_off + 849 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yyyy = cbuffer.data(fg_geom_1010_off + 850 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yyyz = cbuffer.data(fg_geom_1010_off + 851 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yyzz = cbuffer.data(fg_geom_1010_off + 852 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yzzz = cbuffer.data(fg_geom_1010_off + 853 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_zzzz = cbuffer.data(fg_geom_1010_off + 854 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxxx = cbuffer.data(fg_geom_1010_off + 855 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxxy = cbuffer.data(fg_geom_1010_off + 856 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxxz = cbuffer.data(fg_geom_1010_off + 857 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxyy = cbuffer.data(fg_geom_1010_off + 858 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxyz = cbuffer.data(fg_geom_1010_off + 859 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxzz = cbuffer.data(fg_geom_1010_off + 860 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xyyy = cbuffer.data(fg_geom_1010_off + 861 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xyyz = cbuffer.data(fg_geom_1010_off + 862 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xyzz = cbuffer.data(fg_geom_1010_off + 863 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xzzz = cbuffer.data(fg_geom_1010_off + 864 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yyyy = cbuffer.data(fg_geom_1010_off + 865 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yyyz = cbuffer.data(fg_geom_1010_off + 866 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yyzz = cbuffer.data(fg_geom_1010_off + 867 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yzzz = cbuffer.data(fg_geom_1010_off + 868 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_zzzz = cbuffer.data(fg_geom_1010_off + 869 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxxx = cbuffer.data(fg_geom_1010_off + 870 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxxy = cbuffer.data(fg_geom_1010_off + 871 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxxz = cbuffer.data(fg_geom_1010_off + 872 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxyy = cbuffer.data(fg_geom_1010_off + 873 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxyz = cbuffer.data(fg_geom_1010_off + 874 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxzz = cbuffer.data(fg_geom_1010_off + 875 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xyyy = cbuffer.data(fg_geom_1010_off + 876 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xyyz = cbuffer.data(fg_geom_1010_off + 877 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xyzz = cbuffer.data(fg_geom_1010_off + 878 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xzzz = cbuffer.data(fg_geom_1010_off + 879 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yyyy = cbuffer.data(fg_geom_1010_off + 880 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yyyz = cbuffer.data(fg_geom_1010_off + 881 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yyzz = cbuffer.data(fg_geom_1010_off + 882 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yzzz = cbuffer.data(fg_geom_1010_off + 883 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_zzzz = cbuffer.data(fg_geom_1010_off + 884 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxxx = cbuffer.data(fg_geom_1010_off + 885 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxxy = cbuffer.data(fg_geom_1010_off + 886 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxxz = cbuffer.data(fg_geom_1010_off + 887 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxyy = cbuffer.data(fg_geom_1010_off + 888 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxyz = cbuffer.data(fg_geom_1010_off + 889 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxzz = cbuffer.data(fg_geom_1010_off + 890 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xyyy = cbuffer.data(fg_geom_1010_off + 891 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xyyz = cbuffer.data(fg_geom_1010_off + 892 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xyzz = cbuffer.data(fg_geom_1010_off + 893 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xzzz = cbuffer.data(fg_geom_1010_off + 894 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yyyy = cbuffer.data(fg_geom_1010_off + 895 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yyyz = cbuffer.data(fg_geom_1010_off + 896 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yyzz = cbuffer.data(fg_geom_1010_off + 897 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yzzz = cbuffer.data(fg_geom_1010_off + 898 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_zzzz = cbuffer.data(fg_geom_1010_off + 899 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxxx = cbuffer.data(fg_geom_1010_off + 900 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxxy = cbuffer.data(fg_geom_1010_off + 901 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxxz = cbuffer.data(fg_geom_1010_off + 902 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxyy = cbuffer.data(fg_geom_1010_off + 903 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxyz = cbuffer.data(fg_geom_1010_off + 904 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxzz = cbuffer.data(fg_geom_1010_off + 905 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xyyy = cbuffer.data(fg_geom_1010_off + 906 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xyyz = cbuffer.data(fg_geom_1010_off + 907 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xyzz = cbuffer.data(fg_geom_1010_off + 908 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xzzz = cbuffer.data(fg_geom_1010_off + 909 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yyyy = cbuffer.data(fg_geom_1010_off + 910 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yyyz = cbuffer.data(fg_geom_1010_off + 911 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yyzz = cbuffer.data(fg_geom_1010_off + 912 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yzzz = cbuffer.data(fg_geom_1010_off + 913 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_zzzz = cbuffer.data(fg_geom_1010_off + 914 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxxx = cbuffer.data(fg_geom_1010_off + 915 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxxy = cbuffer.data(fg_geom_1010_off + 916 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxxz = cbuffer.data(fg_geom_1010_off + 917 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxyy = cbuffer.data(fg_geom_1010_off + 918 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxyz = cbuffer.data(fg_geom_1010_off + 919 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxzz = cbuffer.data(fg_geom_1010_off + 920 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xyyy = cbuffer.data(fg_geom_1010_off + 921 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xyyz = cbuffer.data(fg_geom_1010_off + 922 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xyzz = cbuffer.data(fg_geom_1010_off + 923 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xzzz = cbuffer.data(fg_geom_1010_off + 924 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yyyy = cbuffer.data(fg_geom_1010_off + 925 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yyyz = cbuffer.data(fg_geom_1010_off + 926 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yyzz = cbuffer.data(fg_geom_1010_off + 927 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yzzz = cbuffer.data(fg_geom_1010_off + 928 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_zzzz = cbuffer.data(fg_geom_1010_off + 929 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxxx = cbuffer.data(fg_geom_1010_off + 930 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxxy = cbuffer.data(fg_geom_1010_off + 931 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxxz = cbuffer.data(fg_geom_1010_off + 932 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxyy = cbuffer.data(fg_geom_1010_off + 933 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxyz = cbuffer.data(fg_geom_1010_off + 934 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxzz = cbuffer.data(fg_geom_1010_off + 935 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xyyy = cbuffer.data(fg_geom_1010_off + 936 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xyyz = cbuffer.data(fg_geom_1010_off + 937 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xyzz = cbuffer.data(fg_geom_1010_off + 938 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xzzz = cbuffer.data(fg_geom_1010_off + 939 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yyyy = cbuffer.data(fg_geom_1010_off + 940 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yyyz = cbuffer.data(fg_geom_1010_off + 941 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yyzz = cbuffer.data(fg_geom_1010_off + 942 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yzzz = cbuffer.data(fg_geom_1010_off + 943 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_zzzz = cbuffer.data(fg_geom_1010_off + 944 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxxx = cbuffer.data(fg_geom_1010_off + 945 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxxy = cbuffer.data(fg_geom_1010_off + 946 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxxz = cbuffer.data(fg_geom_1010_off + 947 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxyy = cbuffer.data(fg_geom_1010_off + 948 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxyz = cbuffer.data(fg_geom_1010_off + 949 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxzz = cbuffer.data(fg_geom_1010_off + 950 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xyyy = cbuffer.data(fg_geom_1010_off + 951 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xyyz = cbuffer.data(fg_geom_1010_off + 952 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xyzz = cbuffer.data(fg_geom_1010_off + 953 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xzzz = cbuffer.data(fg_geom_1010_off + 954 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yyyy = cbuffer.data(fg_geom_1010_off + 955 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yyyz = cbuffer.data(fg_geom_1010_off + 956 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yyzz = cbuffer.data(fg_geom_1010_off + 957 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yzzz = cbuffer.data(fg_geom_1010_off + 958 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_zzzz = cbuffer.data(fg_geom_1010_off + 959 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxxx = cbuffer.data(fg_geom_1010_off + 960 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxxy = cbuffer.data(fg_geom_1010_off + 961 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxxz = cbuffer.data(fg_geom_1010_off + 962 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxyy = cbuffer.data(fg_geom_1010_off + 963 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxyz = cbuffer.data(fg_geom_1010_off + 964 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxzz = cbuffer.data(fg_geom_1010_off + 965 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xyyy = cbuffer.data(fg_geom_1010_off + 966 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xyyz = cbuffer.data(fg_geom_1010_off + 967 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xyzz = cbuffer.data(fg_geom_1010_off + 968 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xzzz = cbuffer.data(fg_geom_1010_off + 969 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yyyy = cbuffer.data(fg_geom_1010_off + 970 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yyyz = cbuffer.data(fg_geom_1010_off + 971 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yyzz = cbuffer.data(fg_geom_1010_off + 972 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yzzz = cbuffer.data(fg_geom_1010_off + 973 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_zzzz = cbuffer.data(fg_geom_1010_off + 974 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxxx = cbuffer.data(fg_geom_1010_off + 975 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxxy = cbuffer.data(fg_geom_1010_off + 976 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxxz = cbuffer.data(fg_geom_1010_off + 977 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxyy = cbuffer.data(fg_geom_1010_off + 978 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxyz = cbuffer.data(fg_geom_1010_off + 979 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxzz = cbuffer.data(fg_geom_1010_off + 980 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xyyy = cbuffer.data(fg_geom_1010_off + 981 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xyyz = cbuffer.data(fg_geom_1010_off + 982 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xyzz = cbuffer.data(fg_geom_1010_off + 983 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xzzz = cbuffer.data(fg_geom_1010_off + 984 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yyyy = cbuffer.data(fg_geom_1010_off + 985 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yyyz = cbuffer.data(fg_geom_1010_off + 986 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yyzz = cbuffer.data(fg_geom_1010_off + 987 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yzzz = cbuffer.data(fg_geom_1010_off + 988 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_zzzz = cbuffer.data(fg_geom_1010_off + 989 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxxx = cbuffer.data(fg_geom_1010_off + 990 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxxy = cbuffer.data(fg_geom_1010_off + 991 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxxz = cbuffer.data(fg_geom_1010_off + 992 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxyy = cbuffer.data(fg_geom_1010_off + 993 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxyz = cbuffer.data(fg_geom_1010_off + 994 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxzz = cbuffer.data(fg_geom_1010_off + 995 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xyyy = cbuffer.data(fg_geom_1010_off + 996 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xyyz = cbuffer.data(fg_geom_1010_off + 997 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xyzz = cbuffer.data(fg_geom_1010_off + 998 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xzzz = cbuffer.data(fg_geom_1010_off + 999 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yyyy = cbuffer.data(fg_geom_1010_off + 1000 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yyyz = cbuffer.data(fg_geom_1010_off + 1001 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yyzz = cbuffer.data(fg_geom_1010_off + 1002 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yzzz = cbuffer.data(fg_geom_1010_off + 1003 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_zzzz = cbuffer.data(fg_geom_1010_off + 1004 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxxx = cbuffer.data(fg_geom_1010_off + 1005 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxxy = cbuffer.data(fg_geom_1010_off + 1006 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxxz = cbuffer.data(fg_geom_1010_off + 1007 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxyy = cbuffer.data(fg_geom_1010_off + 1008 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxyz = cbuffer.data(fg_geom_1010_off + 1009 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxzz = cbuffer.data(fg_geom_1010_off + 1010 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xyyy = cbuffer.data(fg_geom_1010_off + 1011 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xyyz = cbuffer.data(fg_geom_1010_off + 1012 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xyzz = cbuffer.data(fg_geom_1010_off + 1013 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xzzz = cbuffer.data(fg_geom_1010_off + 1014 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yyyy = cbuffer.data(fg_geom_1010_off + 1015 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yyyz = cbuffer.data(fg_geom_1010_off + 1016 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yyzz = cbuffer.data(fg_geom_1010_off + 1017 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yzzz = cbuffer.data(fg_geom_1010_off + 1018 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_zzzz = cbuffer.data(fg_geom_1010_off + 1019 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxxx = cbuffer.data(fg_geom_1010_off + 1020 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxxy = cbuffer.data(fg_geom_1010_off + 1021 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxxz = cbuffer.data(fg_geom_1010_off + 1022 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxyy = cbuffer.data(fg_geom_1010_off + 1023 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxyz = cbuffer.data(fg_geom_1010_off + 1024 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxzz = cbuffer.data(fg_geom_1010_off + 1025 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xyyy = cbuffer.data(fg_geom_1010_off + 1026 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xyyz = cbuffer.data(fg_geom_1010_off + 1027 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xyzz = cbuffer.data(fg_geom_1010_off + 1028 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xzzz = cbuffer.data(fg_geom_1010_off + 1029 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yyyy = cbuffer.data(fg_geom_1010_off + 1030 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yyyz = cbuffer.data(fg_geom_1010_off + 1031 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yyzz = cbuffer.data(fg_geom_1010_off + 1032 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yzzz = cbuffer.data(fg_geom_1010_off + 1033 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_zzzz = cbuffer.data(fg_geom_1010_off + 1034 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxxx = cbuffer.data(fg_geom_1010_off + 1035 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxxy = cbuffer.data(fg_geom_1010_off + 1036 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxxz = cbuffer.data(fg_geom_1010_off + 1037 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxyy = cbuffer.data(fg_geom_1010_off + 1038 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxyz = cbuffer.data(fg_geom_1010_off + 1039 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxzz = cbuffer.data(fg_geom_1010_off + 1040 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xyyy = cbuffer.data(fg_geom_1010_off + 1041 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xyyz = cbuffer.data(fg_geom_1010_off + 1042 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xyzz = cbuffer.data(fg_geom_1010_off + 1043 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xzzz = cbuffer.data(fg_geom_1010_off + 1044 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yyyy = cbuffer.data(fg_geom_1010_off + 1045 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yyyz = cbuffer.data(fg_geom_1010_off + 1046 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yyzz = cbuffer.data(fg_geom_1010_off + 1047 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yzzz = cbuffer.data(fg_geom_1010_off + 1048 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_zzzz = cbuffer.data(fg_geom_1010_off + 1049 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxxx = cbuffer.data(fg_geom_1010_off + 1050 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxxy = cbuffer.data(fg_geom_1010_off + 1051 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxxz = cbuffer.data(fg_geom_1010_off + 1052 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxyy = cbuffer.data(fg_geom_1010_off + 1053 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxyz = cbuffer.data(fg_geom_1010_off + 1054 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxzz = cbuffer.data(fg_geom_1010_off + 1055 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xyyy = cbuffer.data(fg_geom_1010_off + 1056 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xyyz = cbuffer.data(fg_geom_1010_off + 1057 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xyzz = cbuffer.data(fg_geom_1010_off + 1058 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xzzz = cbuffer.data(fg_geom_1010_off + 1059 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yyyy = cbuffer.data(fg_geom_1010_off + 1060 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yyyz = cbuffer.data(fg_geom_1010_off + 1061 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yyzz = cbuffer.data(fg_geom_1010_off + 1062 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yzzz = cbuffer.data(fg_geom_1010_off + 1063 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_zzzz = cbuffer.data(fg_geom_1010_off + 1064 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxxx = cbuffer.data(fg_geom_1010_off + 1065 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxxy = cbuffer.data(fg_geom_1010_off + 1066 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxxz = cbuffer.data(fg_geom_1010_off + 1067 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxyy = cbuffer.data(fg_geom_1010_off + 1068 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxyz = cbuffer.data(fg_geom_1010_off + 1069 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxzz = cbuffer.data(fg_geom_1010_off + 1070 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xyyy = cbuffer.data(fg_geom_1010_off + 1071 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xyyz = cbuffer.data(fg_geom_1010_off + 1072 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xyzz = cbuffer.data(fg_geom_1010_off + 1073 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xzzz = cbuffer.data(fg_geom_1010_off + 1074 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yyyy = cbuffer.data(fg_geom_1010_off + 1075 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yyyz = cbuffer.data(fg_geom_1010_off + 1076 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yyzz = cbuffer.data(fg_geom_1010_off + 1077 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yzzz = cbuffer.data(fg_geom_1010_off + 1078 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_zzzz = cbuffer.data(fg_geom_1010_off + 1079 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxxx = cbuffer.data(fg_geom_1010_off + 1080 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxxy = cbuffer.data(fg_geom_1010_off + 1081 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxxz = cbuffer.data(fg_geom_1010_off + 1082 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxyy = cbuffer.data(fg_geom_1010_off + 1083 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxyz = cbuffer.data(fg_geom_1010_off + 1084 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxzz = cbuffer.data(fg_geom_1010_off + 1085 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xyyy = cbuffer.data(fg_geom_1010_off + 1086 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xyyz = cbuffer.data(fg_geom_1010_off + 1087 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xyzz = cbuffer.data(fg_geom_1010_off + 1088 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xzzz = cbuffer.data(fg_geom_1010_off + 1089 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yyyy = cbuffer.data(fg_geom_1010_off + 1090 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yyyz = cbuffer.data(fg_geom_1010_off + 1091 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yyzz = cbuffer.data(fg_geom_1010_off + 1092 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yzzz = cbuffer.data(fg_geom_1010_off + 1093 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_zzzz = cbuffer.data(fg_geom_1010_off + 1094 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxxx = cbuffer.data(fg_geom_1010_off + 1095 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxxy = cbuffer.data(fg_geom_1010_off + 1096 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxxz = cbuffer.data(fg_geom_1010_off + 1097 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxyy = cbuffer.data(fg_geom_1010_off + 1098 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxyz = cbuffer.data(fg_geom_1010_off + 1099 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxzz = cbuffer.data(fg_geom_1010_off + 1100 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xyyy = cbuffer.data(fg_geom_1010_off + 1101 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xyyz = cbuffer.data(fg_geom_1010_off + 1102 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xyzz = cbuffer.data(fg_geom_1010_off + 1103 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xzzz = cbuffer.data(fg_geom_1010_off + 1104 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yyyy = cbuffer.data(fg_geom_1010_off + 1105 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yyyz = cbuffer.data(fg_geom_1010_off + 1106 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yyzz = cbuffer.data(fg_geom_1010_off + 1107 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yzzz = cbuffer.data(fg_geom_1010_off + 1108 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_zzzz = cbuffer.data(fg_geom_1010_off + 1109 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxxx = cbuffer.data(fg_geom_1010_off + 1110 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxxy = cbuffer.data(fg_geom_1010_off + 1111 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxxz = cbuffer.data(fg_geom_1010_off + 1112 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxyy = cbuffer.data(fg_geom_1010_off + 1113 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxyz = cbuffer.data(fg_geom_1010_off + 1114 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxzz = cbuffer.data(fg_geom_1010_off + 1115 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xyyy = cbuffer.data(fg_geom_1010_off + 1116 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xyyz = cbuffer.data(fg_geom_1010_off + 1117 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xyzz = cbuffer.data(fg_geom_1010_off + 1118 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xzzz = cbuffer.data(fg_geom_1010_off + 1119 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yyyy = cbuffer.data(fg_geom_1010_off + 1120 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yyyz = cbuffer.data(fg_geom_1010_off + 1121 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yyzz = cbuffer.data(fg_geom_1010_off + 1122 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yzzz = cbuffer.data(fg_geom_1010_off + 1123 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_zzzz = cbuffer.data(fg_geom_1010_off + 1124 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxxx = cbuffer.data(fg_geom_1010_off + 1125 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxxy = cbuffer.data(fg_geom_1010_off + 1126 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxxz = cbuffer.data(fg_geom_1010_off + 1127 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxyy = cbuffer.data(fg_geom_1010_off + 1128 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxyz = cbuffer.data(fg_geom_1010_off + 1129 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxzz = cbuffer.data(fg_geom_1010_off + 1130 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xyyy = cbuffer.data(fg_geom_1010_off + 1131 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xyyz = cbuffer.data(fg_geom_1010_off + 1132 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xyzz = cbuffer.data(fg_geom_1010_off + 1133 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xzzz = cbuffer.data(fg_geom_1010_off + 1134 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yyyy = cbuffer.data(fg_geom_1010_off + 1135 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yyyz = cbuffer.data(fg_geom_1010_off + 1136 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yyzz = cbuffer.data(fg_geom_1010_off + 1137 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yzzz = cbuffer.data(fg_geom_1010_off + 1138 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_zzzz = cbuffer.data(fg_geom_1010_off + 1139 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxxx = cbuffer.data(fg_geom_1010_off + 1140 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxxy = cbuffer.data(fg_geom_1010_off + 1141 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxxz = cbuffer.data(fg_geom_1010_off + 1142 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxyy = cbuffer.data(fg_geom_1010_off + 1143 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxyz = cbuffer.data(fg_geom_1010_off + 1144 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxzz = cbuffer.data(fg_geom_1010_off + 1145 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xyyy = cbuffer.data(fg_geom_1010_off + 1146 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xyyz = cbuffer.data(fg_geom_1010_off + 1147 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xyzz = cbuffer.data(fg_geom_1010_off + 1148 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xzzz = cbuffer.data(fg_geom_1010_off + 1149 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yyyy = cbuffer.data(fg_geom_1010_off + 1150 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yyyz = cbuffer.data(fg_geom_1010_off + 1151 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yyzz = cbuffer.data(fg_geom_1010_off + 1152 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yzzz = cbuffer.data(fg_geom_1010_off + 1153 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_zzzz = cbuffer.data(fg_geom_1010_off + 1154 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxxx = cbuffer.data(fg_geom_1010_off + 1155 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxxy = cbuffer.data(fg_geom_1010_off + 1156 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxxz = cbuffer.data(fg_geom_1010_off + 1157 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxyy = cbuffer.data(fg_geom_1010_off + 1158 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxyz = cbuffer.data(fg_geom_1010_off + 1159 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxzz = cbuffer.data(fg_geom_1010_off + 1160 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xyyy = cbuffer.data(fg_geom_1010_off + 1161 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xyyz = cbuffer.data(fg_geom_1010_off + 1162 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xyzz = cbuffer.data(fg_geom_1010_off + 1163 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xzzz = cbuffer.data(fg_geom_1010_off + 1164 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yyyy = cbuffer.data(fg_geom_1010_off + 1165 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yyyz = cbuffer.data(fg_geom_1010_off + 1166 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yyzz = cbuffer.data(fg_geom_1010_off + 1167 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yzzz = cbuffer.data(fg_geom_1010_off + 1168 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_zzzz = cbuffer.data(fg_geom_1010_off + 1169 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxxx = cbuffer.data(fg_geom_1010_off + 1170 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxxy = cbuffer.data(fg_geom_1010_off + 1171 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxxz = cbuffer.data(fg_geom_1010_off + 1172 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxyy = cbuffer.data(fg_geom_1010_off + 1173 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxyz = cbuffer.data(fg_geom_1010_off + 1174 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxzz = cbuffer.data(fg_geom_1010_off + 1175 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xyyy = cbuffer.data(fg_geom_1010_off + 1176 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xyyz = cbuffer.data(fg_geom_1010_off + 1177 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xyzz = cbuffer.data(fg_geom_1010_off + 1178 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xzzz = cbuffer.data(fg_geom_1010_off + 1179 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yyyy = cbuffer.data(fg_geom_1010_off + 1180 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yyyz = cbuffer.data(fg_geom_1010_off + 1181 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yyzz = cbuffer.data(fg_geom_1010_off + 1182 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yzzz = cbuffer.data(fg_geom_1010_off + 1183 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_zzzz = cbuffer.data(fg_geom_1010_off + 1184 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxxx = cbuffer.data(fg_geom_1010_off + 1185 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxxy = cbuffer.data(fg_geom_1010_off + 1186 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxxz = cbuffer.data(fg_geom_1010_off + 1187 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxyy = cbuffer.data(fg_geom_1010_off + 1188 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxyz = cbuffer.data(fg_geom_1010_off + 1189 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxzz = cbuffer.data(fg_geom_1010_off + 1190 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xyyy = cbuffer.data(fg_geom_1010_off + 1191 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xyyz = cbuffer.data(fg_geom_1010_off + 1192 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xyzz = cbuffer.data(fg_geom_1010_off + 1193 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xzzz = cbuffer.data(fg_geom_1010_off + 1194 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yyyy = cbuffer.data(fg_geom_1010_off + 1195 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yyyz = cbuffer.data(fg_geom_1010_off + 1196 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yyzz = cbuffer.data(fg_geom_1010_off + 1197 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yzzz = cbuffer.data(fg_geom_1010_off + 1198 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_zzzz = cbuffer.data(fg_geom_1010_off + 1199 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxxx = cbuffer.data(fg_geom_1010_off + 1200 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxxy = cbuffer.data(fg_geom_1010_off + 1201 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxxz = cbuffer.data(fg_geom_1010_off + 1202 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxyy = cbuffer.data(fg_geom_1010_off + 1203 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxyz = cbuffer.data(fg_geom_1010_off + 1204 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxzz = cbuffer.data(fg_geom_1010_off + 1205 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xyyy = cbuffer.data(fg_geom_1010_off + 1206 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xyyz = cbuffer.data(fg_geom_1010_off + 1207 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xyzz = cbuffer.data(fg_geom_1010_off + 1208 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xzzz = cbuffer.data(fg_geom_1010_off + 1209 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yyyy = cbuffer.data(fg_geom_1010_off + 1210 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yyyz = cbuffer.data(fg_geom_1010_off + 1211 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yyzz = cbuffer.data(fg_geom_1010_off + 1212 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yzzz = cbuffer.data(fg_geom_1010_off + 1213 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_zzzz = cbuffer.data(fg_geom_1010_off + 1214 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxxx = cbuffer.data(fg_geom_1010_off + 1215 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxxy = cbuffer.data(fg_geom_1010_off + 1216 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxxz = cbuffer.data(fg_geom_1010_off + 1217 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxyy = cbuffer.data(fg_geom_1010_off + 1218 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxyz = cbuffer.data(fg_geom_1010_off + 1219 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxzz = cbuffer.data(fg_geom_1010_off + 1220 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xyyy = cbuffer.data(fg_geom_1010_off + 1221 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xyyz = cbuffer.data(fg_geom_1010_off + 1222 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xyzz = cbuffer.data(fg_geom_1010_off + 1223 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xzzz = cbuffer.data(fg_geom_1010_off + 1224 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yyyy = cbuffer.data(fg_geom_1010_off + 1225 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yyyz = cbuffer.data(fg_geom_1010_off + 1226 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yyzz = cbuffer.data(fg_geom_1010_off + 1227 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yzzz = cbuffer.data(fg_geom_1010_off + 1228 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_zzzz = cbuffer.data(fg_geom_1010_off + 1229 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxxx = cbuffer.data(fg_geom_1010_off + 1230 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxxy = cbuffer.data(fg_geom_1010_off + 1231 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxxz = cbuffer.data(fg_geom_1010_off + 1232 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxyy = cbuffer.data(fg_geom_1010_off + 1233 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxyz = cbuffer.data(fg_geom_1010_off + 1234 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxzz = cbuffer.data(fg_geom_1010_off + 1235 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xyyy = cbuffer.data(fg_geom_1010_off + 1236 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xyyz = cbuffer.data(fg_geom_1010_off + 1237 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xyzz = cbuffer.data(fg_geom_1010_off + 1238 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xzzz = cbuffer.data(fg_geom_1010_off + 1239 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yyyy = cbuffer.data(fg_geom_1010_off + 1240 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yyyz = cbuffer.data(fg_geom_1010_off + 1241 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yyzz = cbuffer.data(fg_geom_1010_off + 1242 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yzzz = cbuffer.data(fg_geom_1010_off + 1243 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_zzzz = cbuffer.data(fg_geom_1010_off + 1244 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxxx = cbuffer.data(fg_geom_1010_off + 1245 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxxy = cbuffer.data(fg_geom_1010_off + 1246 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxxz = cbuffer.data(fg_geom_1010_off + 1247 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxyy = cbuffer.data(fg_geom_1010_off + 1248 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxyz = cbuffer.data(fg_geom_1010_off + 1249 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxzz = cbuffer.data(fg_geom_1010_off + 1250 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xyyy = cbuffer.data(fg_geom_1010_off + 1251 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xyyz = cbuffer.data(fg_geom_1010_off + 1252 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xyzz = cbuffer.data(fg_geom_1010_off + 1253 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xzzz = cbuffer.data(fg_geom_1010_off + 1254 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yyyy = cbuffer.data(fg_geom_1010_off + 1255 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yyyz = cbuffer.data(fg_geom_1010_off + 1256 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yyzz = cbuffer.data(fg_geom_1010_off + 1257 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yzzz = cbuffer.data(fg_geom_1010_off + 1258 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_zzzz = cbuffer.data(fg_geom_1010_off + 1259 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxxx = cbuffer.data(fg_geom_1010_off + 1260 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxxy = cbuffer.data(fg_geom_1010_off + 1261 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxxz = cbuffer.data(fg_geom_1010_off + 1262 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxyy = cbuffer.data(fg_geom_1010_off + 1263 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxyz = cbuffer.data(fg_geom_1010_off + 1264 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxzz = cbuffer.data(fg_geom_1010_off + 1265 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xyyy = cbuffer.data(fg_geom_1010_off + 1266 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xyyz = cbuffer.data(fg_geom_1010_off + 1267 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xyzz = cbuffer.data(fg_geom_1010_off + 1268 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xzzz = cbuffer.data(fg_geom_1010_off + 1269 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yyyy = cbuffer.data(fg_geom_1010_off + 1270 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yyyz = cbuffer.data(fg_geom_1010_off + 1271 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yyzz = cbuffer.data(fg_geom_1010_off + 1272 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yzzz = cbuffer.data(fg_geom_1010_off + 1273 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_zzzz = cbuffer.data(fg_geom_1010_off + 1274 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxxx = cbuffer.data(fg_geom_1010_off + 1275 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxxy = cbuffer.data(fg_geom_1010_off + 1276 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxxz = cbuffer.data(fg_geom_1010_off + 1277 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxyy = cbuffer.data(fg_geom_1010_off + 1278 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxyz = cbuffer.data(fg_geom_1010_off + 1279 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxzz = cbuffer.data(fg_geom_1010_off + 1280 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xyyy = cbuffer.data(fg_geom_1010_off + 1281 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xyyz = cbuffer.data(fg_geom_1010_off + 1282 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xyzz = cbuffer.data(fg_geom_1010_off + 1283 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xzzz = cbuffer.data(fg_geom_1010_off + 1284 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yyyy = cbuffer.data(fg_geom_1010_off + 1285 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yyyz = cbuffer.data(fg_geom_1010_off + 1286 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yyzz = cbuffer.data(fg_geom_1010_off + 1287 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yzzz = cbuffer.data(fg_geom_1010_off + 1288 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_zzzz = cbuffer.data(fg_geom_1010_off + 1289 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxxx = cbuffer.data(fg_geom_1010_off + 1290 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxxy = cbuffer.data(fg_geom_1010_off + 1291 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxxz = cbuffer.data(fg_geom_1010_off + 1292 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxyy = cbuffer.data(fg_geom_1010_off + 1293 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxyz = cbuffer.data(fg_geom_1010_off + 1294 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxzz = cbuffer.data(fg_geom_1010_off + 1295 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xyyy = cbuffer.data(fg_geom_1010_off + 1296 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xyyz = cbuffer.data(fg_geom_1010_off + 1297 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xyzz = cbuffer.data(fg_geom_1010_off + 1298 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xzzz = cbuffer.data(fg_geom_1010_off + 1299 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yyyy = cbuffer.data(fg_geom_1010_off + 1300 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yyyz = cbuffer.data(fg_geom_1010_off + 1301 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yyzz = cbuffer.data(fg_geom_1010_off + 1302 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yzzz = cbuffer.data(fg_geom_1010_off + 1303 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_zzzz = cbuffer.data(fg_geom_1010_off + 1304 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxxx = cbuffer.data(fg_geom_1010_off + 1305 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxxy = cbuffer.data(fg_geom_1010_off + 1306 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxxz = cbuffer.data(fg_geom_1010_off + 1307 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxyy = cbuffer.data(fg_geom_1010_off + 1308 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxyz = cbuffer.data(fg_geom_1010_off + 1309 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxzz = cbuffer.data(fg_geom_1010_off + 1310 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xyyy = cbuffer.data(fg_geom_1010_off + 1311 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xyyz = cbuffer.data(fg_geom_1010_off + 1312 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xyzz = cbuffer.data(fg_geom_1010_off + 1313 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xzzz = cbuffer.data(fg_geom_1010_off + 1314 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yyyy = cbuffer.data(fg_geom_1010_off + 1315 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yyyz = cbuffer.data(fg_geom_1010_off + 1316 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yyzz = cbuffer.data(fg_geom_1010_off + 1317 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yzzz = cbuffer.data(fg_geom_1010_off + 1318 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_zzzz = cbuffer.data(fg_geom_1010_off + 1319 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxxx = cbuffer.data(fg_geom_1010_off + 1320 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxxy = cbuffer.data(fg_geom_1010_off + 1321 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxxz = cbuffer.data(fg_geom_1010_off + 1322 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxyy = cbuffer.data(fg_geom_1010_off + 1323 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxyz = cbuffer.data(fg_geom_1010_off + 1324 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxzz = cbuffer.data(fg_geom_1010_off + 1325 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xyyy = cbuffer.data(fg_geom_1010_off + 1326 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xyyz = cbuffer.data(fg_geom_1010_off + 1327 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xyzz = cbuffer.data(fg_geom_1010_off + 1328 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xzzz = cbuffer.data(fg_geom_1010_off + 1329 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yyyy = cbuffer.data(fg_geom_1010_off + 1330 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yyyz = cbuffer.data(fg_geom_1010_off + 1331 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yyzz = cbuffer.data(fg_geom_1010_off + 1332 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yzzz = cbuffer.data(fg_geom_1010_off + 1333 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_zzzz = cbuffer.data(fg_geom_1010_off + 1334 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxxx = cbuffer.data(fg_geom_1010_off + 1335 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxxy = cbuffer.data(fg_geom_1010_off + 1336 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxxz = cbuffer.data(fg_geom_1010_off + 1337 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxyy = cbuffer.data(fg_geom_1010_off + 1338 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxyz = cbuffer.data(fg_geom_1010_off + 1339 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxzz = cbuffer.data(fg_geom_1010_off + 1340 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xyyy = cbuffer.data(fg_geom_1010_off + 1341 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xyyz = cbuffer.data(fg_geom_1010_off + 1342 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xyzz = cbuffer.data(fg_geom_1010_off + 1343 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xzzz = cbuffer.data(fg_geom_1010_off + 1344 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yyyy = cbuffer.data(fg_geom_1010_off + 1345 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yyyz = cbuffer.data(fg_geom_1010_off + 1346 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yyzz = cbuffer.data(fg_geom_1010_off + 1347 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yzzz = cbuffer.data(fg_geom_1010_off + 1348 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_zzzz = cbuffer.data(fg_geom_1010_off + 1349 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FHSS

            const auto fh_geom_1010_off = idx_geom_1010_fhxx + i * dcomps + j;

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

            /// set up bra offset for contr_buffer_ggxx

            const auto gg_geom_1010_off = idx_geom_1010_ggxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxx_xxxx = cbuffer.data(gg_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_xxxy = cbuffer.data(gg_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_xxxz = cbuffer.data(gg_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_xxyy = cbuffer.data(gg_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_xxyz = cbuffer.data(gg_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_xxzz = cbuffer.data(gg_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_xyyy = cbuffer.data(gg_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_xyyz = cbuffer.data(gg_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_xyzz = cbuffer.data(gg_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_xzzz = cbuffer.data(gg_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_yyyy = cbuffer.data(gg_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_yyyz = cbuffer.data(gg_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_yyzz = cbuffer.data(gg_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_yzzz = cbuffer.data(gg_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_zzzz = cbuffer.data(gg_geom_1010_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_xxx_xxxx, g_0_0_x_0_xxx_xxxy, g_0_0_x_0_xxx_xxxz, g_0_0_x_0_xxx_xxyy, g_0_0_x_0_xxx_xxyz, g_0_0_x_0_xxx_xxzz, g_0_0_x_0_xxx_xyyy, g_0_0_x_0_xxx_xyyz, g_0_0_x_0_xxx_xyzz, g_0_0_x_0_xxx_xzzz, g_0_0_x_0_xxx_yyyy, g_0_0_x_0_xxx_yyyz, g_0_0_x_0_xxx_yyzz, g_0_0_x_0_xxx_yzzz, g_0_0_x_0_xxx_zzzz, g_x_0_x_0_xxx_xxxx, g_x_0_x_0_xxx_xxxxx, g_x_0_x_0_xxx_xxxxy, g_x_0_x_0_xxx_xxxxz, g_x_0_x_0_xxx_xxxy, g_x_0_x_0_xxx_xxxyy, g_x_0_x_0_xxx_xxxyz, g_x_0_x_0_xxx_xxxz, g_x_0_x_0_xxx_xxxzz, g_x_0_x_0_xxx_xxyy, g_x_0_x_0_xxx_xxyyy, g_x_0_x_0_xxx_xxyyz, g_x_0_x_0_xxx_xxyz, g_x_0_x_0_xxx_xxyzz, g_x_0_x_0_xxx_xxzz, g_x_0_x_0_xxx_xxzzz, g_x_0_x_0_xxx_xyyy, g_x_0_x_0_xxx_xyyyy, g_x_0_x_0_xxx_xyyyz, g_x_0_x_0_xxx_xyyz, g_x_0_x_0_xxx_xyyzz, g_x_0_x_0_xxx_xyzz, g_x_0_x_0_xxx_xyzzz, g_x_0_x_0_xxx_xzzz, g_x_0_x_0_xxx_xzzzz, g_x_0_x_0_xxx_yyyy, g_x_0_x_0_xxx_yyyz, g_x_0_x_0_xxx_yyzz, g_x_0_x_0_xxx_yzzz, g_x_0_x_0_xxx_zzzz, g_x_0_x_0_xxxx_xxxx, g_x_0_x_0_xxxx_xxxy, g_x_0_x_0_xxxx_xxxz, g_x_0_x_0_xxxx_xxyy, g_x_0_x_0_xxxx_xxyz, g_x_0_x_0_xxxx_xxzz, g_x_0_x_0_xxxx_xyyy, g_x_0_x_0_xxxx_xyyz, g_x_0_x_0_xxxx_xyzz, g_x_0_x_0_xxxx_xzzz, g_x_0_x_0_xxxx_yyyy, g_x_0_x_0_xxxx_yyyz, g_x_0_x_0_xxxx_yyzz, g_x_0_x_0_xxxx_yzzz, g_x_0_x_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxx_xxxx[k] = -g_0_0_x_0_xxx_xxxx[k] - g_x_0_x_0_xxx_xxxx[k] * ab_x + g_x_0_x_0_xxx_xxxxx[k];

                g_x_0_x_0_xxxx_xxxy[k] = -g_0_0_x_0_xxx_xxxy[k] - g_x_0_x_0_xxx_xxxy[k] * ab_x + g_x_0_x_0_xxx_xxxxy[k];

                g_x_0_x_0_xxxx_xxxz[k] = -g_0_0_x_0_xxx_xxxz[k] - g_x_0_x_0_xxx_xxxz[k] * ab_x + g_x_0_x_0_xxx_xxxxz[k];

                g_x_0_x_0_xxxx_xxyy[k] = -g_0_0_x_0_xxx_xxyy[k] - g_x_0_x_0_xxx_xxyy[k] * ab_x + g_x_0_x_0_xxx_xxxyy[k];

                g_x_0_x_0_xxxx_xxyz[k] = -g_0_0_x_0_xxx_xxyz[k] - g_x_0_x_0_xxx_xxyz[k] * ab_x + g_x_0_x_0_xxx_xxxyz[k];

                g_x_0_x_0_xxxx_xxzz[k] = -g_0_0_x_0_xxx_xxzz[k] - g_x_0_x_0_xxx_xxzz[k] * ab_x + g_x_0_x_0_xxx_xxxzz[k];

                g_x_0_x_0_xxxx_xyyy[k] = -g_0_0_x_0_xxx_xyyy[k] - g_x_0_x_0_xxx_xyyy[k] * ab_x + g_x_0_x_0_xxx_xxyyy[k];

                g_x_0_x_0_xxxx_xyyz[k] = -g_0_0_x_0_xxx_xyyz[k] - g_x_0_x_0_xxx_xyyz[k] * ab_x + g_x_0_x_0_xxx_xxyyz[k];

                g_x_0_x_0_xxxx_xyzz[k] = -g_0_0_x_0_xxx_xyzz[k] - g_x_0_x_0_xxx_xyzz[k] * ab_x + g_x_0_x_0_xxx_xxyzz[k];

                g_x_0_x_0_xxxx_xzzz[k] = -g_0_0_x_0_xxx_xzzz[k] - g_x_0_x_0_xxx_xzzz[k] * ab_x + g_x_0_x_0_xxx_xxzzz[k];

                g_x_0_x_0_xxxx_yyyy[k] = -g_0_0_x_0_xxx_yyyy[k] - g_x_0_x_0_xxx_yyyy[k] * ab_x + g_x_0_x_0_xxx_xyyyy[k];

                g_x_0_x_0_xxxx_yyyz[k] = -g_0_0_x_0_xxx_yyyz[k] - g_x_0_x_0_xxx_yyyz[k] * ab_x + g_x_0_x_0_xxx_xyyyz[k];

                g_x_0_x_0_xxxx_yyzz[k] = -g_0_0_x_0_xxx_yyzz[k] - g_x_0_x_0_xxx_yyzz[k] * ab_x + g_x_0_x_0_xxx_xyyzz[k];

                g_x_0_x_0_xxxx_yzzz[k] = -g_0_0_x_0_xxx_yzzz[k] - g_x_0_x_0_xxx_yzzz[k] * ab_x + g_x_0_x_0_xxx_xyzzz[k];

                g_x_0_x_0_xxxx_zzzz[k] = -g_0_0_x_0_xxx_zzzz[k] - g_x_0_x_0_xxx_zzzz[k] * ab_x + g_x_0_x_0_xxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxy_xxxx = cbuffer.data(gg_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_xxxy = cbuffer.data(gg_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_xxxz = cbuffer.data(gg_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_xxyy = cbuffer.data(gg_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_xxyz = cbuffer.data(gg_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_xxzz = cbuffer.data(gg_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_xyyy = cbuffer.data(gg_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_xyyz = cbuffer.data(gg_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_xyzz = cbuffer.data(gg_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_xzzz = cbuffer.data(gg_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_yyyy = cbuffer.data(gg_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_yyyz = cbuffer.data(gg_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_yyzz = cbuffer.data(gg_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_yzzz = cbuffer.data(gg_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_zzzz = cbuffer.data(gg_geom_1010_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxx_xxxx, g_x_0_x_0_xxx_xxxxy, g_x_0_x_0_xxx_xxxy, g_x_0_x_0_xxx_xxxyy, g_x_0_x_0_xxx_xxxyz, g_x_0_x_0_xxx_xxxz, g_x_0_x_0_xxx_xxyy, g_x_0_x_0_xxx_xxyyy, g_x_0_x_0_xxx_xxyyz, g_x_0_x_0_xxx_xxyz, g_x_0_x_0_xxx_xxyzz, g_x_0_x_0_xxx_xxzz, g_x_0_x_0_xxx_xyyy, g_x_0_x_0_xxx_xyyyy, g_x_0_x_0_xxx_xyyyz, g_x_0_x_0_xxx_xyyz, g_x_0_x_0_xxx_xyyzz, g_x_0_x_0_xxx_xyzz, g_x_0_x_0_xxx_xyzzz, g_x_0_x_0_xxx_xzzz, g_x_0_x_0_xxx_yyyy, g_x_0_x_0_xxx_yyyyy, g_x_0_x_0_xxx_yyyyz, g_x_0_x_0_xxx_yyyz, g_x_0_x_0_xxx_yyyzz, g_x_0_x_0_xxx_yyzz, g_x_0_x_0_xxx_yyzzz, g_x_0_x_0_xxx_yzzz, g_x_0_x_0_xxx_yzzzz, g_x_0_x_0_xxx_zzzz, g_x_0_x_0_xxxy_xxxx, g_x_0_x_0_xxxy_xxxy, g_x_0_x_0_xxxy_xxxz, g_x_0_x_0_xxxy_xxyy, g_x_0_x_0_xxxy_xxyz, g_x_0_x_0_xxxy_xxzz, g_x_0_x_0_xxxy_xyyy, g_x_0_x_0_xxxy_xyyz, g_x_0_x_0_xxxy_xyzz, g_x_0_x_0_xxxy_xzzz, g_x_0_x_0_xxxy_yyyy, g_x_0_x_0_xxxy_yyyz, g_x_0_x_0_xxxy_yyzz, g_x_0_x_0_xxxy_yzzz, g_x_0_x_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxy_xxxx[k] = -g_x_0_x_0_xxx_xxxx[k] * ab_y + g_x_0_x_0_xxx_xxxxy[k];

                g_x_0_x_0_xxxy_xxxy[k] = -g_x_0_x_0_xxx_xxxy[k] * ab_y + g_x_0_x_0_xxx_xxxyy[k];

                g_x_0_x_0_xxxy_xxxz[k] = -g_x_0_x_0_xxx_xxxz[k] * ab_y + g_x_0_x_0_xxx_xxxyz[k];

                g_x_0_x_0_xxxy_xxyy[k] = -g_x_0_x_0_xxx_xxyy[k] * ab_y + g_x_0_x_0_xxx_xxyyy[k];

                g_x_0_x_0_xxxy_xxyz[k] = -g_x_0_x_0_xxx_xxyz[k] * ab_y + g_x_0_x_0_xxx_xxyyz[k];

                g_x_0_x_0_xxxy_xxzz[k] = -g_x_0_x_0_xxx_xxzz[k] * ab_y + g_x_0_x_0_xxx_xxyzz[k];

                g_x_0_x_0_xxxy_xyyy[k] = -g_x_0_x_0_xxx_xyyy[k] * ab_y + g_x_0_x_0_xxx_xyyyy[k];

                g_x_0_x_0_xxxy_xyyz[k] = -g_x_0_x_0_xxx_xyyz[k] * ab_y + g_x_0_x_0_xxx_xyyyz[k];

                g_x_0_x_0_xxxy_xyzz[k] = -g_x_0_x_0_xxx_xyzz[k] * ab_y + g_x_0_x_0_xxx_xyyzz[k];

                g_x_0_x_0_xxxy_xzzz[k] = -g_x_0_x_0_xxx_xzzz[k] * ab_y + g_x_0_x_0_xxx_xyzzz[k];

                g_x_0_x_0_xxxy_yyyy[k] = -g_x_0_x_0_xxx_yyyy[k] * ab_y + g_x_0_x_0_xxx_yyyyy[k];

                g_x_0_x_0_xxxy_yyyz[k] = -g_x_0_x_0_xxx_yyyz[k] * ab_y + g_x_0_x_0_xxx_yyyyz[k];

                g_x_0_x_0_xxxy_yyzz[k] = -g_x_0_x_0_xxx_yyzz[k] * ab_y + g_x_0_x_0_xxx_yyyzz[k];

                g_x_0_x_0_xxxy_yzzz[k] = -g_x_0_x_0_xxx_yzzz[k] * ab_y + g_x_0_x_0_xxx_yyzzz[k];

                g_x_0_x_0_xxxy_zzzz[k] = -g_x_0_x_0_xxx_zzzz[k] * ab_y + g_x_0_x_0_xxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxz_xxxx = cbuffer.data(gg_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_xxxy = cbuffer.data(gg_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_xxxz = cbuffer.data(gg_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_xxyy = cbuffer.data(gg_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_xxyz = cbuffer.data(gg_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_xxzz = cbuffer.data(gg_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_xyyy = cbuffer.data(gg_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_xyyz = cbuffer.data(gg_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_xyzz = cbuffer.data(gg_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_xzzz = cbuffer.data(gg_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_yyyy = cbuffer.data(gg_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_yyyz = cbuffer.data(gg_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_yyzz = cbuffer.data(gg_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_yzzz = cbuffer.data(gg_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_zzzz = cbuffer.data(gg_geom_1010_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxx_xxxx, g_x_0_x_0_xxx_xxxxz, g_x_0_x_0_xxx_xxxy, g_x_0_x_0_xxx_xxxyz, g_x_0_x_0_xxx_xxxz, g_x_0_x_0_xxx_xxxzz, g_x_0_x_0_xxx_xxyy, g_x_0_x_0_xxx_xxyyz, g_x_0_x_0_xxx_xxyz, g_x_0_x_0_xxx_xxyzz, g_x_0_x_0_xxx_xxzz, g_x_0_x_0_xxx_xxzzz, g_x_0_x_0_xxx_xyyy, g_x_0_x_0_xxx_xyyyz, g_x_0_x_0_xxx_xyyz, g_x_0_x_0_xxx_xyyzz, g_x_0_x_0_xxx_xyzz, g_x_0_x_0_xxx_xyzzz, g_x_0_x_0_xxx_xzzz, g_x_0_x_0_xxx_xzzzz, g_x_0_x_0_xxx_yyyy, g_x_0_x_0_xxx_yyyyz, g_x_0_x_0_xxx_yyyz, g_x_0_x_0_xxx_yyyzz, g_x_0_x_0_xxx_yyzz, g_x_0_x_0_xxx_yyzzz, g_x_0_x_0_xxx_yzzz, g_x_0_x_0_xxx_yzzzz, g_x_0_x_0_xxx_zzzz, g_x_0_x_0_xxx_zzzzz, g_x_0_x_0_xxxz_xxxx, g_x_0_x_0_xxxz_xxxy, g_x_0_x_0_xxxz_xxxz, g_x_0_x_0_xxxz_xxyy, g_x_0_x_0_xxxz_xxyz, g_x_0_x_0_xxxz_xxzz, g_x_0_x_0_xxxz_xyyy, g_x_0_x_0_xxxz_xyyz, g_x_0_x_0_xxxz_xyzz, g_x_0_x_0_xxxz_xzzz, g_x_0_x_0_xxxz_yyyy, g_x_0_x_0_xxxz_yyyz, g_x_0_x_0_xxxz_yyzz, g_x_0_x_0_xxxz_yzzz, g_x_0_x_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxz_xxxx[k] = -g_x_0_x_0_xxx_xxxx[k] * ab_z + g_x_0_x_0_xxx_xxxxz[k];

                g_x_0_x_0_xxxz_xxxy[k] = -g_x_0_x_0_xxx_xxxy[k] * ab_z + g_x_0_x_0_xxx_xxxyz[k];

                g_x_0_x_0_xxxz_xxxz[k] = -g_x_0_x_0_xxx_xxxz[k] * ab_z + g_x_0_x_0_xxx_xxxzz[k];

                g_x_0_x_0_xxxz_xxyy[k] = -g_x_0_x_0_xxx_xxyy[k] * ab_z + g_x_0_x_0_xxx_xxyyz[k];

                g_x_0_x_0_xxxz_xxyz[k] = -g_x_0_x_0_xxx_xxyz[k] * ab_z + g_x_0_x_0_xxx_xxyzz[k];

                g_x_0_x_0_xxxz_xxzz[k] = -g_x_0_x_0_xxx_xxzz[k] * ab_z + g_x_0_x_0_xxx_xxzzz[k];

                g_x_0_x_0_xxxz_xyyy[k] = -g_x_0_x_0_xxx_xyyy[k] * ab_z + g_x_0_x_0_xxx_xyyyz[k];

                g_x_0_x_0_xxxz_xyyz[k] = -g_x_0_x_0_xxx_xyyz[k] * ab_z + g_x_0_x_0_xxx_xyyzz[k];

                g_x_0_x_0_xxxz_xyzz[k] = -g_x_0_x_0_xxx_xyzz[k] * ab_z + g_x_0_x_0_xxx_xyzzz[k];

                g_x_0_x_0_xxxz_xzzz[k] = -g_x_0_x_0_xxx_xzzz[k] * ab_z + g_x_0_x_0_xxx_xzzzz[k];

                g_x_0_x_0_xxxz_yyyy[k] = -g_x_0_x_0_xxx_yyyy[k] * ab_z + g_x_0_x_0_xxx_yyyyz[k];

                g_x_0_x_0_xxxz_yyyz[k] = -g_x_0_x_0_xxx_yyyz[k] * ab_z + g_x_0_x_0_xxx_yyyzz[k];

                g_x_0_x_0_xxxz_yyzz[k] = -g_x_0_x_0_xxx_yyzz[k] * ab_z + g_x_0_x_0_xxx_yyzzz[k];

                g_x_0_x_0_xxxz_yzzz[k] = -g_x_0_x_0_xxx_yzzz[k] * ab_z + g_x_0_x_0_xxx_yzzzz[k];

                g_x_0_x_0_xxxz_zzzz[k] = -g_x_0_x_0_xxx_zzzz[k] * ab_z + g_x_0_x_0_xxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyy_xxxx = cbuffer.data(gg_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_xxxy = cbuffer.data(gg_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_xxxz = cbuffer.data(gg_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_xxyy = cbuffer.data(gg_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_xxyz = cbuffer.data(gg_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_xxzz = cbuffer.data(gg_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_xyyy = cbuffer.data(gg_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_xyyz = cbuffer.data(gg_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_xyzz = cbuffer.data(gg_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_xzzz = cbuffer.data(gg_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_yyyy = cbuffer.data(gg_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_yyyz = cbuffer.data(gg_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_yyzz = cbuffer.data(gg_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_yzzz = cbuffer.data(gg_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_zzzz = cbuffer.data(gg_geom_1010_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxy_xxxx, g_x_0_x_0_xxy_xxxxy, g_x_0_x_0_xxy_xxxy, g_x_0_x_0_xxy_xxxyy, g_x_0_x_0_xxy_xxxyz, g_x_0_x_0_xxy_xxxz, g_x_0_x_0_xxy_xxyy, g_x_0_x_0_xxy_xxyyy, g_x_0_x_0_xxy_xxyyz, g_x_0_x_0_xxy_xxyz, g_x_0_x_0_xxy_xxyzz, g_x_0_x_0_xxy_xxzz, g_x_0_x_0_xxy_xyyy, g_x_0_x_0_xxy_xyyyy, g_x_0_x_0_xxy_xyyyz, g_x_0_x_0_xxy_xyyz, g_x_0_x_0_xxy_xyyzz, g_x_0_x_0_xxy_xyzz, g_x_0_x_0_xxy_xyzzz, g_x_0_x_0_xxy_xzzz, g_x_0_x_0_xxy_yyyy, g_x_0_x_0_xxy_yyyyy, g_x_0_x_0_xxy_yyyyz, g_x_0_x_0_xxy_yyyz, g_x_0_x_0_xxy_yyyzz, g_x_0_x_0_xxy_yyzz, g_x_0_x_0_xxy_yyzzz, g_x_0_x_0_xxy_yzzz, g_x_0_x_0_xxy_yzzzz, g_x_0_x_0_xxy_zzzz, g_x_0_x_0_xxyy_xxxx, g_x_0_x_0_xxyy_xxxy, g_x_0_x_0_xxyy_xxxz, g_x_0_x_0_xxyy_xxyy, g_x_0_x_0_xxyy_xxyz, g_x_0_x_0_xxyy_xxzz, g_x_0_x_0_xxyy_xyyy, g_x_0_x_0_xxyy_xyyz, g_x_0_x_0_xxyy_xyzz, g_x_0_x_0_xxyy_xzzz, g_x_0_x_0_xxyy_yyyy, g_x_0_x_0_xxyy_yyyz, g_x_0_x_0_xxyy_yyzz, g_x_0_x_0_xxyy_yzzz, g_x_0_x_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyy_xxxx[k] = -g_x_0_x_0_xxy_xxxx[k] * ab_y + g_x_0_x_0_xxy_xxxxy[k];

                g_x_0_x_0_xxyy_xxxy[k] = -g_x_0_x_0_xxy_xxxy[k] * ab_y + g_x_0_x_0_xxy_xxxyy[k];

                g_x_0_x_0_xxyy_xxxz[k] = -g_x_0_x_0_xxy_xxxz[k] * ab_y + g_x_0_x_0_xxy_xxxyz[k];

                g_x_0_x_0_xxyy_xxyy[k] = -g_x_0_x_0_xxy_xxyy[k] * ab_y + g_x_0_x_0_xxy_xxyyy[k];

                g_x_0_x_0_xxyy_xxyz[k] = -g_x_0_x_0_xxy_xxyz[k] * ab_y + g_x_0_x_0_xxy_xxyyz[k];

                g_x_0_x_0_xxyy_xxzz[k] = -g_x_0_x_0_xxy_xxzz[k] * ab_y + g_x_0_x_0_xxy_xxyzz[k];

                g_x_0_x_0_xxyy_xyyy[k] = -g_x_0_x_0_xxy_xyyy[k] * ab_y + g_x_0_x_0_xxy_xyyyy[k];

                g_x_0_x_0_xxyy_xyyz[k] = -g_x_0_x_0_xxy_xyyz[k] * ab_y + g_x_0_x_0_xxy_xyyyz[k];

                g_x_0_x_0_xxyy_xyzz[k] = -g_x_0_x_0_xxy_xyzz[k] * ab_y + g_x_0_x_0_xxy_xyyzz[k];

                g_x_0_x_0_xxyy_xzzz[k] = -g_x_0_x_0_xxy_xzzz[k] * ab_y + g_x_0_x_0_xxy_xyzzz[k];

                g_x_0_x_0_xxyy_yyyy[k] = -g_x_0_x_0_xxy_yyyy[k] * ab_y + g_x_0_x_0_xxy_yyyyy[k];

                g_x_0_x_0_xxyy_yyyz[k] = -g_x_0_x_0_xxy_yyyz[k] * ab_y + g_x_0_x_0_xxy_yyyyz[k];

                g_x_0_x_0_xxyy_yyzz[k] = -g_x_0_x_0_xxy_yyzz[k] * ab_y + g_x_0_x_0_xxy_yyyzz[k];

                g_x_0_x_0_xxyy_yzzz[k] = -g_x_0_x_0_xxy_yzzz[k] * ab_y + g_x_0_x_0_xxy_yyzzz[k];

                g_x_0_x_0_xxyy_zzzz[k] = -g_x_0_x_0_xxy_zzzz[k] * ab_y + g_x_0_x_0_xxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyz_xxxx = cbuffer.data(gg_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_xxxy = cbuffer.data(gg_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_xxxz = cbuffer.data(gg_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_xxyy = cbuffer.data(gg_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_xxyz = cbuffer.data(gg_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_xxzz = cbuffer.data(gg_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_xyyy = cbuffer.data(gg_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_xyyz = cbuffer.data(gg_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_xyzz = cbuffer.data(gg_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_xzzz = cbuffer.data(gg_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_yyyy = cbuffer.data(gg_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_yyyz = cbuffer.data(gg_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_yyzz = cbuffer.data(gg_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_yzzz = cbuffer.data(gg_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_zzzz = cbuffer.data(gg_geom_1010_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyz_xxxx, g_x_0_x_0_xxyz_xxxy, g_x_0_x_0_xxyz_xxxz, g_x_0_x_0_xxyz_xxyy, g_x_0_x_0_xxyz_xxyz, g_x_0_x_0_xxyz_xxzz, g_x_0_x_0_xxyz_xyyy, g_x_0_x_0_xxyz_xyyz, g_x_0_x_0_xxyz_xyzz, g_x_0_x_0_xxyz_xzzz, g_x_0_x_0_xxyz_yyyy, g_x_0_x_0_xxyz_yyyz, g_x_0_x_0_xxyz_yyzz, g_x_0_x_0_xxyz_yzzz, g_x_0_x_0_xxyz_zzzz, g_x_0_x_0_xxz_xxxx, g_x_0_x_0_xxz_xxxxy, g_x_0_x_0_xxz_xxxy, g_x_0_x_0_xxz_xxxyy, g_x_0_x_0_xxz_xxxyz, g_x_0_x_0_xxz_xxxz, g_x_0_x_0_xxz_xxyy, g_x_0_x_0_xxz_xxyyy, g_x_0_x_0_xxz_xxyyz, g_x_0_x_0_xxz_xxyz, g_x_0_x_0_xxz_xxyzz, g_x_0_x_0_xxz_xxzz, g_x_0_x_0_xxz_xyyy, g_x_0_x_0_xxz_xyyyy, g_x_0_x_0_xxz_xyyyz, g_x_0_x_0_xxz_xyyz, g_x_0_x_0_xxz_xyyzz, g_x_0_x_0_xxz_xyzz, g_x_0_x_0_xxz_xyzzz, g_x_0_x_0_xxz_xzzz, g_x_0_x_0_xxz_yyyy, g_x_0_x_0_xxz_yyyyy, g_x_0_x_0_xxz_yyyyz, g_x_0_x_0_xxz_yyyz, g_x_0_x_0_xxz_yyyzz, g_x_0_x_0_xxz_yyzz, g_x_0_x_0_xxz_yyzzz, g_x_0_x_0_xxz_yzzz, g_x_0_x_0_xxz_yzzzz, g_x_0_x_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyz_xxxx[k] = -g_x_0_x_0_xxz_xxxx[k] * ab_y + g_x_0_x_0_xxz_xxxxy[k];

                g_x_0_x_0_xxyz_xxxy[k] = -g_x_0_x_0_xxz_xxxy[k] * ab_y + g_x_0_x_0_xxz_xxxyy[k];

                g_x_0_x_0_xxyz_xxxz[k] = -g_x_0_x_0_xxz_xxxz[k] * ab_y + g_x_0_x_0_xxz_xxxyz[k];

                g_x_0_x_0_xxyz_xxyy[k] = -g_x_0_x_0_xxz_xxyy[k] * ab_y + g_x_0_x_0_xxz_xxyyy[k];

                g_x_0_x_0_xxyz_xxyz[k] = -g_x_0_x_0_xxz_xxyz[k] * ab_y + g_x_0_x_0_xxz_xxyyz[k];

                g_x_0_x_0_xxyz_xxzz[k] = -g_x_0_x_0_xxz_xxzz[k] * ab_y + g_x_0_x_0_xxz_xxyzz[k];

                g_x_0_x_0_xxyz_xyyy[k] = -g_x_0_x_0_xxz_xyyy[k] * ab_y + g_x_0_x_0_xxz_xyyyy[k];

                g_x_0_x_0_xxyz_xyyz[k] = -g_x_0_x_0_xxz_xyyz[k] * ab_y + g_x_0_x_0_xxz_xyyyz[k];

                g_x_0_x_0_xxyz_xyzz[k] = -g_x_0_x_0_xxz_xyzz[k] * ab_y + g_x_0_x_0_xxz_xyyzz[k];

                g_x_0_x_0_xxyz_xzzz[k] = -g_x_0_x_0_xxz_xzzz[k] * ab_y + g_x_0_x_0_xxz_xyzzz[k];

                g_x_0_x_0_xxyz_yyyy[k] = -g_x_0_x_0_xxz_yyyy[k] * ab_y + g_x_0_x_0_xxz_yyyyy[k];

                g_x_0_x_0_xxyz_yyyz[k] = -g_x_0_x_0_xxz_yyyz[k] * ab_y + g_x_0_x_0_xxz_yyyyz[k];

                g_x_0_x_0_xxyz_yyzz[k] = -g_x_0_x_0_xxz_yyzz[k] * ab_y + g_x_0_x_0_xxz_yyyzz[k];

                g_x_0_x_0_xxyz_yzzz[k] = -g_x_0_x_0_xxz_yzzz[k] * ab_y + g_x_0_x_0_xxz_yyzzz[k];

                g_x_0_x_0_xxyz_zzzz[k] = -g_x_0_x_0_xxz_zzzz[k] * ab_y + g_x_0_x_0_xxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxzz_xxxx = cbuffer.data(gg_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_xxxy = cbuffer.data(gg_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_xxxz = cbuffer.data(gg_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_xxyy = cbuffer.data(gg_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_xxyz = cbuffer.data(gg_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_xxzz = cbuffer.data(gg_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_xyyy = cbuffer.data(gg_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_xyyz = cbuffer.data(gg_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_xyzz = cbuffer.data(gg_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_xzzz = cbuffer.data(gg_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_yyyy = cbuffer.data(gg_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_yyyz = cbuffer.data(gg_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_yyzz = cbuffer.data(gg_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_yzzz = cbuffer.data(gg_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_zzzz = cbuffer.data(gg_geom_1010_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxz_xxxx, g_x_0_x_0_xxz_xxxxz, g_x_0_x_0_xxz_xxxy, g_x_0_x_0_xxz_xxxyz, g_x_0_x_0_xxz_xxxz, g_x_0_x_0_xxz_xxxzz, g_x_0_x_0_xxz_xxyy, g_x_0_x_0_xxz_xxyyz, g_x_0_x_0_xxz_xxyz, g_x_0_x_0_xxz_xxyzz, g_x_0_x_0_xxz_xxzz, g_x_0_x_0_xxz_xxzzz, g_x_0_x_0_xxz_xyyy, g_x_0_x_0_xxz_xyyyz, g_x_0_x_0_xxz_xyyz, g_x_0_x_0_xxz_xyyzz, g_x_0_x_0_xxz_xyzz, g_x_0_x_0_xxz_xyzzz, g_x_0_x_0_xxz_xzzz, g_x_0_x_0_xxz_xzzzz, g_x_0_x_0_xxz_yyyy, g_x_0_x_0_xxz_yyyyz, g_x_0_x_0_xxz_yyyz, g_x_0_x_0_xxz_yyyzz, g_x_0_x_0_xxz_yyzz, g_x_0_x_0_xxz_yyzzz, g_x_0_x_0_xxz_yzzz, g_x_0_x_0_xxz_yzzzz, g_x_0_x_0_xxz_zzzz, g_x_0_x_0_xxz_zzzzz, g_x_0_x_0_xxzz_xxxx, g_x_0_x_0_xxzz_xxxy, g_x_0_x_0_xxzz_xxxz, g_x_0_x_0_xxzz_xxyy, g_x_0_x_0_xxzz_xxyz, g_x_0_x_0_xxzz_xxzz, g_x_0_x_0_xxzz_xyyy, g_x_0_x_0_xxzz_xyyz, g_x_0_x_0_xxzz_xyzz, g_x_0_x_0_xxzz_xzzz, g_x_0_x_0_xxzz_yyyy, g_x_0_x_0_xxzz_yyyz, g_x_0_x_0_xxzz_yyzz, g_x_0_x_0_xxzz_yzzz, g_x_0_x_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxzz_xxxx[k] = -g_x_0_x_0_xxz_xxxx[k] * ab_z + g_x_0_x_0_xxz_xxxxz[k];

                g_x_0_x_0_xxzz_xxxy[k] = -g_x_0_x_0_xxz_xxxy[k] * ab_z + g_x_0_x_0_xxz_xxxyz[k];

                g_x_0_x_0_xxzz_xxxz[k] = -g_x_0_x_0_xxz_xxxz[k] * ab_z + g_x_0_x_0_xxz_xxxzz[k];

                g_x_0_x_0_xxzz_xxyy[k] = -g_x_0_x_0_xxz_xxyy[k] * ab_z + g_x_0_x_0_xxz_xxyyz[k];

                g_x_0_x_0_xxzz_xxyz[k] = -g_x_0_x_0_xxz_xxyz[k] * ab_z + g_x_0_x_0_xxz_xxyzz[k];

                g_x_0_x_0_xxzz_xxzz[k] = -g_x_0_x_0_xxz_xxzz[k] * ab_z + g_x_0_x_0_xxz_xxzzz[k];

                g_x_0_x_0_xxzz_xyyy[k] = -g_x_0_x_0_xxz_xyyy[k] * ab_z + g_x_0_x_0_xxz_xyyyz[k];

                g_x_0_x_0_xxzz_xyyz[k] = -g_x_0_x_0_xxz_xyyz[k] * ab_z + g_x_0_x_0_xxz_xyyzz[k];

                g_x_0_x_0_xxzz_xyzz[k] = -g_x_0_x_0_xxz_xyzz[k] * ab_z + g_x_0_x_0_xxz_xyzzz[k];

                g_x_0_x_0_xxzz_xzzz[k] = -g_x_0_x_0_xxz_xzzz[k] * ab_z + g_x_0_x_0_xxz_xzzzz[k];

                g_x_0_x_0_xxzz_yyyy[k] = -g_x_0_x_0_xxz_yyyy[k] * ab_z + g_x_0_x_0_xxz_yyyyz[k];

                g_x_0_x_0_xxzz_yyyz[k] = -g_x_0_x_0_xxz_yyyz[k] * ab_z + g_x_0_x_0_xxz_yyyzz[k];

                g_x_0_x_0_xxzz_yyzz[k] = -g_x_0_x_0_xxz_yyzz[k] * ab_z + g_x_0_x_0_xxz_yyzzz[k];

                g_x_0_x_0_xxzz_yzzz[k] = -g_x_0_x_0_xxz_yzzz[k] * ab_z + g_x_0_x_0_xxz_yzzzz[k];

                g_x_0_x_0_xxzz_zzzz[k] = -g_x_0_x_0_xxz_zzzz[k] * ab_z + g_x_0_x_0_xxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyy_xxxx = cbuffer.data(gg_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_xxxy = cbuffer.data(gg_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_xxxz = cbuffer.data(gg_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_xxyy = cbuffer.data(gg_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_xxyz = cbuffer.data(gg_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_xxzz = cbuffer.data(gg_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_xyyy = cbuffer.data(gg_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_xyyz = cbuffer.data(gg_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_xyzz = cbuffer.data(gg_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_xzzz = cbuffer.data(gg_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_yyyy = cbuffer.data(gg_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_yyyz = cbuffer.data(gg_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_yyzz = cbuffer.data(gg_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_yzzz = cbuffer.data(gg_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_zzzz = cbuffer.data(gg_geom_1010_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyy_xxxx, g_x_0_x_0_xyy_xxxxy, g_x_0_x_0_xyy_xxxy, g_x_0_x_0_xyy_xxxyy, g_x_0_x_0_xyy_xxxyz, g_x_0_x_0_xyy_xxxz, g_x_0_x_0_xyy_xxyy, g_x_0_x_0_xyy_xxyyy, g_x_0_x_0_xyy_xxyyz, g_x_0_x_0_xyy_xxyz, g_x_0_x_0_xyy_xxyzz, g_x_0_x_0_xyy_xxzz, g_x_0_x_0_xyy_xyyy, g_x_0_x_0_xyy_xyyyy, g_x_0_x_0_xyy_xyyyz, g_x_0_x_0_xyy_xyyz, g_x_0_x_0_xyy_xyyzz, g_x_0_x_0_xyy_xyzz, g_x_0_x_0_xyy_xyzzz, g_x_0_x_0_xyy_xzzz, g_x_0_x_0_xyy_yyyy, g_x_0_x_0_xyy_yyyyy, g_x_0_x_0_xyy_yyyyz, g_x_0_x_0_xyy_yyyz, g_x_0_x_0_xyy_yyyzz, g_x_0_x_0_xyy_yyzz, g_x_0_x_0_xyy_yyzzz, g_x_0_x_0_xyy_yzzz, g_x_0_x_0_xyy_yzzzz, g_x_0_x_0_xyy_zzzz, g_x_0_x_0_xyyy_xxxx, g_x_0_x_0_xyyy_xxxy, g_x_0_x_0_xyyy_xxxz, g_x_0_x_0_xyyy_xxyy, g_x_0_x_0_xyyy_xxyz, g_x_0_x_0_xyyy_xxzz, g_x_0_x_0_xyyy_xyyy, g_x_0_x_0_xyyy_xyyz, g_x_0_x_0_xyyy_xyzz, g_x_0_x_0_xyyy_xzzz, g_x_0_x_0_xyyy_yyyy, g_x_0_x_0_xyyy_yyyz, g_x_0_x_0_xyyy_yyzz, g_x_0_x_0_xyyy_yzzz, g_x_0_x_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyy_xxxx[k] = -g_x_0_x_0_xyy_xxxx[k] * ab_y + g_x_0_x_0_xyy_xxxxy[k];

                g_x_0_x_0_xyyy_xxxy[k] = -g_x_0_x_0_xyy_xxxy[k] * ab_y + g_x_0_x_0_xyy_xxxyy[k];

                g_x_0_x_0_xyyy_xxxz[k] = -g_x_0_x_0_xyy_xxxz[k] * ab_y + g_x_0_x_0_xyy_xxxyz[k];

                g_x_0_x_0_xyyy_xxyy[k] = -g_x_0_x_0_xyy_xxyy[k] * ab_y + g_x_0_x_0_xyy_xxyyy[k];

                g_x_0_x_0_xyyy_xxyz[k] = -g_x_0_x_0_xyy_xxyz[k] * ab_y + g_x_0_x_0_xyy_xxyyz[k];

                g_x_0_x_0_xyyy_xxzz[k] = -g_x_0_x_0_xyy_xxzz[k] * ab_y + g_x_0_x_0_xyy_xxyzz[k];

                g_x_0_x_0_xyyy_xyyy[k] = -g_x_0_x_0_xyy_xyyy[k] * ab_y + g_x_0_x_0_xyy_xyyyy[k];

                g_x_0_x_0_xyyy_xyyz[k] = -g_x_0_x_0_xyy_xyyz[k] * ab_y + g_x_0_x_0_xyy_xyyyz[k];

                g_x_0_x_0_xyyy_xyzz[k] = -g_x_0_x_0_xyy_xyzz[k] * ab_y + g_x_0_x_0_xyy_xyyzz[k];

                g_x_0_x_0_xyyy_xzzz[k] = -g_x_0_x_0_xyy_xzzz[k] * ab_y + g_x_0_x_0_xyy_xyzzz[k];

                g_x_0_x_0_xyyy_yyyy[k] = -g_x_0_x_0_xyy_yyyy[k] * ab_y + g_x_0_x_0_xyy_yyyyy[k];

                g_x_0_x_0_xyyy_yyyz[k] = -g_x_0_x_0_xyy_yyyz[k] * ab_y + g_x_0_x_0_xyy_yyyyz[k];

                g_x_0_x_0_xyyy_yyzz[k] = -g_x_0_x_0_xyy_yyzz[k] * ab_y + g_x_0_x_0_xyy_yyyzz[k];

                g_x_0_x_0_xyyy_yzzz[k] = -g_x_0_x_0_xyy_yzzz[k] * ab_y + g_x_0_x_0_xyy_yyzzz[k];

                g_x_0_x_0_xyyy_zzzz[k] = -g_x_0_x_0_xyy_zzzz[k] * ab_y + g_x_0_x_0_xyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyz_xxxx = cbuffer.data(gg_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_xxxy = cbuffer.data(gg_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_xxxz = cbuffer.data(gg_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_xxyy = cbuffer.data(gg_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_xxyz = cbuffer.data(gg_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_xxzz = cbuffer.data(gg_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_xyyy = cbuffer.data(gg_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_xyyz = cbuffer.data(gg_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_xyzz = cbuffer.data(gg_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_xzzz = cbuffer.data(gg_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_yyyy = cbuffer.data(gg_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_yyyz = cbuffer.data(gg_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_yyzz = cbuffer.data(gg_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_yzzz = cbuffer.data(gg_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_zzzz = cbuffer.data(gg_geom_1010_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyz_xxxx, g_x_0_x_0_xyyz_xxxy, g_x_0_x_0_xyyz_xxxz, g_x_0_x_0_xyyz_xxyy, g_x_0_x_0_xyyz_xxyz, g_x_0_x_0_xyyz_xxzz, g_x_0_x_0_xyyz_xyyy, g_x_0_x_0_xyyz_xyyz, g_x_0_x_0_xyyz_xyzz, g_x_0_x_0_xyyz_xzzz, g_x_0_x_0_xyyz_yyyy, g_x_0_x_0_xyyz_yyyz, g_x_0_x_0_xyyz_yyzz, g_x_0_x_0_xyyz_yzzz, g_x_0_x_0_xyyz_zzzz, g_x_0_x_0_xyz_xxxx, g_x_0_x_0_xyz_xxxxy, g_x_0_x_0_xyz_xxxy, g_x_0_x_0_xyz_xxxyy, g_x_0_x_0_xyz_xxxyz, g_x_0_x_0_xyz_xxxz, g_x_0_x_0_xyz_xxyy, g_x_0_x_0_xyz_xxyyy, g_x_0_x_0_xyz_xxyyz, g_x_0_x_0_xyz_xxyz, g_x_0_x_0_xyz_xxyzz, g_x_0_x_0_xyz_xxzz, g_x_0_x_0_xyz_xyyy, g_x_0_x_0_xyz_xyyyy, g_x_0_x_0_xyz_xyyyz, g_x_0_x_0_xyz_xyyz, g_x_0_x_0_xyz_xyyzz, g_x_0_x_0_xyz_xyzz, g_x_0_x_0_xyz_xyzzz, g_x_0_x_0_xyz_xzzz, g_x_0_x_0_xyz_yyyy, g_x_0_x_0_xyz_yyyyy, g_x_0_x_0_xyz_yyyyz, g_x_0_x_0_xyz_yyyz, g_x_0_x_0_xyz_yyyzz, g_x_0_x_0_xyz_yyzz, g_x_0_x_0_xyz_yyzzz, g_x_0_x_0_xyz_yzzz, g_x_0_x_0_xyz_yzzzz, g_x_0_x_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyz_xxxx[k] = -g_x_0_x_0_xyz_xxxx[k] * ab_y + g_x_0_x_0_xyz_xxxxy[k];

                g_x_0_x_0_xyyz_xxxy[k] = -g_x_0_x_0_xyz_xxxy[k] * ab_y + g_x_0_x_0_xyz_xxxyy[k];

                g_x_0_x_0_xyyz_xxxz[k] = -g_x_0_x_0_xyz_xxxz[k] * ab_y + g_x_0_x_0_xyz_xxxyz[k];

                g_x_0_x_0_xyyz_xxyy[k] = -g_x_0_x_0_xyz_xxyy[k] * ab_y + g_x_0_x_0_xyz_xxyyy[k];

                g_x_0_x_0_xyyz_xxyz[k] = -g_x_0_x_0_xyz_xxyz[k] * ab_y + g_x_0_x_0_xyz_xxyyz[k];

                g_x_0_x_0_xyyz_xxzz[k] = -g_x_0_x_0_xyz_xxzz[k] * ab_y + g_x_0_x_0_xyz_xxyzz[k];

                g_x_0_x_0_xyyz_xyyy[k] = -g_x_0_x_0_xyz_xyyy[k] * ab_y + g_x_0_x_0_xyz_xyyyy[k];

                g_x_0_x_0_xyyz_xyyz[k] = -g_x_0_x_0_xyz_xyyz[k] * ab_y + g_x_0_x_0_xyz_xyyyz[k];

                g_x_0_x_0_xyyz_xyzz[k] = -g_x_0_x_0_xyz_xyzz[k] * ab_y + g_x_0_x_0_xyz_xyyzz[k];

                g_x_0_x_0_xyyz_xzzz[k] = -g_x_0_x_0_xyz_xzzz[k] * ab_y + g_x_0_x_0_xyz_xyzzz[k];

                g_x_0_x_0_xyyz_yyyy[k] = -g_x_0_x_0_xyz_yyyy[k] * ab_y + g_x_0_x_0_xyz_yyyyy[k];

                g_x_0_x_0_xyyz_yyyz[k] = -g_x_0_x_0_xyz_yyyz[k] * ab_y + g_x_0_x_0_xyz_yyyyz[k];

                g_x_0_x_0_xyyz_yyzz[k] = -g_x_0_x_0_xyz_yyzz[k] * ab_y + g_x_0_x_0_xyz_yyyzz[k];

                g_x_0_x_0_xyyz_yzzz[k] = -g_x_0_x_0_xyz_yzzz[k] * ab_y + g_x_0_x_0_xyz_yyzzz[k];

                g_x_0_x_0_xyyz_zzzz[k] = -g_x_0_x_0_xyz_zzzz[k] * ab_y + g_x_0_x_0_xyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyzz_xxxx = cbuffer.data(gg_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_xxxy = cbuffer.data(gg_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_xxxz = cbuffer.data(gg_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_xxyy = cbuffer.data(gg_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_xxyz = cbuffer.data(gg_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_xxzz = cbuffer.data(gg_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_xyyy = cbuffer.data(gg_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_xyyz = cbuffer.data(gg_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_xyzz = cbuffer.data(gg_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_xzzz = cbuffer.data(gg_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_yyyy = cbuffer.data(gg_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_yyyz = cbuffer.data(gg_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_yyzz = cbuffer.data(gg_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_yzzz = cbuffer.data(gg_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_zzzz = cbuffer.data(gg_geom_1010_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyzz_xxxx, g_x_0_x_0_xyzz_xxxy, g_x_0_x_0_xyzz_xxxz, g_x_0_x_0_xyzz_xxyy, g_x_0_x_0_xyzz_xxyz, g_x_0_x_0_xyzz_xxzz, g_x_0_x_0_xyzz_xyyy, g_x_0_x_0_xyzz_xyyz, g_x_0_x_0_xyzz_xyzz, g_x_0_x_0_xyzz_xzzz, g_x_0_x_0_xyzz_yyyy, g_x_0_x_0_xyzz_yyyz, g_x_0_x_0_xyzz_yyzz, g_x_0_x_0_xyzz_yzzz, g_x_0_x_0_xyzz_zzzz, g_x_0_x_0_xzz_xxxx, g_x_0_x_0_xzz_xxxxy, g_x_0_x_0_xzz_xxxy, g_x_0_x_0_xzz_xxxyy, g_x_0_x_0_xzz_xxxyz, g_x_0_x_0_xzz_xxxz, g_x_0_x_0_xzz_xxyy, g_x_0_x_0_xzz_xxyyy, g_x_0_x_0_xzz_xxyyz, g_x_0_x_0_xzz_xxyz, g_x_0_x_0_xzz_xxyzz, g_x_0_x_0_xzz_xxzz, g_x_0_x_0_xzz_xyyy, g_x_0_x_0_xzz_xyyyy, g_x_0_x_0_xzz_xyyyz, g_x_0_x_0_xzz_xyyz, g_x_0_x_0_xzz_xyyzz, g_x_0_x_0_xzz_xyzz, g_x_0_x_0_xzz_xyzzz, g_x_0_x_0_xzz_xzzz, g_x_0_x_0_xzz_yyyy, g_x_0_x_0_xzz_yyyyy, g_x_0_x_0_xzz_yyyyz, g_x_0_x_0_xzz_yyyz, g_x_0_x_0_xzz_yyyzz, g_x_0_x_0_xzz_yyzz, g_x_0_x_0_xzz_yyzzz, g_x_0_x_0_xzz_yzzz, g_x_0_x_0_xzz_yzzzz, g_x_0_x_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyzz_xxxx[k] = -g_x_0_x_0_xzz_xxxx[k] * ab_y + g_x_0_x_0_xzz_xxxxy[k];

                g_x_0_x_0_xyzz_xxxy[k] = -g_x_0_x_0_xzz_xxxy[k] * ab_y + g_x_0_x_0_xzz_xxxyy[k];

                g_x_0_x_0_xyzz_xxxz[k] = -g_x_0_x_0_xzz_xxxz[k] * ab_y + g_x_0_x_0_xzz_xxxyz[k];

                g_x_0_x_0_xyzz_xxyy[k] = -g_x_0_x_0_xzz_xxyy[k] * ab_y + g_x_0_x_0_xzz_xxyyy[k];

                g_x_0_x_0_xyzz_xxyz[k] = -g_x_0_x_0_xzz_xxyz[k] * ab_y + g_x_0_x_0_xzz_xxyyz[k];

                g_x_0_x_0_xyzz_xxzz[k] = -g_x_0_x_0_xzz_xxzz[k] * ab_y + g_x_0_x_0_xzz_xxyzz[k];

                g_x_0_x_0_xyzz_xyyy[k] = -g_x_0_x_0_xzz_xyyy[k] * ab_y + g_x_0_x_0_xzz_xyyyy[k];

                g_x_0_x_0_xyzz_xyyz[k] = -g_x_0_x_0_xzz_xyyz[k] * ab_y + g_x_0_x_0_xzz_xyyyz[k];

                g_x_0_x_0_xyzz_xyzz[k] = -g_x_0_x_0_xzz_xyzz[k] * ab_y + g_x_0_x_0_xzz_xyyzz[k];

                g_x_0_x_0_xyzz_xzzz[k] = -g_x_0_x_0_xzz_xzzz[k] * ab_y + g_x_0_x_0_xzz_xyzzz[k];

                g_x_0_x_0_xyzz_yyyy[k] = -g_x_0_x_0_xzz_yyyy[k] * ab_y + g_x_0_x_0_xzz_yyyyy[k];

                g_x_0_x_0_xyzz_yyyz[k] = -g_x_0_x_0_xzz_yyyz[k] * ab_y + g_x_0_x_0_xzz_yyyyz[k];

                g_x_0_x_0_xyzz_yyzz[k] = -g_x_0_x_0_xzz_yyzz[k] * ab_y + g_x_0_x_0_xzz_yyyzz[k];

                g_x_0_x_0_xyzz_yzzz[k] = -g_x_0_x_0_xzz_yzzz[k] * ab_y + g_x_0_x_0_xzz_yyzzz[k];

                g_x_0_x_0_xyzz_zzzz[k] = -g_x_0_x_0_xzz_zzzz[k] * ab_y + g_x_0_x_0_xzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xzzz_xxxx = cbuffer.data(gg_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_xxxy = cbuffer.data(gg_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_xxxz = cbuffer.data(gg_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_xxyy = cbuffer.data(gg_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_xxyz = cbuffer.data(gg_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_xxzz = cbuffer.data(gg_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_xyyy = cbuffer.data(gg_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_xyyz = cbuffer.data(gg_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_xyzz = cbuffer.data(gg_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_xzzz = cbuffer.data(gg_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_yyyy = cbuffer.data(gg_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_yyyz = cbuffer.data(gg_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_yyzz = cbuffer.data(gg_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_yzzz = cbuffer.data(gg_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_zzzz = cbuffer.data(gg_geom_1010_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xzz_xxxx, g_x_0_x_0_xzz_xxxxz, g_x_0_x_0_xzz_xxxy, g_x_0_x_0_xzz_xxxyz, g_x_0_x_0_xzz_xxxz, g_x_0_x_0_xzz_xxxzz, g_x_0_x_0_xzz_xxyy, g_x_0_x_0_xzz_xxyyz, g_x_0_x_0_xzz_xxyz, g_x_0_x_0_xzz_xxyzz, g_x_0_x_0_xzz_xxzz, g_x_0_x_0_xzz_xxzzz, g_x_0_x_0_xzz_xyyy, g_x_0_x_0_xzz_xyyyz, g_x_0_x_0_xzz_xyyz, g_x_0_x_0_xzz_xyyzz, g_x_0_x_0_xzz_xyzz, g_x_0_x_0_xzz_xyzzz, g_x_0_x_0_xzz_xzzz, g_x_0_x_0_xzz_xzzzz, g_x_0_x_0_xzz_yyyy, g_x_0_x_0_xzz_yyyyz, g_x_0_x_0_xzz_yyyz, g_x_0_x_0_xzz_yyyzz, g_x_0_x_0_xzz_yyzz, g_x_0_x_0_xzz_yyzzz, g_x_0_x_0_xzz_yzzz, g_x_0_x_0_xzz_yzzzz, g_x_0_x_0_xzz_zzzz, g_x_0_x_0_xzz_zzzzz, g_x_0_x_0_xzzz_xxxx, g_x_0_x_0_xzzz_xxxy, g_x_0_x_0_xzzz_xxxz, g_x_0_x_0_xzzz_xxyy, g_x_0_x_0_xzzz_xxyz, g_x_0_x_0_xzzz_xxzz, g_x_0_x_0_xzzz_xyyy, g_x_0_x_0_xzzz_xyyz, g_x_0_x_0_xzzz_xyzz, g_x_0_x_0_xzzz_xzzz, g_x_0_x_0_xzzz_yyyy, g_x_0_x_0_xzzz_yyyz, g_x_0_x_0_xzzz_yyzz, g_x_0_x_0_xzzz_yzzz, g_x_0_x_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xzzz_xxxx[k] = -g_x_0_x_0_xzz_xxxx[k] * ab_z + g_x_0_x_0_xzz_xxxxz[k];

                g_x_0_x_0_xzzz_xxxy[k] = -g_x_0_x_0_xzz_xxxy[k] * ab_z + g_x_0_x_0_xzz_xxxyz[k];

                g_x_0_x_0_xzzz_xxxz[k] = -g_x_0_x_0_xzz_xxxz[k] * ab_z + g_x_0_x_0_xzz_xxxzz[k];

                g_x_0_x_0_xzzz_xxyy[k] = -g_x_0_x_0_xzz_xxyy[k] * ab_z + g_x_0_x_0_xzz_xxyyz[k];

                g_x_0_x_0_xzzz_xxyz[k] = -g_x_0_x_0_xzz_xxyz[k] * ab_z + g_x_0_x_0_xzz_xxyzz[k];

                g_x_0_x_0_xzzz_xxzz[k] = -g_x_0_x_0_xzz_xxzz[k] * ab_z + g_x_0_x_0_xzz_xxzzz[k];

                g_x_0_x_0_xzzz_xyyy[k] = -g_x_0_x_0_xzz_xyyy[k] * ab_z + g_x_0_x_0_xzz_xyyyz[k];

                g_x_0_x_0_xzzz_xyyz[k] = -g_x_0_x_0_xzz_xyyz[k] * ab_z + g_x_0_x_0_xzz_xyyzz[k];

                g_x_0_x_0_xzzz_xyzz[k] = -g_x_0_x_0_xzz_xyzz[k] * ab_z + g_x_0_x_0_xzz_xyzzz[k];

                g_x_0_x_0_xzzz_xzzz[k] = -g_x_0_x_0_xzz_xzzz[k] * ab_z + g_x_0_x_0_xzz_xzzzz[k];

                g_x_0_x_0_xzzz_yyyy[k] = -g_x_0_x_0_xzz_yyyy[k] * ab_z + g_x_0_x_0_xzz_yyyyz[k];

                g_x_0_x_0_xzzz_yyyz[k] = -g_x_0_x_0_xzz_yyyz[k] * ab_z + g_x_0_x_0_xzz_yyyzz[k];

                g_x_0_x_0_xzzz_yyzz[k] = -g_x_0_x_0_xzz_yyzz[k] * ab_z + g_x_0_x_0_xzz_yyzzz[k];

                g_x_0_x_0_xzzz_yzzz[k] = -g_x_0_x_0_xzz_yzzz[k] * ab_z + g_x_0_x_0_xzz_yzzzz[k];

                g_x_0_x_0_xzzz_zzzz[k] = -g_x_0_x_0_xzz_zzzz[k] * ab_z + g_x_0_x_0_xzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyy_xxxx = cbuffer.data(gg_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_xxxy = cbuffer.data(gg_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_xxxz = cbuffer.data(gg_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_xxyy = cbuffer.data(gg_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_xxyz = cbuffer.data(gg_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_xxzz = cbuffer.data(gg_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_xyyy = cbuffer.data(gg_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_xyyz = cbuffer.data(gg_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_xyzz = cbuffer.data(gg_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_xzzz = cbuffer.data(gg_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_yyyy = cbuffer.data(gg_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_yyyz = cbuffer.data(gg_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_yyzz = cbuffer.data(gg_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_yzzz = cbuffer.data(gg_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_zzzz = cbuffer.data(gg_geom_1010_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyy_xxxx, g_x_0_x_0_yyy_xxxxy, g_x_0_x_0_yyy_xxxy, g_x_0_x_0_yyy_xxxyy, g_x_0_x_0_yyy_xxxyz, g_x_0_x_0_yyy_xxxz, g_x_0_x_0_yyy_xxyy, g_x_0_x_0_yyy_xxyyy, g_x_0_x_0_yyy_xxyyz, g_x_0_x_0_yyy_xxyz, g_x_0_x_0_yyy_xxyzz, g_x_0_x_0_yyy_xxzz, g_x_0_x_0_yyy_xyyy, g_x_0_x_0_yyy_xyyyy, g_x_0_x_0_yyy_xyyyz, g_x_0_x_0_yyy_xyyz, g_x_0_x_0_yyy_xyyzz, g_x_0_x_0_yyy_xyzz, g_x_0_x_0_yyy_xyzzz, g_x_0_x_0_yyy_xzzz, g_x_0_x_0_yyy_yyyy, g_x_0_x_0_yyy_yyyyy, g_x_0_x_0_yyy_yyyyz, g_x_0_x_0_yyy_yyyz, g_x_0_x_0_yyy_yyyzz, g_x_0_x_0_yyy_yyzz, g_x_0_x_0_yyy_yyzzz, g_x_0_x_0_yyy_yzzz, g_x_0_x_0_yyy_yzzzz, g_x_0_x_0_yyy_zzzz, g_x_0_x_0_yyyy_xxxx, g_x_0_x_0_yyyy_xxxy, g_x_0_x_0_yyyy_xxxz, g_x_0_x_0_yyyy_xxyy, g_x_0_x_0_yyyy_xxyz, g_x_0_x_0_yyyy_xxzz, g_x_0_x_0_yyyy_xyyy, g_x_0_x_0_yyyy_xyyz, g_x_0_x_0_yyyy_xyzz, g_x_0_x_0_yyyy_xzzz, g_x_0_x_0_yyyy_yyyy, g_x_0_x_0_yyyy_yyyz, g_x_0_x_0_yyyy_yyzz, g_x_0_x_0_yyyy_yzzz, g_x_0_x_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyy_xxxx[k] = -g_x_0_x_0_yyy_xxxx[k] * ab_y + g_x_0_x_0_yyy_xxxxy[k];

                g_x_0_x_0_yyyy_xxxy[k] = -g_x_0_x_0_yyy_xxxy[k] * ab_y + g_x_0_x_0_yyy_xxxyy[k];

                g_x_0_x_0_yyyy_xxxz[k] = -g_x_0_x_0_yyy_xxxz[k] * ab_y + g_x_0_x_0_yyy_xxxyz[k];

                g_x_0_x_0_yyyy_xxyy[k] = -g_x_0_x_0_yyy_xxyy[k] * ab_y + g_x_0_x_0_yyy_xxyyy[k];

                g_x_0_x_0_yyyy_xxyz[k] = -g_x_0_x_0_yyy_xxyz[k] * ab_y + g_x_0_x_0_yyy_xxyyz[k];

                g_x_0_x_0_yyyy_xxzz[k] = -g_x_0_x_0_yyy_xxzz[k] * ab_y + g_x_0_x_0_yyy_xxyzz[k];

                g_x_0_x_0_yyyy_xyyy[k] = -g_x_0_x_0_yyy_xyyy[k] * ab_y + g_x_0_x_0_yyy_xyyyy[k];

                g_x_0_x_0_yyyy_xyyz[k] = -g_x_0_x_0_yyy_xyyz[k] * ab_y + g_x_0_x_0_yyy_xyyyz[k];

                g_x_0_x_0_yyyy_xyzz[k] = -g_x_0_x_0_yyy_xyzz[k] * ab_y + g_x_0_x_0_yyy_xyyzz[k];

                g_x_0_x_0_yyyy_xzzz[k] = -g_x_0_x_0_yyy_xzzz[k] * ab_y + g_x_0_x_0_yyy_xyzzz[k];

                g_x_0_x_0_yyyy_yyyy[k] = -g_x_0_x_0_yyy_yyyy[k] * ab_y + g_x_0_x_0_yyy_yyyyy[k];

                g_x_0_x_0_yyyy_yyyz[k] = -g_x_0_x_0_yyy_yyyz[k] * ab_y + g_x_0_x_0_yyy_yyyyz[k];

                g_x_0_x_0_yyyy_yyzz[k] = -g_x_0_x_0_yyy_yyzz[k] * ab_y + g_x_0_x_0_yyy_yyyzz[k];

                g_x_0_x_0_yyyy_yzzz[k] = -g_x_0_x_0_yyy_yzzz[k] * ab_y + g_x_0_x_0_yyy_yyzzz[k];

                g_x_0_x_0_yyyy_zzzz[k] = -g_x_0_x_0_yyy_zzzz[k] * ab_y + g_x_0_x_0_yyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyz_xxxx = cbuffer.data(gg_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_xxxy = cbuffer.data(gg_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_xxxz = cbuffer.data(gg_geom_1010_off + 167 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_xxyy = cbuffer.data(gg_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_xxyz = cbuffer.data(gg_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_xxzz = cbuffer.data(gg_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_xyyy = cbuffer.data(gg_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_xyyz = cbuffer.data(gg_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_xyzz = cbuffer.data(gg_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_xzzz = cbuffer.data(gg_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_yyyy = cbuffer.data(gg_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_yyyz = cbuffer.data(gg_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_yyzz = cbuffer.data(gg_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_yzzz = cbuffer.data(gg_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_zzzz = cbuffer.data(gg_geom_1010_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyz_xxxx, g_x_0_x_0_yyyz_xxxy, g_x_0_x_0_yyyz_xxxz, g_x_0_x_0_yyyz_xxyy, g_x_0_x_0_yyyz_xxyz, g_x_0_x_0_yyyz_xxzz, g_x_0_x_0_yyyz_xyyy, g_x_0_x_0_yyyz_xyyz, g_x_0_x_0_yyyz_xyzz, g_x_0_x_0_yyyz_xzzz, g_x_0_x_0_yyyz_yyyy, g_x_0_x_0_yyyz_yyyz, g_x_0_x_0_yyyz_yyzz, g_x_0_x_0_yyyz_yzzz, g_x_0_x_0_yyyz_zzzz, g_x_0_x_0_yyz_xxxx, g_x_0_x_0_yyz_xxxxy, g_x_0_x_0_yyz_xxxy, g_x_0_x_0_yyz_xxxyy, g_x_0_x_0_yyz_xxxyz, g_x_0_x_0_yyz_xxxz, g_x_0_x_0_yyz_xxyy, g_x_0_x_0_yyz_xxyyy, g_x_0_x_0_yyz_xxyyz, g_x_0_x_0_yyz_xxyz, g_x_0_x_0_yyz_xxyzz, g_x_0_x_0_yyz_xxzz, g_x_0_x_0_yyz_xyyy, g_x_0_x_0_yyz_xyyyy, g_x_0_x_0_yyz_xyyyz, g_x_0_x_0_yyz_xyyz, g_x_0_x_0_yyz_xyyzz, g_x_0_x_0_yyz_xyzz, g_x_0_x_0_yyz_xyzzz, g_x_0_x_0_yyz_xzzz, g_x_0_x_0_yyz_yyyy, g_x_0_x_0_yyz_yyyyy, g_x_0_x_0_yyz_yyyyz, g_x_0_x_0_yyz_yyyz, g_x_0_x_0_yyz_yyyzz, g_x_0_x_0_yyz_yyzz, g_x_0_x_0_yyz_yyzzz, g_x_0_x_0_yyz_yzzz, g_x_0_x_0_yyz_yzzzz, g_x_0_x_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyz_xxxx[k] = -g_x_0_x_0_yyz_xxxx[k] * ab_y + g_x_0_x_0_yyz_xxxxy[k];

                g_x_0_x_0_yyyz_xxxy[k] = -g_x_0_x_0_yyz_xxxy[k] * ab_y + g_x_0_x_0_yyz_xxxyy[k];

                g_x_0_x_0_yyyz_xxxz[k] = -g_x_0_x_0_yyz_xxxz[k] * ab_y + g_x_0_x_0_yyz_xxxyz[k];

                g_x_0_x_0_yyyz_xxyy[k] = -g_x_0_x_0_yyz_xxyy[k] * ab_y + g_x_0_x_0_yyz_xxyyy[k];

                g_x_0_x_0_yyyz_xxyz[k] = -g_x_0_x_0_yyz_xxyz[k] * ab_y + g_x_0_x_0_yyz_xxyyz[k];

                g_x_0_x_0_yyyz_xxzz[k] = -g_x_0_x_0_yyz_xxzz[k] * ab_y + g_x_0_x_0_yyz_xxyzz[k];

                g_x_0_x_0_yyyz_xyyy[k] = -g_x_0_x_0_yyz_xyyy[k] * ab_y + g_x_0_x_0_yyz_xyyyy[k];

                g_x_0_x_0_yyyz_xyyz[k] = -g_x_0_x_0_yyz_xyyz[k] * ab_y + g_x_0_x_0_yyz_xyyyz[k];

                g_x_0_x_0_yyyz_xyzz[k] = -g_x_0_x_0_yyz_xyzz[k] * ab_y + g_x_0_x_0_yyz_xyyzz[k];

                g_x_0_x_0_yyyz_xzzz[k] = -g_x_0_x_0_yyz_xzzz[k] * ab_y + g_x_0_x_0_yyz_xyzzz[k];

                g_x_0_x_0_yyyz_yyyy[k] = -g_x_0_x_0_yyz_yyyy[k] * ab_y + g_x_0_x_0_yyz_yyyyy[k];

                g_x_0_x_0_yyyz_yyyz[k] = -g_x_0_x_0_yyz_yyyz[k] * ab_y + g_x_0_x_0_yyz_yyyyz[k];

                g_x_0_x_0_yyyz_yyzz[k] = -g_x_0_x_0_yyz_yyzz[k] * ab_y + g_x_0_x_0_yyz_yyyzz[k];

                g_x_0_x_0_yyyz_yzzz[k] = -g_x_0_x_0_yyz_yzzz[k] * ab_y + g_x_0_x_0_yyz_yyzzz[k];

                g_x_0_x_0_yyyz_zzzz[k] = -g_x_0_x_0_yyz_zzzz[k] * ab_y + g_x_0_x_0_yyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyzz_xxxx = cbuffer.data(gg_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_xxxy = cbuffer.data(gg_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_xxxz = cbuffer.data(gg_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_xxyy = cbuffer.data(gg_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_xxyz = cbuffer.data(gg_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_xxzz = cbuffer.data(gg_geom_1010_off + 185 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_xyyy = cbuffer.data(gg_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_xyyz = cbuffer.data(gg_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_xyzz = cbuffer.data(gg_geom_1010_off + 188 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_xzzz = cbuffer.data(gg_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_yyyy = cbuffer.data(gg_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_yyyz = cbuffer.data(gg_geom_1010_off + 191 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_yyzz = cbuffer.data(gg_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_yzzz = cbuffer.data(gg_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_zzzz = cbuffer.data(gg_geom_1010_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyzz_xxxx, g_x_0_x_0_yyzz_xxxy, g_x_0_x_0_yyzz_xxxz, g_x_0_x_0_yyzz_xxyy, g_x_0_x_0_yyzz_xxyz, g_x_0_x_0_yyzz_xxzz, g_x_0_x_0_yyzz_xyyy, g_x_0_x_0_yyzz_xyyz, g_x_0_x_0_yyzz_xyzz, g_x_0_x_0_yyzz_xzzz, g_x_0_x_0_yyzz_yyyy, g_x_0_x_0_yyzz_yyyz, g_x_0_x_0_yyzz_yyzz, g_x_0_x_0_yyzz_yzzz, g_x_0_x_0_yyzz_zzzz, g_x_0_x_0_yzz_xxxx, g_x_0_x_0_yzz_xxxxy, g_x_0_x_0_yzz_xxxy, g_x_0_x_0_yzz_xxxyy, g_x_0_x_0_yzz_xxxyz, g_x_0_x_0_yzz_xxxz, g_x_0_x_0_yzz_xxyy, g_x_0_x_0_yzz_xxyyy, g_x_0_x_0_yzz_xxyyz, g_x_0_x_0_yzz_xxyz, g_x_0_x_0_yzz_xxyzz, g_x_0_x_0_yzz_xxzz, g_x_0_x_0_yzz_xyyy, g_x_0_x_0_yzz_xyyyy, g_x_0_x_0_yzz_xyyyz, g_x_0_x_0_yzz_xyyz, g_x_0_x_0_yzz_xyyzz, g_x_0_x_0_yzz_xyzz, g_x_0_x_0_yzz_xyzzz, g_x_0_x_0_yzz_xzzz, g_x_0_x_0_yzz_yyyy, g_x_0_x_0_yzz_yyyyy, g_x_0_x_0_yzz_yyyyz, g_x_0_x_0_yzz_yyyz, g_x_0_x_0_yzz_yyyzz, g_x_0_x_0_yzz_yyzz, g_x_0_x_0_yzz_yyzzz, g_x_0_x_0_yzz_yzzz, g_x_0_x_0_yzz_yzzzz, g_x_0_x_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyzz_xxxx[k] = -g_x_0_x_0_yzz_xxxx[k] * ab_y + g_x_0_x_0_yzz_xxxxy[k];

                g_x_0_x_0_yyzz_xxxy[k] = -g_x_0_x_0_yzz_xxxy[k] * ab_y + g_x_0_x_0_yzz_xxxyy[k];

                g_x_0_x_0_yyzz_xxxz[k] = -g_x_0_x_0_yzz_xxxz[k] * ab_y + g_x_0_x_0_yzz_xxxyz[k];

                g_x_0_x_0_yyzz_xxyy[k] = -g_x_0_x_0_yzz_xxyy[k] * ab_y + g_x_0_x_0_yzz_xxyyy[k];

                g_x_0_x_0_yyzz_xxyz[k] = -g_x_0_x_0_yzz_xxyz[k] * ab_y + g_x_0_x_0_yzz_xxyyz[k];

                g_x_0_x_0_yyzz_xxzz[k] = -g_x_0_x_0_yzz_xxzz[k] * ab_y + g_x_0_x_0_yzz_xxyzz[k];

                g_x_0_x_0_yyzz_xyyy[k] = -g_x_0_x_0_yzz_xyyy[k] * ab_y + g_x_0_x_0_yzz_xyyyy[k];

                g_x_0_x_0_yyzz_xyyz[k] = -g_x_0_x_0_yzz_xyyz[k] * ab_y + g_x_0_x_0_yzz_xyyyz[k];

                g_x_0_x_0_yyzz_xyzz[k] = -g_x_0_x_0_yzz_xyzz[k] * ab_y + g_x_0_x_0_yzz_xyyzz[k];

                g_x_0_x_0_yyzz_xzzz[k] = -g_x_0_x_0_yzz_xzzz[k] * ab_y + g_x_0_x_0_yzz_xyzzz[k];

                g_x_0_x_0_yyzz_yyyy[k] = -g_x_0_x_0_yzz_yyyy[k] * ab_y + g_x_0_x_0_yzz_yyyyy[k];

                g_x_0_x_0_yyzz_yyyz[k] = -g_x_0_x_0_yzz_yyyz[k] * ab_y + g_x_0_x_0_yzz_yyyyz[k];

                g_x_0_x_0_yyzz_yyzz[k] = -g_x_0_x_0_yzz_yyzz[k] * ab_y + g_x_0_x_0_yzz_yyyzz[k];

                g_x_0_x_0_yyzz_yzzz[k] = -g_x_0_x_0_yzz_yzzz[k] * ab_y + g_x_0_x_0_yzz_yyzzz[k];

                g_x_0_x_0_yyzz_zzzz[k] = -g_x_0_x_0_yzz_zzzz[k] * ab_y + g_x_0_x_0_yzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yzzz_xxxx = cbuffer.data(gg_geom_1010_off + 195 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_xxxy = cbuffer.data(gg_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_xxxz = cbuffer.data(gg_geom_1010_off + 197 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_xxyy = cbuffer.data(gg_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_xxyz = cbuffer.data(gg_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_xxzz = cbuffer.data(gg_geom_1010_off + 200 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_xyyy = cbuffer.data(gg_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_xyyz = cbuffer.data(gg_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_xyzz = cbuffer.data(gg_geom_1010_off + 203 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_xzzz = cbuffer.data(gg_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_yyyy = cbuffer.data(gg_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_yyyz = cbuffer.data(gg_geom_1010_off + 206 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_yyzz = cbuffer.data(gg_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_yzzz = cbuffer.data(gg_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_zzzz = cbuffer.data(gg_geom_1010_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yzzz_xxxx, g_x_0_x_0_yzzz_xxxy, g_x_0_x_0_yzzz_xxxz, g_x_0_x_0_yzzz_xxyy, g_x_0_x_0_yzzz_xxyz, g_x_0_x_0_yzzz_xxzz, g_x_0_x_0_yzzz_xyyy, g_x_0_x_0_yzzz_xyyz, g_x_0_x_0_yzzz_xyzz, g_x_0_x_0_yzzz_xzzz, g_x_0_x_0_yzzz_yyyy, g_x_0_x_0_yzzz_yyyz, g_x_0_x_0_yzzz_yyzz, g_x_0_x_0_yzzz_yzzz, g_x_0_x_0_yzzz_zzzz, g_x_0_x_0_zzz_xxxx, g_x_0_x_0_zzz_xxxxy, g_x_0_x_0_zzz_xxxy, g_x_0_x_0_zzz_xxxyy, g_x_0_x_0_zzz_xxxyz, g_x_0_x_0_zzz_xxxz, g_x_0_x_0_zzz_xxyy, g_x_0_x_0_zzz_xxyyy, g_x_0_x_0_zzz_xxyyz, g_x_0_x_0_zzz_xxyz, g_x_0_x_0_zzz_xxyzz, g_x_0_x_0_zzz_xxzz, g_x_0_x_0_zzz_xyyy, g_x_0_x_0_zzz_xyyyy, g_x_0_x_0_zzz_xyyyz, g_x_0_x_0_zzz_xyyz, g_x_0_x_0_zzz_xyyzz, g_x_0_x_0_zzz_xyzz, g_x_0_x_0_zzz_xyzzz, g_x_0_x_0_zzz_xzzz, g_x_0_x_0_zzz_yyyy, g_x_0_x_0_zzz_yyyyy, g_x_0_x_0_zzz_yyyyz, g_x_0_x_0_zzz_yyyz, g_x_0_x_0_zzz_yyyzz, g_x_0_x_0_zzz_yyzz, g_x_0_x_0_zzz_yyzzz, g_x_0_x_0_zzz_yzzz, g_x_0_x_0_zzz_yzzzz, g_x_0_x_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yzzz_xxxx[k] = -g_x_0_x_0_zzz_xxxx[k] * ab_y + g_x_0_x_0_zzz_xxxxy[k];

                g_x_0_x_0_yzzz_xxxy[k] = -g_x_0_x_0_zzz_xxxy[k] * ab_y + g_x_0_x_0_zzz_xxxyy[k];

                g_x_0_x_0_yzzz_xxxz[k] = -g_x_0_x_0_zzz_xxxz[k] * ab_y + g_x_0_x_0_zzz_xxxyz[k];

                g_x_0_x_0_yzzz_xxyy[k] = -g_x_0_x_0_zzz_xxyy[k] * ab_y + g_x_0_x_0_zzz_xxyyy[k];

                g_x_0_x_0_yzzz_xxyz[k] = -g_x_0_x_0_zzz_xxyz[k] * ab_y + g_x_0_x_0_zzz_xxyyz[k];

                g_x_0_x_0_yzzz_xxzz[k] = -g_x_0_x_0_zzz_xxzz[k] * ab_y + g_x_0_x_0_zzz_xxyzz[k];

                g_x_0_x_0_yzzz_xyyy[k] = -g_x_0_x_0_zzz_xyyy[k] * ab_y + g_x_0_x_0_zzz_xyyyy[k];

                g_x_0_x_0_yzzz_xyyz[k] = -g_x_0_x_0_zzz_xyyz[k] * ab_y + g_x_0_x_0_zzz_xyyyz[k];

                g_x_0_x_0_yzzz_xyzz[k] = -g_x_0_x_0_zzz_xyzz[k] * ab_y + g_x_0_x_0_zzz_xyyzz[k];

                g_x_0_x_0_yzzz_xzzz[k] = -g_x_0_x_0_zzz_xzzz[k] * ab_y + g_x_0_x_0_zzz_xyzzz[k];

                g_x_0_x_0_yzzz_yyyy[k] = -g_x_0_x_0_zzz_yyyy[k] * ab_y + g_x_0_x_0_zzz_yyyyy[k];

                g_x_0_x_0_yzzz_yyyz[k] = -g_x_0_x_0_zzz_yyyz[k] * ab_y + g_x_0_x_0_zzz_yyyyz[k];

                g_x_0_x_0_yzzz_yyzz[k] = -g_x_0_x_0_zzz_yyzz[k] * ab_y + g_x_0_x_0_zzz_yyyzz[k];

                g_x_0_x_0_yzzz_yzzz[k] = -g_x_0_x_0_zzz_yzzz[k] * ab_y + g_x_0_x_0_zzz_yyzzz[k];

                g_x_0_x_0_yzzz_zzzz[k] = -g_x_0_x_0_zzz_zzzz[k] * ab_y + g_x_0_x_0_zzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_zzzz_xxxx = cbuffer.data(gg_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_xxxy = cbuffer.data(gg_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_xxxz = cbuffer.data(gg_geom_1010_off + 212 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_xxyy = cbuffer.data(gg_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_xxyz = cbuffer.data(gg_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_xxzz = cbuffer.data(gg_geom_1010_off + 215 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_xyyy = cbuffer.data(gg_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_xyyz = cbuffer.data(gg_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_xyzz = cbuffer.data(gg_geom_1010_off + 218 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_xzzz = cbuffer.data(gg_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_yyyy = cbuffer.data(gg_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_yyyz = cbuffer.data(gg_geom_1010_off + 221 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_yyzz = cbuffer.data(gg_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_yzzz = cbuffer.data(gg_geom_1010_off + 223 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_zzzz = cbuffer.data(gg_geom_1010_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_zzz_xxxx, g_x_0_x_0_zzz_xxxxz, g_x_0_x_0_zzz_xxxy, g_x_0_x_0_zzz_xxxyz, g_x_0_x_0_zzz_xxxz, g_x_0_x_0_zzz_xxxzz, g_x_0_x_0_zzz_xxyy, g_x_0_x_0_zzz_xxyyz, g_x_0_x_0_zzz_xxyz, g_x_0_x_0_zzz_xxyzz, g_x_0_x_0_zzz_xxzz, g_x_0_x_0_zzz_xxzzz, g_x_0_x_0_zzz_xyyy, g_x_0_x_0_zzz_xyyyz, g_x_0_x_0_zzz_xyyz, g_x_0_x_0_zzz_xyyzz, g_x_0_x_0_zzz_xyzz, g_x_0_x_0_zzz_xyzzz, g_x_0_x_0_zzz_xzzz, g_x_0_x_0_zzz_xzzzz, g_x_0_x_0_zzz_yyyy, g_x_0_x_0_zzz_yyyyz, g_x_0_x_0_zzz_yyyz, g_x_0_x_0_zzz_yyyzz, g_x_0_x_0_zzz_yyzz, g_x_0_x_0_zzz_yyzzz, g_x_0_x_0_zzz_yzzz, g_x_0_x_0_zzz_yzzzz, g_x_0_x_0_zzz_zzzz, g_x_0_x_0_zzz_zzzzz, g_x_0_x_0_zzzz_xxxx, g_x_0_x_0_zzzz_xxxy, g_x_0_x_0_zzzz_xxxz, g_x_0_x_0_zzzz_xxyy, g_x_0_x_0_zzzz_xxyz, g_x_0_x_0_zzzz_xxzz, g_x_0_x_0_zzzz_xyyy, g_x_0_x_0_zzzz_xyyz, g_x_0_x_0_zzzz_xyzz, g_x_0_x_0_zzzz_xzzz, g_x_0_x_0_zzzz_yyyy, g_x_0_x_0_zzzz_yyyz, g_x_0_x_0_zzzz_yyzz, g_x_0_x_0_zzzz_yzzz, g_x_0_x_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_zzzz_xxxx[k] = -g_x_0_x_0_zzz_xxxx[k] * ab_z + g_x_0_x_0_zzz_xxxxz[k];

                g_x_0_x_0_zzzz_xxxy[k] = -g_x_0_x_0_zzz_xxxy[k] * ab_z + g_x_0_x_0_zzz_xxxyz[k];

                g_x_0_x_0_zzzz_xxxz[k] = -g_x_0_x_0_zzz_xxxz[k] * ab_z + g_x_0_x_0_zzz_xxxzz[k];

                g_x_0_x_0_zzzz_xxyy[k] = -g_x_0_x_0_zzz_xxyy[k] * ab_z + g_x_0_x_0_zzz_xxyyz[k];

                g_x_0_x_0_zzzz_xxyz[k] = -g_x_0_x_0_zzz_xxyz[k] * ab_z + g_x_0_x_0_zzz_xxyzz[k];

                g_x_0_x_0_zzzz_xxzz[k] = -g_x_0_x_0_zzz_xxzz[k] * ab_z + g_x_0_x_0_zzz_xxzzz[k];

                g_x_0_x_0_zzzz_xyyy[k] = -g_x_0_x_0_zzz_xyyy[k] * ab_z + g_x_0_x_0_zzz_xyyyz[k];

                g_x_0_x_0_zzzz_xyyz[k] = -g_x_0_x_0_zzz_xyyz[k] * ab_z + g_x_0_x_0_zzz_xyyzz[k];

                g_x_0_x_0_zzzz_xyzz[k] = -g_x_0_x_0_zzz_xyzz[k] * ab_z + g_x_0_x_0_zzz_xyzzz[k];

                g_x_0_x_0_zzzz_xzzz[k] = -g_x_0_x_0_zzz_xzzz[k] * ab_z + g_x_0_x_0_zzz_xzzzz[k];

                g_x_0_x_0_zzzz_yyyy[k] = -g_x_0_x_0_zzz_yyyy[k] * ab_z + g_x_0_x_0_zzz_yyyyz[k];

                g_x_0_x_0_zzzz_yyyz[k] = -g_x_0_x_0_zzz_yyyz[k] * ab_z + g_x_0_x_0_zzz_yyyzz[k];

                g_x_0_x_0_zzzz_yyzz[k] = -g_x_0_x_0_zzz_yyzz[k] * ab_z + g_x_0_x_0_zzz_yyzzz[k];

                g_x_0_x_0_zzzz_yzzz[k] = -g_x_0_x_0_zzz_yzzz[k] * ab_z + g_x_0_x_0_zzz_yzzzz[k];

                g_x_0_x_0_zzzz_zzzz[k] = -g_x_0_x_0_zzz_zzzz[k] * ab_z + g_x_0_x_0_zzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxx_xxxx = cbuffer.data(gg_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_xxxy = cbuffer.data(gg_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_xxxz = cbuffer.data(gg_geom_1010_off + 227 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_xxyy = cbuffer.data(gg_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_xxyz = cbuffer.data(gg_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_xxzz = cbuffer.data(gg_geom_1010_off + 230 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_xyyy = cbuffer.data(gg_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_xyyz = cbuffer.data(gg_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_xyzz = cbuffer.data(gg_geom_1010_off + 233 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_xzzz = cbuffer.data(gg_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_yyyy = cbuffer.data(gg_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_yyyz = cbuffer.data(gg_geom_1010_off + 236 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_yyzz = cbuffer.data(gg_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_yzzz = cbuffer.data(gg_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_zzzz = cbuffer.data(gg_geom_1010_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_xxx_xxxx, g_0_0_y_0_xxx_xxxy, g_0_0_y_0_xxx_xxxz, g_0_0_y_0_xxx_xxyy, g_0_0_y_0_xxx_xxyz, g_0_0_y_0_xxx_xxzz, g_0_0_y_0_xxx_xyyy, g_0_0_y_0_xxx_xyyz, g_0_0_y_0_xxx_xyzz, g_0_0_y_0_xxx_xzzz, g_0_0_y_0_xxx_yyyy, g_0_0_y_0_xxx_yyyz, g_0_0_y_0_xxx_yyzz, g_0_0_y_0_xxx_yzzz, g_0_0_y_0_xxx_zzzz, g_x_0_y_0_xxx_xxxx, g_x_0_y_0_xxx_xxxxx, g_x_0_y_0_xxx_xxxxy, g_x_0_y_0_xxx_xxxxz, g_x_0_y_0_xxx_xxxy, g_x_0_y_0_xxx_xxxyy, g_x_0_y_0_xxx_xxxyz, g_x_0_y_0_xxx_xxxz, g_x_0_y_0_xxx_xxxzz, g_x_0_y_0_xxx_xxyy, g_x_0_y_0_xxx_xxyyy, g_x_0_y_0_xxx_xxyyz, g_x_0_y_0_xxx_xxyz, g_x_0_y_0_xxx_xxyzz, g_x_0_y_0_xxx_xxzz, g_x_0_y_0_xxx_xxzzz, g_x_0_y_0_xxx_xyyy, g_x_0_y_0_xxx_xyyyy, g_x_0_y_0_xxx_xyyyz, g_x_0_y_0_xxx_xyyz, g_x_0_y_0_xxx_xyyzz, g_x_0_y_0_xxx_xyzz, g_x_0_y_0_xxx_xyzzz, g_x_0_y_0_xxx_xzzz, g_x_0_y_0_xxx_xzzzz, g_x_0_y_0_xxx_yyyy, g_x_0_y_0_xxx_yyyz, g_x_0_y_0_xxx_yyzz, g_x_0_y_0_xxx_yzzz, g_x_0_y_0_xxx_zzzz, g_x_0_y_0_xxxx_xxxx, g_x_0_y_0_xxxx_xxxy, g_x_0_y_0_xxxx_xxxz, g_x_0_y_0_xxxx_xxyy, g_x_0_y_0_xxxx_xxyz, g_x_0_y_0_xxxx_xxzz, g_x_0_y_0_xxxx_xyyy, g_x_0_y_0_xxxx_xyyz, g_x_0_y_0_xxxx_xyzz, g_x_0_y_0_xxxx_xzzz, g_x_0_y_0_xxxx_yyyy, g_x_0_y_0_xxxx_yyyz, g_x_0_y_0_xxxx_yyzz, g_x_0_y_0_xxxx_yzzz, g_x_0_y_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxx_xxxx[k] = -g_0_0_y_0_xxx_xxxx[k] - g_x_0_y_0_xxx_xxxx[k] * ab_x + g_x_0_y_0_xxx_xxxxx[k];

                g_x_0_y_0_xxxx_xxxy[k] = -g_0_0_y_0_xxx_xxxy[k] - g_x_0_y_0_xxx_xxxy[k] * ab_x + g_x_0_y_0_xxx_xxxxy[k];

                g_x_0_y_0_xxxx_xxxz[k] = -g_0_0_y_0_xxx_xxxz[k] - g_x_0_y_0_xxx_xxxz[k] * ab_x + g_x_0_y_0_xxx_xxxxz[k];

                g_x_0_y_0_xxxx_xxyy[k] = -g_0_0_y_0_xxx_xxyy[k] - g_x_0_y_0_xxx_xxyy[k] * ab_x + g_x_0_y_0_xxx_xxxyy[k];

                g_x_0_y_0_xxxx_xxyz[k] = -g_0_0_y_0_xxx_xxyz[k] - g_x_0_y_0_xxx_xxyz[k] * ab_x + g_x_0_y_0_xxx_xxxyz[k];

                g_x_0_y_0_xxxx_xxzz[k] = -g_0_0_y_0_xxx_xxzz[k] - g_x_0_y_0_xxx_xxzz[k] * ab_x + g_x_0_y_0_xxx_xxxzz[k];

                g_x_0_y_0_xxxx_xyyy[k] = -g_0_0_y_0_xxx_xyyy[k] - g_x_0_y_0_xxx_xyyy[k] * ab_x + g_x_0_y_0_xxx_xxyyy[k];

                g_x_0_y_0_xxxx_xyyz[k] = -g_0_0_y_0_xxx_xyyz[k] - g_x_0_y_0_xxx_xyyz[k] * ab_x + g_x_0_y_0_xxx_xxyyz[k];

                g_x_0_y_0_xxxx_xyzz[k] = -g_0_0_y_0_xxx_xyzz[k] - g_x_0_y_0_xxx_xyzz[k] * ab_x + g_x_0_y_0_xxx_xxyzz[k];

                g_x_0_y_0_xxxx_xzzz[k] = -g_0_0_y_0_xxx_xzzz[k] - g_x_0_y_0_xxx_xzzz[k] * ab_x + g_x_0_y_0_xxx_xxzzz[k];

                g_x_0_y_0_xxxx_yyyy[k] = -g_0_0_y_0_xxx_yyyy[k] - g_x_0_y_0_xxx_yyyy[k] * ab_x + g_x_0_y_0_xxx_xyyyy[k];

                g_x_0_y_0_xxxx_yyyz[k] = -g_0_0_y_0_xxx_yyyz[k] - g_x_0_y_0_xxx_yyyz[k] * ab_x + g_x_0_y_0_xxx_xyyyz[k];

                g_x_0_y_0_xxxx_yyzz[k] = -g_0_0_y_0_xxx_yyzz[k] - g_x_0_y_0_xxx_yyzz[k] * ab_x + g_x_0_y_0_xxx_xyyzz[k];

                g_x_0_y_0_xxxx_yzzz[k] = -g_0_0_y_0_xxx_yzzz[k] - g_x_0_y_0_xxx_yzzz[k] * ab_x + g_x_0_y_0_xxx_xyzzz[k];

                g_x_0_y_0_xxxx_zzzz[k] = -g_0_0_y_0_xxx_zzzz[k] - g_x_0_y_0_xxx_zzzz[k] * ab_x + g_x_0_y_0_xxx_xzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxy_xxxx = cbuffer.data(gg_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_xxxy = cbuffer.data(gg_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_xxxz = cbuffer.data(gg_geom_1010_off + 242 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_xxyy = cbuffer.data(gg_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_xxyz = cbuffer.data(gg_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_xxzz = cbuffer.data(gg_geom_1010_off + 245 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_xyyy = cbuffer.data(gg_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_xyyz = cbuffer.data(gg_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_xyzz = cbuffer.data(gg_geom_1010_off + 248 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_xzzz = cbuffer.data(gg_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_yyyy = cbuffer.data(gg_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_yyyz = cbuffer.data(gg_geom_1010_off + 251 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_yyzz = cbuffer.data(gg_geom_1010_off + 252 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_yzzz = cbuffer.data(gg_geom_1010_off + 253 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_zzzz = cbuffer.data(gg_geom_1010_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxx_xxxx, g_x_0_y_0_xxx_xxxxy, g_x_0_y_0_xxx_xxxy, g_x_0_y_0_xxx_xxxyy, g_x_0_y_0_xxx_xxxyz, g_x_0_y_0_xxx_xxxz, g_x_0_y_0_xxx_xxyy, g_x_0_y_0_xxx_xxyyy, g_x_0_y_0_xxx_xxyyz, g_x_0_y_0_xxx_xxyz, g_x_0_y_0_xxx_xxyzz, g_x_0_y_0_xxx_xxzz, g_x_0_y_0_xxx_xyyy, g_x_0_y_0_xxx_xyyyy, g_x_0_y_0_xxx_xyyyz, g_x_0_y_0_xxx_xyyz, g_x_0_y_0_xxx_xyyzz, g_x_0_y_0_xxx_xyzz, g_x_0_y_0_xxx_xyzzz, g_x_0_y_0_xxx_xzzz, g_x_0_y_0_xxx_yyyy, g_x_0_y_0_xxx_yyyyy, g_x_0_y_0_xxx_yyyyz, g_x_0_y_0_xxx_yyyz, g_x_0_y_0_xxx_yyyzz, g_x_0_y_0_xxx_yyzz, g_x_0_y_0_xxx_yyzzz, g_x_0_y_0_xxx_yzzz, g_x_0_y_0_xxx_yzzzz, g_x_0_y_0_xxx_zzzz, g_x_0_y_0_xxxy_xxxx, g_x_0_y_0_xxxy_xxxy, g_x_0_y_0_xxxy_xxxz, g_x_0_y_0_xxxy_xxyy, g_x_0_y_0_xxxy_xxyz, g_x_0_y_0_xxxy_xxzz, g_x_0_y_0_xxxy_xyyy, g_x_0_y_0_xxxy_xyyz, g_x_0_y_0_xxxy_xyzz, g_x_0_y_0_xxxy_xzzz, g_x_0_y_0_xxxy_yyyy, g_x_0_y_0_xxxy_yyyz, g_x_0_y_0_xxxy_yyzz, g_x_0_y_0_xxxy_yzzz, g_x_0_y_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxy_xxxx[k] = -g_x_0_y_0_xxx_xxxx[k] * ab_y + g_x_0_y_0_xxx_xxxxy[k];

                g_x_0_y_0_xxxy_xxxy[k] = -g_x_0_y_0_xxx_xxxy[k] * ab_y + g_x_0_y_0_xxx_xxxyy[k];

                g_x_0_y_0_xxxy_xxxz[k] = -g_x_0_y_0_xxx_xxxz[k] * ab_y + g_x_0_y_0_xxx_xxxyz[k];

                g_x_0_y_0_xxxy_xxyy[k] = -g_x_0_y_0_xxx_xxyy[k] * ab_y + g_x_0_y_0_xxx_xxyyy[k];

                g_x_0_y_0_xxxy_xxyz[k] = -g_x_0_y_0_xxx_xxyz[k] * ab_y + g_x_0_y_0_xxx_xxyyz[k];

                g_x_0_y_0_xxxy_xxzz[k] = -g_x_0_y_0_xxx_xxzz[k] * ab_y + g_x_0_y_0_xxx_xxyzz[k];

                g_x_0_y_0_xxxy_xyyy[k] = -g_x_0_y_0_xxx_xyyy[k] * ab_y + g_x_0_y_0_xxx_xyyyy[k];

                g_x_0_y_0_xxxy_xyyz[k] = -g_x_0_y_0_xxx_xyyz[k] * ab_y + g_x_0_y_0_xxx_xyyyz[k];

                g_x_0_y_0_xxxy_xyzz[k] = -g_x_0_y_0_xxx_xyzz[k] * ab_y + g_x_0_y_0_xxx_xyyzz[k];

                g_x_0_y_0_xxxy_xzzz[k] = -g_x_0_y_0_xxx_xzzz[k] * ab_y + g_x_0_y_0_xxx_xyzzz[k];

                g_x_0_y_0_xxxy_yyyy[k] = -g_x_0_y_0_xxx_yyyy[k] * ab_y + g_x_0_y_0_xxx_yyyyy[k];

                g_x_0_y_0_xxxy_yyyz[k] = -g_x_0_y_0_xxx_yyyz[k] * ab_y + g_x_0_y_0_xxx_yyyyz[k];

                g_x_0_y_0_xxxy_yyzz[k] = -g_x_0_y_0_xxx_yyzz[k] * ab_y + g_x_0_y_0_xxx_yyyzz[k];

                g_x_0_y_0_xxxy_yzzz[k] = -g_x_0_y_0_xxx_yzzz[k] * ab_y + g_x_0_y_0_xxx_yyzzz[k];

                g_x_0_y_0_xxxy_zzzz[k] = -g_x_0_y_0_xxx_zzzz[k] * ab_y + g_x_0_y_0_xxx_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxz_xxxx = cbuffer.data(gg_geom_1010_off + 255 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_xxxy = cbuffer.data(gg_geom_1010_off + 256 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_xxxz = cbuffer.data(gg_geom_1010_off + 257 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_xxyy = cbuffer.data(gg_geom_1010_off + 258 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_xxyz = cbuffer.data(gg_geom_1010_off + 259 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_xxzz = cbuffer.data(gg_geom_1010_off + 260 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_xyyy = cbuffer.data(gg_geom_1010_off + 261 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_xyyz = cbuffer.data(gg_geom_1010_off + 262 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_xyzz = cbuffer.data(gg_geom_1010_off + 263 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_xzzz = cbuffer.data(gg_geom_1010_off + 264 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_yyyy = cbuffer.data(gg_geom_1010_off + 265 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_yyyz = cbuffer.data(gg_geom_1010_off + 266 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_yyzz = cbuffer.data(gg_geom_1010_off + 267 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_yzzz = cbuffer.data(gg_geom_1010_off + 268 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_zzzz = cbuffer.data(gg_geom_1010_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxx_xxxx, g_x_0_y_0_xxx_xxxxz, g_x_0_y_0_xxx_xxxy, g_x_0_y_0_xxx_xxxyz, g_x_0_y_0_xxx_xxxz, g_x_0_y_0_xxx_xxxzz, g_x_0_y_0_xxx_xxyy, g_x_0_y_0_xxx_xxyyz, g_x_0_y_0_xxx_xxyz, g_x_0_y_0_xxx_xxyzz, g_x_0_y_0_xxx_xxzz, g_x_0_y_0_xxx_xxzzz, g_x_0_y_0_xxx_xyyy, g_x_0_y_0_xxx_xyyyz, g_x_0_y_0_xxx_xyyz, g_x_0_y_0_xxx_xyyzz, g_x_0_y_0_xxx_xyzz, g_x_0_y_0_xxx_xyzzz, g_x_0_y_0_xxx_xzzz, g_x_0_y_0_xxx_xzzzz, g_x_0_y_0_xxx_yyyy, g_x_0_y_0_xxx_yyyyz, g_x_0_y_0_xxx_yyyz, g_x_0_y_0_xxx_yyyzz, g_x_0_y_0_xxx_yyzz, g_x_0_y_0_xxx_yyzzz, g_x_0_y_0_xxx_yzzz, g_x_0_y_0_xxx_yzzzz, g_x_0_y_0_xxx_zzzz, g_x_0_y_0_xxx_zzzzz, g_x_0_y_0_xxxz_xxxx, g_x_0_y_0_xxxz_xxxy, g_x_0_y_0_xxxz_xxxz, g_x_0_y_0_xxxz_xxyy, g_x_0_y_0_xxxz_xxyz, g_x_0_y_0_xxxz_xxzz, g_x_0_y_0_xxxz_xyyy, g_x_0_y_0_xxxz_xyyz, g_x_0_y_0_xxxz_xyzz, g_x_0_y_0_xxxz_xzzz, g_x_0_y_0_xxxz_yyyy, g_x_0_y_0_xxxz_yyyz, g_x_0_y_0_xxxz_yyzz, g_x_0_y_0_xxxz_yzzz, g_x_0_y_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxz_xxxx[k] = -g_x_0_y_0_xxx_xxxx[k] * ab_z + g_x_0_y_0_xxx_xxxxz[k];

                g_x_0_y_0_xxxz_xxxy[k] = -g_x_0_y_0_xxx_xxxy[k] * ab_z + g_x_0_y_0_xxx_xxxyz[k];

                g_x_0_y_0_xxxz_xxxz[k] = -g_x_0_y_0_xxx_xxxz[k] * ab_z + g_x_0_y_0_xxx_xxxzz[k];

                g_x_0_y_0_xxxz_xxyy[k] = -g_x_0_y_0_xxx_xxyy[k] * ab_z + g_x_0_y_0_xxx_xxyyz[k];

                g_x_0_y_0_xxxz_xxyz[k] = -g_x_0_y_0_xxx_xxyz[k] * ab_z + g_x_0_y_0_xxx_xxyzz[k];

                g_x_0_y_0_xxxz_xxzz[k] = -g_x_0_y_0_xxx_xxzz[k] * ab_z + g_x_0_y_0_xxx_xxzzz[k];

                g_x_0_y_0_xxxz_xyyy[k] = -g_x_0_y_0_xxx_xyyy[k] * ab_z + g_x_0_y_0_xxx_xyyyz[k];

                g_x_0_y_0_xxxz_xyyz[k] = -g_x_0_y_0_xxx_xyyz[k] * ab_z + g_x_0_y_0_xxx_xyyzz[k];

                g_x_0_y_0_xxxz_xyzz[k] = -g_x_0_y_0_xxx_xyzz[k] * ab_z + g_x_0_y_0_xxx_xyzzz[k];

                g_x_0_y_0_xxxz_xzzz[k] = -g_x_0_y_0_xxx_xzzz[k] * ab_z + g_x_0_y_0_xxx_xzzzz[k];

                g_x_0_y_0_xxxz_yyyy[k] = -g_x_0_y_0_xxx_yyyy[k] * ab_z + g_x_0_y_0_xxx_yyyyz[k];

                g_x_0_y_0_xxxz_yyyz[k] = -g_x_0_y_0_xxx_yyyz[k] * ab_z + g_x_0_y_0_xxx_yyyzz[k];

                g_x_0_y_0_xxxz_yyzz[k] = -g_x_0_y_0_xxx_yyzz[k] * ab_z + g_x_0_y_0_xxx_yyzzz[k];

                g_x_0_y_0_xxxz_yzzz[k] = -g_x_0_y_0_xxx_yzzz[k] * ab_z + g_x_0_y_0_xxx_yzzzz[k];

                g_x_0_y_0_xxxz_zzzz[k] = -g_x_0_y_0_xxx_zzzz[k] * ab_z + g_x_0_y_0_xxx_zzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyy_xxxx = cbuffer.data(gg_geom_1010_off + 270 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_xxxy = cbuffer.data(gg_geom_1010_off + 271 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_xxxz = cbuffer.data(gg_geom_1010_off + 272 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_xxyy = cbuffer.data(gg_geom_1010_off + 273 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_xxyz = cbuffer.data(gg_geom_1010_off + 274 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_xxzz = cbuffer.data(gg_geom_1010_off + 275 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_xyyy = cbuffer.data(gg_geom_1010_off + 276 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_xyyz = cbuffer.data(gg_geom_1010_off + 277 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_xyzz = cbuffer.data(gg_geom_1010_off + 278 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_xzzz = cbuffer.data(gg_geom_1010_off + 279 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_yyyy = cbuffer.data(gg_geom_1010_off + 280 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_yyyz = cbuffer.data(gg_geom_1010_off + 281 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_yyzz = cbuffer.data(gg_geom_1010_off + 282 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_yzzz = cbuffer.data(gg_geom_1010_off + 283 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_zzzz = cbuffer.data(gg_geom_1010_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxy_xxxx, g_x_0_y_0_xxy_xxxxy, g_x_0_y_0_xxy_xxxy, g_x_0_y_0_xxy_xxxyy, g_x_0_y_0_xxy_xxxyz, g_x_0_y_0_xxy_xxxz, g_x_0_y_0_xxy_xxyy, g_x_0_y_0_xxy_xxyyy, g_x_0_y_0_xxy_xxyyz, g_x_0_y_0_xxy_xxyz, g_x_0_y_0_xxy_xxyzz, g_x_0_y_0_xxy_xxzz, g_x_0_y_0_xxy_xyyy, g_x_0_y_0_xxy_xyyyy, g_x_0_y_0_xxy_xyyyz, g_x_0_y_0_xxy_xyyz, g_x_0_y_0_xxy_xyyzz, g_x_0_y_0_xxy_xyzz, g_x_0_y_0_xxy_xyzzz, g_x_0_y_0_xxy_xzzz, g_x_0_y_0_xxy_yyyy, g_x_0_y_0_xxy_yyyyy, g_x_0_y_0_xxy_yyyyz, g_x_0_y_0_xxy_yyyz, g_x_0_y_0_xxy_yyyzz, g_x_0_y_0_xxy_yyzz, g_x_0_y_0_xxy_yyzzz, g_x_0_y_0_xxy_yzzz, g_x_0_y_0_xxy_yzzzz, g_x_0_y_0_xxy_zzzz, g_x_0_y_0_xxyy_xxxx, g_x_0_y_0_xxyy_xxxy, g_x_0_y_0_xxyy_xxxz, g_x_0_y_0_xxyy_xxyy, g_x_0_y_0_xxyy_xxyz, g_x_0_y_0_xxyy_xxzz, g_x_0_y_0_xxyy_xyyy, g_x_0_y_0_xxyy_xyyz, g_x_0_y_0_xxyy_xyzz, g_x_0_y_0_xxyy_xzzz, g_x_0_y_0_xxyy_yyyy, g_x_0_y_0_xxyy_yyyz, g_x_0_y_0_xxyy_yyzz, g_x_0_y_0_xxyy_yzzz, g_x_0_y_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyy_xxxx[k] = -g_x_0_y_0_xxy_xxxx[k] * ab_y + g_x_0_y_0_xxy_xxxxy[k];

                g_x_0_y_0_xxyy_xxxy[k] = -g_x_0_y_0_xxy_xxxy[k] * ab_y + g_x_0_y_0_xxy_xxxyy[k];

                g_x_0_y_0_xxyy_xxxz[k] = -g_x_0_y_0_xxy_xxxz[k] * ab_y + g_x_0_y_0_xxy_xxxyz[k];

                g_x_0_y_0_xxyy_xxyy[k] = -g_x_0_y_0_xxy_xxyy[k] * ab_y + g_x_0_y_0_xxy_xxyyy[k];

                g_x_0_y_0_xxyy_xxyz[k] = -g_x_0_y_0_xxy_xxyz[k] * ab_y + g_x_0_y_0_xxy_xxyyz[k];

                g_x_0_y_0_xxyy_xxzz[k] = -g_x_0_y_0_xxy_xxzz[k] * ab_y + g_x_0_y_0_xxy_xxyzz[k];

                g_x_0_y_0_xxyy_xyyy[k] = -g_x_0_y_0_xxy_xyyy[k] * ab_y + g_x_0_y_0_xxy_xyyyy[k];

                g_x_0_y_0_xxyy_xyyz[k] = -g_x_0_y_0_xxy_xyyz[k] * ab_y + g_x_0_y_0_xxy_xyyyz[k];

                g_x_0_y_0_xxyy_xyzz[k] = -g_x_0_y_0_xxy_xyzz[k] * ab_y + g_x_0_y_0_xxy_xyyzz[k];

                g_x_0_y_0_xxyy_xzzz[k] = -g_x_0_y_0_xxy_xzzz[k] * ab_y + g_x_0_y_0_xxy_xyzzz[k];

                g_x_0_y_0_xxyy_yyyy[k] = -g_x_0_y_0_xxy_yyyy[k] * ab_y + g_x_0_y_0_xxy_yyyyy[k];

                g_x_0_y_0_xxyy_yyyz[k] = -g_x_0_y_0_xxy_yyyz[k] * ab_y + g_x_0_y_0_xxy_yyyyz[k];

                g_x_0_y_0_xxyy_yyzz[k] = -g_x_0_y_0_xxy_yyzz[k] * ab_y + g_x_0_y_0_xxy_yyyzz[k];

                g_x_0_y_0_xxyy_yzzz[k] = -g_x_0_y_0_xxy_yzzz[k] * ab_y + g_x_0_y_0_xxy_yyzzz[k];

                g_x_0_y_0_xxyy_zzzz[k] = -g_x_0_y_0_xxy_zzzz[k] * ab_y + g_x_0_y_0_xxy_yzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyz_xxxx = cbuffer.data(gg_geom_1010_off + 285 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_xxxy = cbuffer.data(gg_geom_1010_off + 286 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_xxxz = cbuffer.data(gg_geom_1010_off + 287 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_xxyy = cbuffer.data(gg_geom_1010_off + 288 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_xxyz = cbuffer.data(gg_geom_1010_off + 289 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_xxzz = cbuffer.data(gg_geom_1010_off + 290 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_xyyy = cbuffer.data(gg_geom_1010_off + 291 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_xyyz = cbuffer.data(gg_geom_1010_off + 292 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_xyzz = cbuffer.data(gg_geom_1010_off + 293 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_xzzz = cbuffer.data(gg_geom_1010_off + 294 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_yyyy = cbuffer.data(gg_geom_1010_off + 295 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_yyyz = cbuffer.data(gg_geom_1010_off + 296 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_yyzz = cbuffer.data(gg_geom_1010_off + 297 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_yzzz = cbuffer.data(gg_geom_1010_off + 298 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_zzzz = cbuffer.data(gg_geom_1010_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyz_xxxx, g_x_0_y_0_xxyz_xxxy, g_x_0_y_0_xxyz_xxxz, g_x_0_y_0_xxyz_xxyy, g_x_0_y_0_xxyz_xxyz, g_x_0_y_0_xxyz_xxzz, g_x_0_y_0_xxyz_xyyy, g_x_0_y_0_xxyz_xyyz, g_x_0_y_0_xxyz_xyzz, g_x_0_y_0_xxyz_xzzz, g_x_0_y_0_xxyz_yyyy, g_x_0_y_0_xxyz_yyyz, g_x_0_y_0_xxyz_yyzz, g_x_0_y_0_xxyz_yzzz, g_x_0_y_0_xxyz_zzzz, g_x_0_y_0_xxz_xxxx, g_x_0_y_0_xxz_xxxxy, g_x_0_y_0_xxz_xxxy, g_x_0_y_0_xxz_xxxyy, g_x_0_y_0_xxz_xxxyz, g_x_0_y_0_xxz_xxxz, g_x_0_y_0_xxz_xxyy, g_x_0_y_0_xxz_xxyyy, g_x_0_y_0_xxz_xxyyz, g_x_0_y_0_xxz_xxyz, g_x_0_y_0_xxz_xxyzz, g_x_0_y_0_xxz_xxzz, g_x_0_y_0_xxz_xyyy, g_x_0_y_0_xxz_xyyyy, g_x_0_y_0_xxz_xyyyz, g_x_0_y_0_xxz_xyyz, g_x_0_y_0_xxz_xyyzz, g_x_0_y_0_xxz_xyzz, g_x_0_y_0_xxz_xyzzz, g_x_0_y_0_xxz_xzzz, g_x_0_y_0_xxz_yyyy, g_x_0_y_0_xxz_yyyyy, g_x_0_y_0_xxz_yyyyz, g_x_0_y_0_xxz_yyyz, g_x_0_y_0_xxz_yyyzz, g_x_0_y_0_xxz_yyzz, g_x_0_y_0_xxz_yyzzz, g_x_0_y_0_xxz_yzzz, g_x_0_y_0_xxz_yzzzz, g_x_0_y_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyz_xxxx[k] = -g_x_0_y_0_xxz_xxxx[k] * ab_y + g_x_0_y_0_xxz_xxxxy[k];

                g_x_0_y_0_xxyz_xxxy[k] = -g_x_0_y_0_xxz_xxxy[k] * ab_y + g_x_0_y_0_xxz_xxxyy[k];

                g_x_0_y_0_xxyz_xxxz[k] = -g_x_0_y_0_xxz_xxxz[k] * ab_y + g_x_0_y_0_xxz_xxxyz[k];

                g_x_0_y_0_xxyz_xxyy[k] = -g_x_0_y_0_xxz_xxyy[k] * ab_y + g_x_0_y_0_xxz_xxyyy[k];

                g_x_0_y_0_xxyz_xxyz[k] = -g_x_0_y_0_xxz_xxyz[k] * ab_y + g_x_0_y_0_xxz_xxyyz[k];

                g_x_0_y_0_xxyz_xxzz[k] = -g_x_0_y_0_xxz_xxzz[k] * ab_y + g_x_0_y_0_xxz_xxyzz[k];

                g_x_0_y_0_xxyz_xyyy[k] = -g_x_0_y_0_xxz_xyyy[k] * ab_y + g_x_0_y_0_xxz_xyyyy[k];

                g_x_0_y_0_xxyz_xyyz[k] = -g_x_0_y_0_xxz_xyyz[k] * ab_y + g_x_0_y_0_xxz_xyyyz[k];

                g_x_0_y_0_xxyz_xyzz[k] = -g_x_0_y_0_xxz_xyzz[k] * ab_y + g_x_0_y_0_xxz_xyyzz[k];

                g_x_0_y_0_xxyz_xzzz[k] = -g_x_0_y_0_xxz_xzzz[k] * ab_y + g_x_0_y_0_xxz_xyzzz[k];

                g_x_0_y_0_xxyz_yyyy[k] = -g_x_0_y_0_xxz_yyyy[k] * ab_y + g_x_0_y_0_xxz_yyyyy[k];

                g_x_0_y_0_xxyz_yyyz[k] = -g_x_0_y_0_xxz_yyyz[k] * ab_y + g_x_0_y_0_xxz_yyyyz[k];

                g_x_0_y_0_xxyz_yyzz[k] = -g_x_0_y_0_xxz_yyzz[k] * ab_y + g_x_0_y_0_xxz_yyyzz[k];

                g_x_0_y_0_xxyz_yzzz[k] = -g_x_0_y_0_xxz_yzzz[k] * ab_y + g_x_0_y_0_xxz_yyzzz[k];

                g_x_0_y_0_xxyz_zzzz[k] = -g_x_0_y_0_xxz_zzzz[k] * ab_y + g_x_0_y_0_xxz_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxzz_xxxx = cbuffer.data(gg_geom_1010_off + 300 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_xxxy = cbuffer.data(gg_geom_1010_off + 301 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_xxxz = cbuffer.data(gg_geom_1010_off + 302 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_xxyy = cbuffer.data(gg_geom_1010_off + 303 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_xxyz = cbuffer.data(gg_geom_1010_off + 304 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_xxzz = cbuffer.data(gg_geom_1010_off + 305 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_xyyy = cbuffer.data(gg_geom_1010_off + 306 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_xyyz = cbuffer.data(gg_geom_1010_off + 307 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_xyzz = cbuffer.data(gg_geom_1010_off + 308 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_xzzz = cbuffer.data(gg_geom_1010_off + 309 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_yyyy = cbuffer.data(gg_geom_1010_off + 310 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_yyyz = cbuffer.data(gg_geom_1010_off + 311 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_yyzz = cbuffer.data(gg_geom_1010_off + 312 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_yzzz = cbuffer.data(gg_geom_1010_off + 313 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_zzzz = cbuffer.data(gg_geom_1010_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxz_xxxx, g_x_0_y_0_xxz_xxxxz, g_x_0_y_0_xxz_xxxy, g_x_0_y_0_xxz_xxxyz, g_x_0_y_0_xxz_xxxz, g_x_0_y_0_xxz_xxxzz, g_x_0_y_0_xxz_xxyy, g_x_0_y_0_xxz_xxyyz, g_x_0_y_0_xxz_xxyz, g_x_0_y_0_xxz_xxyzz, g_x_0_y_0_xxz_xxzz, g_x_0_y_0_xxz_xxzzz, g_x_0_y_0_xxz_xyyy, g_x_0_y_0_xxz_xyyyz, g_x_0_y_0_xxz_xyyz, g_x_0_y_0_xxz_xyyzz, g_x_0_y_0_xxz_xyzz, g_x_0_y_0_xxz_xyzzz, g_x_0_y_0_xxz_xzzz, g_x_0_y_0_xxz_xzzzz, g_x_0_y_0_xxz_yyyy, g_x_0_y_0_xxz_yyyyz, g_x_0_y_0_xxz_yyyz, g_x_0_y_0_xxz_yyyzz, g_x_0_y_0_xxz_yyzz, g_x_0_y_0_xxz_yyzzz, g_x_0_y_0_xxz_yzzz, g_x_0_y_0_xxz_yzzzz, g_x_0_y_0_xxz_zzzz, g_x_0_y_0_xxz_zzzzz, g_x_0_y_0_xxzz_xxxx, g_x_0_y_0_xxzz_xxxy, g_x_0_y_0_xxzz_xxxz, g_x_0_y_0_xxzz_xxyy, g_x_0_y_0_xxzz_xxyz, g_x_0_y_0_xxzz_xxzz, g_x_0_y_0_xxzz_xyyy, g_x_0_y_0_xxzz_xyyz, g_x_0_y_0_xxzz_xyzz, g_x_0_y_0_xxzz_xzzz, g_x_0_y_0_xxzz_yyyy, g_x_0_y_0_xxzz_yyyz, g_x_0_y_0_xxzz_yyzz, g_x_0_y_0_xxzz_yzzz, g_x_0_y_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxzz_xxxx[k] = -g_x_0_y_0_xxz_xxxx[k] * ab_z + g_x_0_y_0_xxz_xxxxz[k];

                g_x_0_y_0_xxzz_xxxy[k] = -g_x_0_y_0_xxz_xxxy[k] * ab_z + g_x_0_y_0_xxz_xxxyz[k];

                g_x_0_y_0_xxzz_xxxz[k] = -g_x_0_y_0_xxz_xxxz[k] * ab_z + g_x_0_y_0_xxz_xxxzz[k];

                g_x_0_y_0_xxzz_xxyy[k] = -g_x_0_y_0_xxz_xxyy[k] * ab_z + g_x_0_y_0_xxz_xxyyz[k];

                g_x_0_y_0_xxzz_xxyz[k] = -g_x_0_y_0_xxz_xxyz[k] * ab_z + g_x_0_y_0_xxz_xxyzz[k];

                g_x_0_y_0_xxzz_xxzz[k] = -g_x_0_y_0_xxz_xxzz[k] * ab_z + g_x_0_y_0_xxz_xxzzz[k];

                g_x_0_y_0_xxzz_xyyy[k] = -g_x_0_y_0_xxz_xyyy[k] * ab_z + g_x_0_y_0_xxz_xyyyz[k];

                g_x_0_y_0_xxzz_xyyz[k] = -g_x_0_y_0_xxz_xyyz[k] * ab_z + g_x_0_y_0_xxz_xyyzz[k];

                g_x_0_y_0_xxzz_xyzz[k] = -g_x_0_y_0_xxz_xyzz[k] * ab_z + g_x_0_y_0_xxz_xyzzz[k];

                g_x_0_y_0_xxzz_xzzz[k] = -g_x_0_y_0_xxz_xzzz[k] * ab_z + g_x_0_y_0_xxz_xzzzz[k];

                g_x_0_y_0_xxzz_yyyy[k] = -g_x_0_y_0_xxz_yyyy[k] * ab_z + g_x_0_y_0_xxz_yyyyz[k];

                g_x_0_y_0_xxzz_yyyz[k] = -g_x_0_y_0_xxz_yyyz[k] * ab_z + g_x_0_y_0_xxz_yyyzz[k];

                g_x_0_y_0_xxzz_yyzz[k] = -g_x_0_y_0_xxz_yyzz[k] * ab_z + g_x_0_y_0_xxz_yyzzz[k];

                g_x_0_y_0_xxzz_yzzz[k] = -g_x_0_y_0_xxz_yzzz[k] * ab_z + g_x_0_y_0_xxz_yzzzz[k];

                g_x_0_y_0_xxzz_zzzz[k] = -g_x_0_y_0_xxz_zzzz[k] * ab_z + g_x_0_y_0_xxz_zzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyy_xxxx = cbuffer.data(gg_geom_1010_off + 315 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_xxxy = cbuffer.data(gg_geom_1010_off + 316 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_xxxz = cbuffer.data(gg_geom_1010_off + 317 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_xxyy = cbuffer.data(gg_geom_1010_off + 318 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_xxyz = cbuffer.data(gg_geom_1010_off + 319 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_xxzz = cbuffer.data(gg_geom_1010_off + 320 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_xyyy = cbuffer.data(gg_geom_1010_off + 321 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_xyyz = cbuffer.data(gg_geom_1010_off + 322 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_xyzz = cbuffer.data(gg_geom_1010_off + 323 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_xzzz = cbuffer.data(gg_geom_1010_off + 324 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_yyyy = cbuffer.data(gg_geom_1010_off + 325 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_yyyz = cbuffer.data(gg_geom_1010_off + 326 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_yyzz = cbuffer.data(gg_geom_1010_off + 327 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_yzzz = cbuffer.data(gg_geom_1010_off + 328 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_zzzz = cbuffer.data(gg_geom_1010_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyy_xxxx, g_x_0_y_0_xyy_xxxxy, g_x_0_y_0_xyy_xxxy, g_x_0_y_0_xyy_xxxyy, g_x_0_y_0_xyy_xxxyz, g_x_0_y_0_xyy_xxxz, g_x_0_y_0_xyy_xxyy, g_x_0_y_0_xyy_xxyyy, g_x_0_y_0_xyy_xxyyz, g_x_0_y_0_xyy_xxyz, g_x_0_y_0_xyy_xxyzz, g_x_0_y_0_xyy_xxzz, g_x_0_y_0_xyy_xyyy, g_x_0_y_0_xyy_xyyyy, g_x_0_y_0_xyy_xyyyz, g_x_0_y_0_xyy_xyyz, g_x_0_y_0_xyy_xyyzz, g_x_0_y_0_xyy_xyzz, g_x_0_y_0_xyy_xyzzz, g_x_0_y_0_xyy_xzzz, g_x_0_y_0_xyy_yyyy, g_x_0_y_0_xyy_yyyyy, g_x_0_y_0_xyy_yyyyz, g_x_0_y_0_xyy_yyyz, g_x_0_y_0_xyy_yyyzz, g_x_0_y_0_xyy_yyzz, g_x_0_y_0_xyy_yyzzz, g_x_0_y_0_xyy_yzzz, g_x_0_y_0_xyy_yzzzz, g_x_0_y_0_xyy_zzzz, g_x_0_y_0_xyyy_xxxx, g_x_0_y_0_xyyy_xxxy, g_x_0_y_0_xyyy_xxxz, g_x_0_y_0_xyyy_xxyy, g_x_0_y_0_xyyy_xxyz, g_x_0_y_0_xyyy_xxzz, g_x_0_y_0_xyyy_xyyy, g_x_0_y_0_xyyy_xyyz, g_x_0_y_0_xyyy_xyzz, g_x_0_y_0_xyyy_xzzz, g_x_0_y_0_xyyy_yyyy, g_x_0_y_0_xyyy_yyyz, g_x_0_y_0_xyyy_yyzz, g_x_0_y_0_xyyy_yzzz, g_x_0_y_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyy_xxxx[k] = -g_x_0_y_0_xyy_xxxx[k] * ab_y + g_x_0_y_0_xyy_xxxxy[k];

                g_x_0_y_0_xyyy_xxxy[k] = -g_x_0_y_0_xyy_xxxy[k] * ab_y + g_x_0_y_0_xyy_xxxyy[k];

                g_x_0_y_0_xyyy_xxxz[k] = -g_x_0_y_0_xyy_xxxz[k] * ab_y + g_x_0_y_0_xyy_xxxyz[k];

                g_x_0_y_0_xyyy_xxyy[k] = -g_x_0_y_0_xyy_xxyy[k] * ab_y + g_x_0_y_0_xyy_xxyyy[k];

                g_x_0_y_0_xyyy_xxyz[k] = -g_x_0_y_0_xyy_xxyz[k] * ab_y + g_x_0_y_0_xyy_xxyyz[k];

                g_x_0_y_0_xyyy_xxzz[k] = -g_x_0_y_0_xyy_xxzz[k] * ab_y + g_x_0_y_0_xyy_xxyzz[k];

                g_x_0_y_0_xyyy_xyyy[k] = -g_x_0_y_0_xyy_xyyy[k] * ab_y + g_x_0_y_0_xyy_xyyyy[k];

                g_x_0_y_0_xyyy_xyyz[k] = -g_x_0_y_0_xyy_xyyz[k] * ab_y + g_x_0_y_0_xyy_xyyyz[k];

                g_x_0_y_0_xyyy_xyzz[k] = -g_x_0_y_0_xyy_xyzz[k] * ab_y + g_x_0_y_0_xyy_xyyzz[k];

                g_x_0_y_0_xyyy_xzzz[k] = -g_x_0_y_0_xyy_xzzz[k] * ab_y + g_x_0_y_0_xyy_xyzzz[k];

                g_x_0_y_0_xyyy_yyyy[k] = -g_x_0_y_0_xyy_yyyy[k] * ab_y + g_x_0_y_0_xyy_yyyyy[k];

                g_x_0_y_0_xyyy_yyyz[k] = -g_x_0_y_0_xyy_yyyz[k] * ab_y + g_x_0_y_0_xyy_yyyyz[k];

                g_x_0_y_0_xyyy_yyzz[k] = -g_x_0_y_0_xyy_yyzz[k] * ab_y + g_x_0_y_0_xyy_yyyzz[k];

                g_x_0_y_0_xyyy_yzzz[k] = -g_x_0_y_0_xyy_yzzz[k] * ab_y + g_x_0_y_0_xyy_yyzzz[k];

                g_x_0_y_0_xyyy_zzzz[k] = -g_x_0_y_0_xyy_zzzz[k] * ab_y + g_x_0_y_0_xyy_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyz_xxxx = cbuffer.data(gg_geom_1010_off + 330 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_xxxy = cbuffer.data(gg_geom_1010_off + 331 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_xxxz = cbuffer.data(gg_geom_1010_off + 332 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_xxyy = cbuffer.data(gg_geom_1010_off + 333 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_xxyz = cbuffer.data(gg_geom_1010_off + 334 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_xxzz = cbuffer.data(gg_geom_1010_off + 335 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_xyyy = cbuffer.data(gg_geom_1010_off + 336 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_xyyz = cbuffer.data(gg_geom_1010_off + 337 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_xyzz = cbuffer.data(gg_geom_1010_off + 338 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_xzzz = cbuffer.data(gg_geom_1010_off + 339 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_yyyy = cbuffer.data(gg_geom_1010_off + 340 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_yyyz = cbuffer.data(gg_geom_1010_off + 341 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_yyzz = cbuffer.data(gg_geom_1010_off + 342 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_yzzz = cbuffer.data(gg_geom_1010_off + 343 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_zzzz = cbuffer.data(gg_geom_1010_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyz_xxxx, g_x_0_y_0_xyyz_xxxy, g_x_0_y_0_xyyz_xxxz, g_x_0_y_0_xyyz_xxyy, g_x_0_y_0_xyyz_xxyz, g_x_0_y_0_xyyz_xxzz, g_x_0_y_0_xyyz_xyyy, g_x_0_y_0_xyyz_xyyz, g_x_0_y_0_xyyz_xyzz, g_x_0_y_0_xyyz_xzzz, g_x_0_y_0_xyyz_yyyy, g_x_0_y_0_xyyz_yyyz, g_x_0_y_0_xyyz_yyzz, g_x_0_y_0_xyyz_yzzz, g_x_0_y_0_xyyz_zzzz, g_x_0_y_0_xyz_xxxx, g_x_0_y_0_xyz_xxxxy, g_x_0_y_0_xyz_xxxy, g_x_0_y_0_xyz_xxxyy, g_x_0_y_0_xyz_xxxyz, g_x_0_y_0_xyz_xxxz, g_x_0_y_0_xyz_xxyy, g_x_0_y_0_xyz_xxyyy, g_x_0_y_0_xyz_xxyyz, g_x_0_y_0_xyz_xxyz, g_x_0_y_0_xyz_xxyzz, g_x_0_y_0_xyz_xxzz, g_x_0_y_0_xyz_xyyy, g_x_0_y_0_xyz_xyyyy, g_x_0_y_0_xyz_xyyyz, g_x_0_y_0_xyz_xyyz, g_x_0_y_0_xyz_xyyzz, g_x_0_y_0_xyz_xyzz, g_x_0_y_0_xyz_xyzzz, g_x_0_y_0_xyz_xzzz, g_x_0_y_0_xyz_yyyy, g_x_0_y_0_xyz_yyyyy, g_x_0_y_0_xyz_yyyyz, g_x_0_y_0_xyz_yyyz, g_x_0_y_0_xyz_yyyzz, g_x_0_y_0_xyz_yyzz, g_x_0_y_0_xyz_yyzzz, g_x_0_y_0_xyz_yzzz, g_x_0_y_0_xyz_yzzzz, g_x_0_y_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyz_xxxx[k] = -g_x_0_y_0_xyz_xxxx[k] * ab_y + g_x_0_y_0_xyz_xxxxy[k];

                g_x_0_y_0_xyyz_xxxy[k] = -g_x_0_y_0_xyz_xxxy[k] * ab_y + g_x_0_y_0_xyz_xxxyy[k];

                g_x_0_y_0_xyyz_xxxz[k] = -g_x_0_y_0_xyz_xxxz[k] * ab_y + g_x_0_y_0_xyz_xxxyz[k];

                g_x_0_y_0_xyyz_xxyy[k] = -g_x_0_y_0_xyz_xxyy[k] * ab_y + g_x_0_y_0_xyz_xxyyy[k];

                g_x_0_y_0_xyyz_xxyz[k] = -g_x_0_y_0_xyz_xxyz[k] * ab_y + g_x_0_y_0_xyz_xxyyz[k];

                g_x_0_y_0_xyyz_xxzz[k] = -g_x_0_y_0_xyz_xxzz[k] * ab_y + g_x_0_y_0_xyz_xxyzz[k];

                g_x_0_y_0_xyyz_xyyy[k] = -g_x_0_y_0_xyz_xyyy[k] * ab_y + g_x_0_y_0_xyz_xyyyy[k];

                g_x_0_y_0_xyyz_xyyz[k] = -g_x_0_y_0_xyz_xyyz[k] * ab_y + g_x_0_y_0_xyz_xyyyz[k];

                g_x_0_y_0_xyyz_xyzz[k] = -g_x_0_y_0_xyz_xyzz[k] * ab_y + g_x_0_y_0_xyz_xyyzz[k];

                g_x_0_y_0_xyyz_xzzz[k] = -g_x_0_y_0_xyz_xzzz[k] * ab_y + g_x_0_y_0_xyz_xyzzz[k];

                g_x_0_y_0_xyyz_yyyy[k] = -g_x_0_y_0_xyz_yyyy[k] * ab_y + g_x_0_y_0_xyz_yyyyy[k];

                g_x_0_y_0_xyyz_yyyz[k] = -g_x_0_y_0_xyz_yyyz[k] * ab_y + g_x_0_y_0_xyz_yyyyz[k];

                g_x_0_y_0_xyyz_yyzz[k] = -g_x_0_y_0_xyz_yyzz[k] * ab_y + g_x_0_y_0_xyz_yyyzz[k];

                g_x_0_y_0_xyyz_yzzz[k] = -g_x_0_y_0_xyz_yzzz[k] * ab_y + g_x_0_y_0_xyz_yyzzz[k];

                g_x_0_y_0_xyyz_zzzz[k] = -g_x_0_y_0_xyz_zzzz[k] * ab_y + g_x_0_y_0_xyz_yzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyzz_xxxx = cbuffer.data(gg_geom_1010_off + 345 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_xxxy = cbuffer.data(gg_geom_1010_off + 346 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_xxxz = cbuffer.data(gg_geom_1010_off + 347 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_xxyy = cbuffer.data(gg_geom_1010_off + 348 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_xxyz = cbuffer.data(gg_geom_1010_off + 349 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_xxzz = cbuffer.data(gg_geom_1010_off + 350 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_xyyy = cbuffer.data(gg_geom_1010_off + 351 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_xyyz = cbuffer.data(gg_geom_1010_off + 352 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_xyzz = cbuffer.data(gg_geom_1010_off + 353 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_xzzz = cbuffer.data(gg_geom_1010_off + 354 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_yyyy = cbuffer.data(gg_geom_1010_off + 355 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_yyyz = cbuffer.data(gg_geom_1010_off + 356 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_yyzz = cbuffer.data(gg_geom_1010_off + 357 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_yzzz = cbuffer.data(gg_geom_1010_off + 358 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_zzzz = cbuffer.data(gg_geom_1010_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyzz_xxxx, g_x_0_y_0_xyzz_xxxy, g_x_0_y_0_xyzz_xxxz, g_x_0_y_0_xyzz_xxyy, g_x_0_y_0_xyzz_xxyz, g_x_0_y_0_xyzz_xxzz, g_x_0_y_0_xyzz_xyyy, g_x_0_y_0_xyzz_xyyz, g_x_0_y_0_xyzz_xyzz, g_x_0_y_0_xyzz_xzzz, g_x_0_y_0_xyzz_yyyy, g_x_0_y_0_xyzz_yyyz, g_x_0_y_0_xyzz_yyzz, g_x_0_y_0_xyzz_yzzz, g_x_0_y_0_xyzz_zzzz, g_x_0_y_0_xzz_xxxx, g_x_0_y_0_xzz_xxxxy, g_x_0_y_0_xzz_xxxy, g_x_0_y_0_xzz_xxxyy, g_x_0_y_0_xzz_xxxyz, g_x_0_y_0_xzz_xxxz, g_x_0_y_0_xzz_xxyy, g_x_0_y_0_xzz_xxyyy, g_x_0_y_0_xzz_xxyyz, g_x_0_y_0_xzz_xxyz, g_x_0_y_0_xzz_xxyzz, g_x_0_y_0_xzz_xxzz, g_x_0_y_0_xzz_xyyy, g_x_0_y_0_xzz_xyyyy, g_x_0_y_0_xzz_xyyyz, g_x_0_y_0_xzz_xyyz, g_x_0_y_0_xzz_xyyzz, g_x_0_y_0_xzz_xyzz, g_x_0_y_0_xzz_xyzzz, g_x_0_y_0_xzz_xzzz, g_x_0_y_0_xzz_yyyy, g_x_0_y_0_xzz_yyyyy, g_x_0_y_0_xzz_yyyyz, g_x_0_y_0_xzz_yyyz, g_x_0_y_0_xzz_yyyzz, g_x_0_y_0_xzz_yyzz, g_x_0_y_0_xzz_yyzzz, g_x_0_y_0_xzz_yzzz, g_x_0_y_0_xzz_yzzzz, g_x_0_y_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyzz_xxxx[k] = -g_x_0_y_0_xzz_xxxx[k] * ab_y + g_x_0_y_0_xzz_xxxxy[k];

                g_x_0_y_0_xyzz_xxxy[k] = -g_x_0_y_0_xzz_xxxy[k] * ab_y + g_x_0_y_0_xzz_xxxyy[k];

                g_x_0_y_0_xyzz_xxxz[k] = -g_x_0_y_0_xzz_xxxz[k] * ab_y + g_x_0_y_0_xzz_xxxyz[k];

                g_x_0_y_0_xyzz_xxyy[k] = -g_x_0_y_0_xzz_xxyy[k] * ab_y + g_x_0_y_0_xzz_xxyyy[k];

                g_x_0_y_0_xyzz_xxyz[k] = -g_x_0_y_0_xzz_xxyz[k] * ab_y + g_x_0_y_0_xzz_xxyyz[k];

                g_x_0_y_0_xyzz_xxzz[k] = -g_x_0_y_0_xzz_xxzz[k] * ab_y + g_x_0_y_0_xzz_xxyzz[k];

                g_x_0_y_0_xyzz_xyyy[k] = -g_x_0_y_0_xzz_xyyy[k] * ab_y + g_x_0_y_0_xzz_xyyyy[k];

                g_x_0_y_0_xyzz_xyyz[k] = -g_x_0_y_0_xzz_xyyz[k] * ab_y + g_x_0_y_0_xzz_xyyyz[k];

                g_x_0_y_0_xyzz_xyzz[k] = -g_x_0_y_0_xzz_xyzz[k] * ab_y + g_x_0_y_0_xzz_xyyzz[k];

                g_x_0_y_0_xyzz_xzzz[k] = -g_x_0_y_0_xzz_xzzz[k] * ab_y + g_x_0_y_0_xzz_xyzzz[k];

                g_x_0_y_0_xyzz_yyyy[k] = -g_x_0_y_0_xzz_yyyy[k] * ab_y + g_x_0_y_0_xzz_yyyyy[k];

                g_x_0_y_0_xyzz_yyyz[k] = -g_x_0_y_0_xzz_yyyz[k] * ab_y + g_x_0_y_0_xzz_yyyyz[k];

                g_x_0_y_0_xyzz_yyzz[k] = -g_x_0_y_0_xzz_yyzz[k] * ab_y + g_x_0_y_0_xzz_yyyzz[k];

                g_x_0_y_0_xyzz_yzzz[k] = -g_x_0_y_0_xzz_yzzz[k] * ab_y + g_x_0_y_0_xzz_yyzzz[k];

                g_x_0_y_0_xyzz_zzzz[k] = -g_x_0_y_0_xzz_zzzz[k] * ab_y + g_x_0_y_0_xzz_yzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xzzz_xxxx = cbuffer.data(gg_geom_1010_off + 360 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_xxxy = cbuffer.data(gg_geom_1010_off + 361 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_xxxz = cbuffer.data(gg_geom_1010_off + 362 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_xxyy = cbuffer.data(gg_geom_1010_off + 363 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_xxyz = cbuffer.data(gg_geom_1010_off + 364 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_xxzz = cbuffer.data(gg_geom_1010_off + 365 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_xyyy = cbuffer.data(gg_geom_1010_off + 366 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_xyyz = cbuffer.data(gg_geom_1010_off + 367 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_xyzz = cbuffer.data(gg_geom_1010_off + 368 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_xzzz = cbuffer.data(gg_geom_1010_off + 369 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_yyyy = cbuffer.data(gg_geom_1010_off + 370 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_yyyz = cbuffer.data(gg_geom_1010_off + 371 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_yyzz = cbuffer.data(gg_geom_1010_off + 372 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_yzzz = cbuffer.data(gg_geom_1010_off + 373 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_zzzz = cbuffer.data(gg_geom_1010_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xzz_xxxx, g_x_0_y_0_xzz_xxxxz, g_x_0_y_0_xzz_xxxy, g_x_0_y_0_xzz_xxxyz, g_x_0_y_0_xzz_xxxz, g_x_0_y_0_xzz_xxxzz, g_x_0_y_0_xzz_xxyy, g_x_0_y_0_xzz_xxyyz, g_x_0_y_0_xzz_xxyz, g_x_0_y_0_xzz_xxyzz, g_x_0_y_0_xzz_xxzz, g_x_0_y_0_xzz_xxzzz, g_x_0_y_0_xzz_xyyy, g_x_0_y_0_xzz_xyyyz, g_x_0_y_0_xzz_xyyz, g_x_0_y_0_xzz_xyyzz, g_x_0_y_0_xzz_xyzz, g_x_0_y_0_xzz_xyzzz, g_x_0_y_0_xzz_xzzz, g_x_0_y_0_xzz_xzzzz, g_x_0_y_0_xzz_yyyy, g_x_0_y_0_xzz_yyyyz, g_x_0_y_0_xzz_yyyz, g_x_0_y_0_xzz_yyyzz, g_x_0_y_0_xzz_yyzz, g_x_0_y_0_xzz_yyzzz, g_x_0_y_0_xzz_yzzz, g_x_0_y_0_xzz_yzzzz, g_x_0_y_0_xzz_zzzz, g_x_0_y_0_xzz_zzzzz, g_x_0_y_0_xzzz_xxxx, g_x_0_y_0_xzzz_xxxy, g_x_0_y_0_xzzz_xxxz, g_x_0_y_0_xzzz_xxyy, g_x_0_y_0_xzzz_xxyz, g_x_0_y_0_xzzz_xxzz, g_x_0_y_0_xzzz_xyyy, g_x_0_y_0_xzzz_xyyz, g_x_0_y_0_xzzz_xyzz, g_x_0_y_0_xzzz_xzzz, g_x_0_y_0_xzzz_yyyy, g_x_0_y_0_xzzz_yyyz, g_x_0_y_0_xzzz_yyzz, g_x_0_y_0_xzzz_yzzz, g_x_0_y_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xzzz_xxxx[k] = -g_x_0_y_0_xzz_xxxx[k] * ab_z + g_x_0_y_0_xzz_xxxxz[k];

                g_x_0_y_0_xzzz_xxxy[k] = -g_x_0_y_0_xzz_xxxy[k] * ab_z + g_x_0_y_0_xzz_xxxyz[k];

                g_x_0_y_0_xzzz_xxxz[k] = -g_x_0_y_0_xzz_xxxz[k] * ab_z + g_x_0_y_0_xzz_xxxzz[k];

                g_x_0_y_0_xzzz_xxyy[k] = -g_x_0_y_0_xzz_xxyy[k] * ab_z + g_x_0_y_0_xzz_xxyyz[k];

                g_x_0_y_0_xzzz_xxyz[k] = -g_x_0_y_0_xzz_xxyz[k] * ab_z + g_x_0_y_0_xzz_xxyzz[k];

                g_x_0_y_0_xzzz_xxzz[k] = -g_x_0_y_0_xzz_xxzz[k] * ab_z + g_x_0_y_0_xzz_xxzzz[k];

                g_x_0_y_0_xzzz_xyyy[k] = -g_x_0_y_0_xzz_xyyy[k] * ab_z + g_x_0_y_0_xzz_xyyyz[k];

                g_x_0_y_0_xzzz_xyyz[k] = -g_x_0_y_0_xzz_xyyz[k] * ab_z + g_x_0_y_0_xzz_xyyzz[k];

                g_x_0_y_0_xzzz_xyzz[k] = -g_x_0_y_0_xzz_xyzz[k] * ab_z + g_x_0_y_0_xzz_xyzzz[k];

                g_x_0_y_0_xzzz_xzzz[k] = -g_x_0_y_0_xzz_xzzz[k] * ab_z + g_x_0_y_0_xzz_xzzzz[k];

                g_x_0_y_0_xzzz_yyyy[k] = -g_x_0_y_0_xzz_yyyy[k] * ab_z + g_x_0_y_0_xzz_yyyyz[k];

                g_x_0_y_0_xzzz_yyyz[k] = -g_x_0_y_0_xzz_yyyz[k] * ab_z + g_x_0_y_0_xzz_yyyzz[k];

                g_x_0_y_0_xzzz_yyzz[k] = -g_x_0_y_0_xzz_yyzz[k] * ab_z + g_x_0_y_0_xzz_yyzzz[k];

                g_x_0_y_0_xzzz_yzzz[k] = -g_x_0_y_0_xzz_yzzz[k] * ab_z + g_x_0_y_0_xzz_yzzzz[k];

                g_x_0_y_0_xzzz_zzzz[k] = -g_x_0_y_0_xzz_zzzz[k] * ab_z + g_x_0_y_0_xzz_zzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyy_xxxx = cbuffer.data(gg_geom_1010_off + 375 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_xxxy = cbuffer.data(gg_geom_1010_off + 376 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_xxxz = cbuffer.data(gg_geom_1010_off + 377 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_xxyy = cbuffer.data(gg_geom_1010_off + 378 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_xxyz = cbuffer.data(gg_geom_1010_off + 379 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_xxzz = cbuffer.data(gg_geom_1010_off + 380 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_xyyy = cbuffer.data(gg_geom_1010_off + 381 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_xyyz = cbuffer.data(gg_geom_1010_off + 382 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_xyzz = cbuffer.data(gg_geom_1010_off + 383 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_xzzz = cbuffer.data(gg_geom_1010_off + 384 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_yyyy = cbuffer.data(gg_geom_1010_off + 385 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_yyyz = cbuffer.data(gg_geom_1010_off + 386 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_yyzz = cbuffer.data(gg_geom_1010_off + 387 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_yzzz = cbuffer.data(gg_geom_1010_off + 388 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_zzzz = cbuffer.data(gg_geom_1010_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyy_xxxx, g_x_0_y_0_yyy_xxxxy, g_x_0_y_0_yyy_xxxy, g_x_0_y_0_yyy_xxxyy, g_x_0_y_0_yyy_xxxyz, g_x_0_y_0_yyy_xxxz, g_x_0_y_0_yyy_xxyy, g_x_0_y_0_yyy_xxyyy, g_x_0_y_0_yyy_xxyyz, g_x_0_y_0_yyy_xxyz, g_x_0_y_0_yyy_xxyzz, g_x_0_y_0_yyy_xxzz, g_x_0_y_0_yyy_xyyy, g_x_0_y_0_yyy_xyyyy, g_x_0_y_0_yyy_xyyyz, g_x_0_y_0_yyy_xyyz, g_x_0_y_0_yyy_xyyzz, g_x_0_y_0_yyy_xyzz, g_x_0_y_0_yyy_xyzzz, g_x_0_y_0_yyy_xzzz, g_x_0_y_0_yyy_yyyy, g_x_0_y_0_yyy_yyyyy, g_x_0_y_0_yyy_yyyyz, g_x_0_y_0_yyy_yyyz, g_x_0_y_0_yyy_yyyzz, g_x_0_y_0_yyy_yyzz, g_x_0_y_0_yyy_yyzzz, g_x_0_y_0_yyy_yzzz, g_x_0_y_0_yyy_yzzzz, g_x_0_y_0_yyy_zzzz, g_x_0_y_0_yyyy_xxxx, g_x_0_y_0_yyyy_xxxy, g_x_0_y_0_yyyy_xxxz, g_x_0_y_0_yyyy_xxyy, g_x_0_y_0_yyyy_xxyz, g_x_0_y_0_yyyy_xxzz, g_x_0_y_0_yyyy_xyyy, g_x_0_y_0_yyyy_xyyz, g_x_0_y_0_yyyy_xyzz, g_x_0_y_0_yyyy_xzzz, g_x_0_y_0_yyyy_yyyy, g_x_0_y_0_yyyy_yyyz, g_x_0_y_0_yyyy_yyzz, g_x_0_y_0_yyyy_yzzz, g_x_0_y_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyy_xxxx[k] = -g_x_0_y_0_yyy_xxxx[k] * ab_y + g_x_0_y_0_yyy_xxxxy[k];

                g_x_0_y_0_yyyy_xxxy[k] = -g_x_0_y_0_yyy_xxxy[k] * ab_y + g_x_0_y_0_yyy_xxxyy[k];

                g_x_0_y_0_yyyy_xxxz[k] = -g_x_0_y_0_yyy_xxxz[k] * ab_y + g_x_0_y_0_yyy_xxxyz[k];

                g_x_0_y_0_yyyy_xxyy[k] = -g_x_0_y_0_yyy_xxyy[k] * ab_y + g_x_0_y_0_yyy_xxyyy[k];

                g_x_0_y_0_yyyy_xxyz[k] = -g_x_0_y_0_yyy_xxyz[k] * ab_y + g_x_0_y_0_yyy_xxyyz[k];

                g_x_0_y_0_yyyy_xxzz[k] = -g_x_0_y_0_yyy_xxzz[k] * ab_y + g_x_0_y_0_yyy_xxyzz[k];

                g_x_0_y_0_yyyy_xyyy[k] = -g_x_0_y_0_yyy_xyyy[k] * ab_y + g_x_0_y_0_yyy_xyyyy[k];

                g_x_0_y_0_yyyy_xyyz[k] = -g_x_0_y_0_yyy_xyyz[k] * ab_y + g_x_0_y_0_yyy_xyyyz[k];

                g_x_0_y_0_yyyy_xyzz[k] = -g_x_0_y_0_yyy_xyzz[k] * ab_y + g_x_0_y_0_yyy_xyyzz[k];

                g_x_0_y_0_yyyy_xzzz[k] = -g_x_0_y_0_yyy_xzzz[k] * ab_y + g_x_0_y_0_yyy_xyzzz[k];

                g_x_0_y_0_yyyy_yyyy[k] = -g_x_0_y_0_yyy_yyyy[k] * ab_y + g_x_0_y_0_yyy_yyyyy[k];

                g_x_0_y_0_yyyy_yyyz[k] = -g_x_0_y_0_yyy_yyyz[k] * ab_y + g_x_0_y_0_yyy_yyyyz[k];

                g_x_0_y_0_yyyy_yyzz[k] = -g_x_0_y_0_yyy_yyzz[k] * ab_y + g_x_0_y_0_yyy_yyyzz[k];

                g_x_0_y_0_yyyy_yzzz[k] = -g_x_0_y_0_yyy_yzzz[k] * ab_y + g_x_0_y_0_yyy_yyzzz[k];

                g_x_0_y_0_yyyy_zzzz[k] = -g_x_0_y_0_yyy_zzzz[k] * ab_y + g_x_0_y_0_yyy_yzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyz_xxxx = cbuffer.data(gg_geom_1010_off + 390 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_xxxy = cbuffer.data(gg_geom_1010_off + 391 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_xxxz = cbuffer.data(gg_geom_1010_off + 392 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_xxyy = cbuffer.data(gg_geom_1010_off + 393 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_xxyz = cbuffer.data(gg_geom_1010_off + 394 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_xxzz = cbuffer.data(gg_geom_1010_off + 395 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_xyyy = cbuffer.data(gg_geom_1010_off + 396 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_xyyz = cbuffer.data(gg_geom_1010_off + 397 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_xyzz = cbuffer.data(gg_geom_1010_off + 398 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_xzzz = cbuffer.data(gg_geom_1010_off + 399 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_yyyy = cbuffer.data(gg_geom_1010_off + 400 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_yyyz = cbuffer.data(gg_geom_1010_off + 401 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_yyzz = cbuffer.data(gg_geom_1010_off + 402 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_yzzz = cbuffer.data(gg_geom_1010_off + 403 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_zzzz = cbuffer.data(gg_geom_1010_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyz_xxxx, g_x_0_y_0_yyyz_xxxy, g_x_0_y_0_yyyz_xxxz, g_x_0_y_0_yyyz_xxyy, g_x_0_y_0_yyyz_xxyz, g_x_0_y_0_yyyz_xxzz, g_x_0_y_0_yyyz_xyyy, g_x_0_y_0_yyyz_xyyz, g_x_0_y_0_yyyz_xyzz, g_x_0_y_0_yyyz_xzzz, g_x_0_y_0_yyyz_yyyy, g_x_0_y_0_yyyz_yyyz, g_x_0_y_0_yyyz_yyzz, g_x_0_y_0_yyyz_yzzz, g_x_0_y_0_yyyz_zzzz, g_x_0_y_0_yyz_xxxx, g_x_0_y_0_yyz_xxxxy, g_x_0_y_0_yyz_xxxy, g_x_0_y_0_yyz_xxxyy, g_x_0_y_0_yyz_xxxyz, g_x_0_y_0_yyz_xxxz, g_x_0_y_0_yyz_xxyy, g_x_0_y_0_yyz_xxyyy, g_x_0_y_0_yyz_xxyyz, g_x_0_y_0_yyz_xxyz, g_x_0_y_0_yyz_xxyzz, g_x_0_y_0_yyz_xxzz, g_x_0_y_0_yyz_xyyy, g_x_0_y_0_yyz_xyyyy, g_x_0_y_0_yyz_xyyyz, g_x_0_y_0_yyz_xyyz, g_x_0_y_0_yyz_xyyzz, g_x_0_y_0_yyz_xyzz, g_x_0_y_0_yyz_xyzzz, g_x_0_y_0_yyz_xzzz, g_x_0_y_0_yyz_yyyy, g_x_0_y_0_yyz_yyyyy, g_x_0_y_0_yyz_yyyyz, g_x_0_y_0_yyz_yyyz, g_x_0_y_0_yyz_yyyzz, g_x_0_y_0_yyz_yyzz, g_x_0_y_0_yyz_yyzzz, g_x_0_y_0_yyz_yzzz, g_x_0_y_0_yyz_yzzzz, g_x_0_y_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyz_xxxx[k] = -g_x_0_y_0_yyz_xxxx[k] * ab_y + g_x_0_y_0_yyz_xxxxy[k];

                g_x_0_y_0_yyyz_xxxy[k] = -g_x_0_y_0_yyz_xxxy[k] * ab_y + g_x_0_y_0_yyz_xxxyy[k];

                g_x_0_y_0_yyyz_xxxz[k] = -g_x_0_y_0_yyz_xxxz[k] * ab_y + g_x_0_y_0_yyz_xxxyz[k];

                g_x_0_y_0_yyyz_xxyy[k] = -g_x_0_y_0_yyz_xxyy[k] * ab_y + g_x_0_y_0_yyz_xxyyy[k];

                g_x_0_y_0_yyyz_xxyz[k] = -g_x_0_y_0_yyz_xxyz[k] * ab_y + g_x_0_y_0_yyz_xxyyz[k];

                g_x_0_y_0_yyyz_xxzz[k] = -g_x_0_y_0_yyz_xxzz[k] * ab_y + g_x_0_y_0_yyz_xxyzz[k];

                g_x_0_y_0_yyyz_xyyy[k] = -g_x_0_y_0_yyz_xyyy[k] * ab_y + g_x_0_y_0_yyz_xyyyy[k];

                g_x_0_y_0_yyyz_xyyz[k] = -g_x_0_y_0_yyz_xyyz[k] * ab_y + g_x_0_y_0_yyz_xyyyz[k];

                g_x_0_y_0_yyyz_xyzz[k] = -g_x_0_y_0_yyz_xyzz[k] * ab_y + g_x_0_y_0_yyz_xyyzz[k];

                g_x_0_y_0_yyyz_xzzz[k] = -g_x_0_y_0_yyz_xzzz[k] * ab_y + g_x_0_y_0_yyz_xyzzz[k];

                g_x_0_y_0_yyyz_yyyy[k] = -g_x_0_y_0_yyz_yyyy[k] * ab_y + g_x_0_y_0_yyz_yyyyy[k];

                g_x_0_y_0_yyyz_yyyz[k] = -g_x_0_y_0_yyz_yyyz[k] * ab_y + g_x_0_y_0_yyz_yyyyz[k];

                g_x_0_y_0_yyyz_yyzz[k] = -g_x_0_y_0_yyz_yyzz[k] * ab_y + g_x_0_y_0_yyz_yyyzz[k];

                g_x_0_y_0_yyyz_yzzz[k] = -g_x_0_y_0_yyz_yzzz[k] * ab_y + g_x_0_y_0_yyz_yyzzz[k];

                g_x_0_y_0_yyyz_zzzz[k] = -g_x_0_y_0_yyz_zzzz[k] * ab_y + g_x_0_y_0_yyz_yzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyzz_xxxx = cbuffer.data(gg_geom_1010_off + 405 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_xxxy = cbuffer.data(gg_geom_1010_off + 406 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_xxxz = cbuffer.data(gg_geom_1010_off + 407 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_xxyy = cbuffer.data(gg_geom_1010_off + 408 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_xxyz = cbuffer.data(gg_geom_1010_off + 409 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_xxzz = cbuffer.data(gg_geom_1010_off + 410 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_xyyy = cbuffer.data(gg_geom_1010_off + 411 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_xyyz = cbuffer.data(gg_geom_1010_off + 412 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_xyzz = cbuffer.data(gg_geom_1010_off + 413 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_xzzz = cbuffer.data(gg_geom_1010_off + 414 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_yyyy = cbuffer.data(gg_geom_1010_off + 415 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_yyyz = cbuffer.data(gg_geom_1010_off + 416 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_yyzz = cbuffer.data(gg_geom_1010_off + 417 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_yzzz = cbuffer.data(gg_geom_1010_off + 418 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_zzzz = cbuffer.data(gg_geom_1010_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyzz_xxxx, g_x_0_y_0_yyzz_xxxy, g_x_0_y_0_yyzz_xxxz, g_x_0_y_0_yyzz_xxyy, g_x_0_y_0_yyzz_xxyz, g_x_0_y_0_yyzz_xxzz, g_x_0_y_0_yyzz_xyyy, g_x_0_y_0_yyzz_xyyz, g_x_0_y_0_yyzz_xyzz, g_x_0_y_0_yyzz_xzzz, g_x_0_y_0_yyzz_yyyy, g_x_0_y_0_yyzz_yyyz, g_x_0_y_0_yyzz_yyzz, g_x_0_y_0_yyzz_yzzz, g_x_0_y_0_yyzz_zzzz, g_x_0_y_0_yzz_xxxx, g_x_0_y_0_yzz_xxxxy, g_x_0_y_0_yzz_xxxy, g_x_0_y_0_yzz_xxxyy, g_x_0_y_0_yzz_xxxyz, g_x_0_y_0_yzz_xxxz, g_x_0_y_0_yzz_xxyy, g_x_0_y_0_yzz_xxyyy, g_x_0_y_0_yzz_xxyyz, g_x_0_y_0_yzz_xxyz, g_x_0_y_0_yzz_xxyzz, g_x_0_y_0_yzz_xxzz, g_x_0_y_0_yzz_xyyy, g_x_0_y_0_yzz_xyyyy, g_x_0_y_0_yzz_xyyyz, g_x_0_y_0_yzz_xyyz, g_x_0_y_0_yzz_xyyzz, g_x_0_y_0_yzz_xyzz, g_x_0_y_0_yzz_xyzzz, g_x_0_y_0_yzz_xzzz, g_x_0_y_0_yzz_yyyy, g_x_0_y_0_yzz_yyyyy, g_x_0_y_0_yzz_yyyyz, g_x_0_y_0_yzz_yyyz, g_x_0_y_0_yzz_yyyzz, g_x_0_y_0_yzz_yyzz, g_x_0_y_0_yzz_yyzzz, g_x_0_y_0_yzz_yzzz, g_x_0_y_0_yzz_yzzzz, g_x_0_y_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyzz_xxxx[k] = -g_x_0_y_0_yzz_xxxx[k] * ab_y + g_x_0_y_0_yzz_xxxxy[k];

                g_x_0_y_0_yyzz_xxxy[k] = -g_x_0_y_0_yzz_xxxy[k] * ab_y + g_x_0_y_0_yzz_xxxyy[k];

                g_x_0_y_0_yyzz_xxxz[k] = -g_x_0_y_0_yzz_xxxz[k] * ab_y + g_x_0_y_0_yzz_xxxyz[k];

                g_x_0_y_0_yyzz_xxyy[k] = -g_x_0_y_0_yzz_xxyy[k] * ab_y + g_x_0_y_0_yzz_xxyyy[k];

                g_x_0_y_0_yyzz_xxyz[k] = -g_x_0_y_0_yzz_xxyz[k] * ab_y + g_x_0_y_0_yzz_xxyyz[k];

                g_x_0_y_0_yyzz_xxzz[k] = -g_x_0_y_0_yzz_xxzz[k] * ab_y + g_x_0_y_0_yzz_xxyzz[k];

                g_x_0_y_0_yyzz_xyyy[k] = -g_x_0_y_0_yzz_xyyy[k] * ab_y + g_x_0_y_0_yzz_xyyyy[k];

                g_x_0_y_0_yyzz_xyyz[k] = -g_x_0_y_0_yzz_xyyz[k] * ab_y + g_x_0_y_0_yzz_xyyyz[k];

                g_x_0_y_0_yyzz_xyzz[k] = -g_x_0_y_0_yzz_xyzz[k] * ab_y + g_x_0_y_0_yzz_xyyzz[k];

                g_x_0_y_0_yyzz_xzzz[k] = -g_x_0_y_0_yzz_xzzz[k] * ab_y + g_x_0_y_0_yzz_xyzzz[k];

                g_x_0_y_0_yyzz_yyyy[k] = -g_x_0_y_0_yzz_yyyy[k] * ab_y + g_x_0_y_0_yzz_yyyyy[k];

                g_x_0_y_0_yyzz_yyyz[k] = -g_x_0_y_0_yzz_yyyz[k] * ab_y + g_x_0_y_0_yzz_yyyyz[k];

                g_x_0_y_0_yyzz_yyzz[k] = -g_x_0_y_0_yzz_yyzz[k] * ab_y + g_x_0_y_0_yzz_yyyzz[k];

                g_x_0_y_0_yyzz_yzzz[k] = -g_x_0_y_0_yzz_yzzz[k] * ab_y + g_x_0_y_0_yzz_yyzzz[k];

                g_x_0_y_0_yyzz_zzzz[k] = -g_x_0_y_0_yzz_zzzz[k] * ab_y + g_x_0_y_0_yzz_yzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yzzz_xxxx = cbuffer.data(gg_geom_1010_off + 420 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_xxxy = cbuffer.data(gg_geom_1010_off + 421 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_xxxz = cbuffer.data(gg_geom_1010_off + 422 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_xxyy = cbuffer.data(gg_geom_1010_off + 423 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_xxyz = cbuffer.data(gg_geom_1010_off + 424 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_xxzz = cbuffer.data(gg_geom_1010_off + 425 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_xyyy = cbuffer.data(gg_geom_1010_off + 426 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_xyyz = cbuffer.data(gg_geom_1010_off + 427 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_xyzz = cbuffer.data(gg_geom_1010_off + 428 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_xzzz = cbuffer.data(gg_geom_1010_off + 429 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_yyyy = cbuffer.data(gg_geom_1010_off + 430 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_yyyz = cbuffer.data(gg_geom_1010_off + 431 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_yyzz = cbuffer.data(gg_geom_1010_off + 432 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_yzzz = cbuffer.data(gg_geom_1010_off + 433 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_zzzz = cbuffer.data(gg_geom_1010_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yzzz_xxxx, g_x_0_y_0_yzzz_xxxy, g_x_0_y_0_yzzz_xxxz, g_x_0_y_0_yzzz_xxyy, g_x_0_y_0_yzzz_xxyz, g_x_0_y_0_yzzz_xxzz, g_x_0_y_0_yzzz_xyyy, g_x_0_y_0_yzzz_xyyz, g_x_0_y_0_yzzz_xyzz, g_x_0_y_0_yzzz_xzzz, g_x_0_y_0_yzzz_yyyy, g_x_0_y_0_yzzz_yyyz, g_x_0_y_0_yzzz_yyzz, g_x_0_y_0_yzzz_yzzz, g_x_0_y_0_yzzz_zzzz, g_x_0_y_0_zzz_xxxx, g_x_0_y_0_zzz_xxxxy, g_x_0_y_0_zzz_xxxy, g_x_0_y_0_zzz_xxxyy, g_x_0_y_0_zzz_xxxyz, g_x_0_y_0_zzz_xxxz, g_x_0_y_0_zzz_xxyy, g_x_0_y_0_zzz_xxyyy, g_x_0_y_0_zzz_xxyyz, g_x_0_y_0_zzz_xxyz, g_x_0_y_0_zzz_xxyzz, g_x_0_y_0_zzz_xxzz, g_x_0_y_0_zzz_xyyy, g_x_0_y_0_zzz_xyyyy, g_x_0_y_0_zzz_xyyyz, g_x_0_y_0_zzz_xyyz, g_x_0_y_0_zzz_xyyzz, g_x_0_y_0_zzz_xyzz, g_x_0_y_0_zzz_xyzzz, g_x_0_y_0_zzz_xzzz, g_x_0_y_0_zzz_yyyy, g_x_0_y_0_zzz_yyyyy, g_x_0_y_0_zzz_yyyyz, g_x_0_y_0_zzz_yyyz, g_x_0_y_0_zzz_yyyzz, g_x_0_y_0_zzz_yyzz, g_x_0_y_0_zzz_yyzzz, g_x_0_y_0_zzz_yzzz, g_x_0_y_0_zzz_yzzzz, g_x_0_y_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yzzz_xxxx[k] = -g_x_0_y_0_zzz_xxxx[k] * ab_y + g_x_0_y_0_zzz_xxxxy[k];

                g_x_0_y_0_yzzz_xxxy[k] = -g_x_0_y_0_zzz_xxxy[k] * ab_y + g_x_0_y_0_zzz_xxxyy[k];

                g_x_0_y_0_yzzz_xxxz[k] = -g_x_0_y_0_zzz_xxxz[k] * ab_y + g_x_0_y_0_zzz_xxxyz[k];

                g_x_0_y_0_yzzz_xxyy[k] = -g_x_0_y_0_zzz_xxyy[k] * ab_y + g_x_0_y_0_zzz_xxyyy[k];

                g_x_0_y_0_yzzz_xxyz[k] = -g_x_0_y_0_zzz_xxyz[k] * ab_y + g_x_0_y_0_zzz_xxyyz[k];

                g_x_0_y_0_yzzz_xxzz[k] = -g_x_0_y_0_zzz_xxzz[k] * ab_y + g_x_0_y_0_zzz_xxyzz[k];

                g_x_0_y_0_yzzz_xyyy[k] = -g_x_0_y_0_zzz_xyyy[k] * ab_y + g_x_0_y_0_zzz_xyyyy[k];

                g_x_0_y_0_yzzz_xyyz[k] = -g_x_0_y_0_zzz_xyyz[k] * ab_y + g_x_0_y_0_zzz_xyyyz[k];

                g_x_0_y_0_yzzz_xyzz[k] = -g_x_0_y_0_zzz_xyzz[k] * ab_y + g_x_0_y_0_zzz_xyyzz[k];

                g_x_0_y_0_yzzz_xzzz[k] = -g_x_0_y_0_zzz_xzzz[k] * ab_y + g_x_0_y_0_zzz_xyzzz[k];

                g_x_0_y_0_yzzz_yyyy[k] = -g_x_0_y_0_zzz_yyyy[k] * ab_y + g_x_0_y_0_zzz_yyyyy[k];

                g_x_0_y_0_yzzz_yyyz[k] = -g_x_0_y_0_zzz_yyyz[k] * ab_y + g_x_0_y_0_zzz_yyyyz[k];

                g_x_0_y_0_yzzz_yyzz[k] = -g_x_0_y_0_zzz_yyzz[k] * ab_y + g_x_0_y_0_zzz_yyyzz[k];

                g_x_0_y_0_yzzz_yzzz[k] = -g_x_0_y_0_zzz_yzzz[k] * ab_y + g_x_0_y_0_zzz_yyzzz[k];

                g_x_0_y_0_yzzz_zzzz[k] = -g_x_0_y_0_zzz_zzzz[k] * ab_y + g_x_0_y_0_zzz_yzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_zzzz_xxxx = cbuffer.data(gg_geom_1010_off + 435 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_xxxy = cbuffer.data(gg_geom_1010_off + 436 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_xxxz = cbuffer.data(gg_geom_1010_off + 437 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_xxyy = cbuffer.data(gg_geom_1010_off + 438 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_xxyz = cbuffer.data(gg_geom_1010_off + 439 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_xxzz = cbuffer.data(gg_geom_1010_off + 440 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_xyyy = cbuffer.data(gg_geom_1010_off + 441 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_xyyz = cbuffer.data(gg_geom_1010_off + 442 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_xyzz = cbuffer.data(gg_geom_1010_off + 443 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_xzzz = cbuffer.data(gg_geom_1010_off + 444 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_yyyy = cbuffer.data(gg_geom_1010_off + 445 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_yyyz = cbuffer.data(gg_geom_1010_off + 446 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_yyzz = cbuffer.data(gg_geom_1010_off + 447 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_yzzz = cbuffer.data(gg_geom_1010_off + 448 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_zzzz = cbuffer.data(gg_geom_1010_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_zzz_xxxx, g_x_0_y_0_zzz_xxxxz, g_x_0_y_0_zzz_xxxy, g_x_0_y_0_zzz_xxxyz, g_x_0_y_0_zzz_xxxz, g_x_0_y_0_zzz_xxxzz, g_x_0_y_0_zzz_xxyy, g_x_0_y_0_zzz_xxyyz, g_x_0_y_0_zzz_xxyz, g_x_0_y_0_zzz_xxyzz, g_x_0_y_0_zzz_xxzz, g_x_0_y_0_zzz_xxzzz, g_x_0_y_0_zzz_xyyy, g_x_0_y_0_zzz_xyyyz, g_x_0_y_0_zzz_xyyz, g_x_0_y_0_zzz_xyyzz, g_x_0_y_0_zzz_xyzz, g_x_0_y_0_zzz_xyzzz, g_x_0_y_0_zzz_xzzz, g_x_0_y_0_zzz_xzzzz, g_x_0_y_0_zzz_yyyy, g_x_0_y_0_zzz_yyyyz, g_x_0_y_0_zzz_yyyz, g_x_0_y_0_zzz_yyyzz, g_x_0_y_0_zzz_yyzz, g_x_0_y_0_zzz_yyzzz, g_x_0_y_0_zzz_yzzz, g_x_0_y_0_zzz_yzzzz, g_x_0_y_0_zzz_zzzz, g_x_0_y_0_zzz_zzzzz, g_x_0_y_0_zzzz_xxxx, g_x_0_y_0_zzzz_xxxy, g_x_0_y_0_zzzz_xxxz, g_x_0_y_0_zzzz_xxyy, g_x_0_y_0_zzzz_xxyz, g_x_0_y_0_zzzz_xxzz, g_x_0_y_0_zzzz_xyyy, g_x_0_y_0_zzzz_xyyz, g_x_0_y_0_zzzz_xyzz, g_x_0_y_0_zzzz_xzzz, g_x_0_y_0_zzzz_yyyy, g_x_0_y_0_zzzz_yyyz, g_x_0_y_0_zzzz_yyzz, g_x_0_y_0_zzzz_yzzz, g_x_0_y_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_zzzz_xxxx[k] = -g_x_0_y_0_zzz_xxxx[k] * ab_z + g_x_0_y_0_zzz_xxxxz[k];

                g_x_0_y_0_zzzz_xxxy[k] = -g_x_0_y_0_zzz_xxxy[k] * ab_z + g_x_0_y_0_zzz_xxxyz[k];

                g_x_0_y_0_zzzz_xxxz[k] = -g_x_0_y_0_zzz_xxxz[k] * ab_z + g_x_0_y_0_zzz_xxxzz[k];

                g_x_0_y_0_zzzz_xxyy[k] = -g_x_0_y_0_zzz_xxyy[k] * ab_z + g_x_0_y_0_zzz_xxyyz[k];

                g_x_0_y_0_zzzz_xxyz[k] = -g_x_0_y_0_zzz_xxyz[k] * ab_z + g_x_0_y_0_zzz_xxyzz[k];

                g_x_0_y_0_zzzz_xxzz[k] = -g_x_0_y_0_zzz_xxzz[k] * ab_z + g_x_0_y_0_zzz_xxzzz[k];

                g_x_0_y_0_zzzz_xyyy[k] = -g_x_0_y_0_zzz_xyyy[k] * ab_z + g_x_0_y_0_zzz_xyyyz[k];

                g_x_0_y_0_zzzz_xyyz[k] = -g_x_0_y_0_zzz_xyyz[k] * ab_z + g_x_0_y_0_zzz_xyyzz[k];

                g_x_0_y_0_zzzz_xyzz[k] = -g_x_0_y_0_zzz_xyzz[k] * ab_z + g_x_0_y_0_zzz_xyzzz[k];

                g_x_0_y_0_zzzz_xzzz[k] = -g_x_0_y_0_zzz_xzzz[k] * ab_z + g_x_0_y_0_zzz_xzzzz[k];

                g_x_0_y_0_zzzz_yyyy[k] = -g_x_0_y_0_zzz_yyyy[k] * ab_z + g_x_0_y_0_zzz_yyyyz[k];

                g_x_0_y_0_zzzz_yyyz[k] = -g_x_0_y_0_zzz_yyyz[k] * ab_z + g_x_0_y_0_zzz_yyyzz[k];

                g_x_0_y_0_zzzz_yyzz[k] = -g_x_0_y_0_zzz_yyzz[k] * ab_z + g_x_0_y_0_zzz_yyzzz[k];

                g_x_0_y_0_zzzz_yzzz[k] = -g_x_0_y_0_zzz_yzzz[k] * ab_z + g_x_0_y_0_zzz_yzzzz[k];

                g_x_0_y_0_zzzz_zzzz[k] = -g_x_0_y_0_zzz_zzzz[k] * ab_z + g_x_0_y_0_zzz_zzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxx_xxxx = cbuffer.data(gg_geom_1010_off + 450 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_xxxy = cbuffer.data(gg_geom_1010_off + 451 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_xxxz = cbuffer.data(gg_geom_1010_off + 452 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_xxyy = cbuffer.data(gg_geom_1010_off + 453 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_xxyz = cbuffer.data(gg_geom_1010_off + 454 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_xxzz = cbuffer.data(gg_geom_1010_off + 455 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_xyyy = cbuffer.data(gg_geom_1010_off + 456 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_xyyz = cbuffer.data(gg_geom_1010_off + 457 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_xyzz = cbuffer.data(gg_geom_1010_off + 458 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_xzzz = cbuffer.data(gg_geom_1010_off + 459 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_yyyy = cbuffer.data(gg_geom_1010_off + 460 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_yyyz = cbuffer.data(gg_geom_1010_off + 461 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_yyzz = cbuffer.data(gg_geom_1010_off + 462 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_yzzz = cbuffer.data(gg_geom_1010_off + 463 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_zzzz = cbuffer.data(gg_geom_1010_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_xxx_xxxx, g_0_0_z_0_xxx_xxxy, g_0_0_z_0_xxx_xxxz, g_0_0_z_0_xxx_xxyy, g_0_0_z_0_xxx_xxyz, g_0_0_z_0_xxx_xxzz, g_0_0_z_0_xxx_xyyy, g_0_0_z_0_xxx_xyyz, g_0_0_z_0_xxx_xyzz, g_0_0_z_0_xxx_xzzz, g_0_0_z_0_xxx_yyyy, g_0_0_z_0_xxx_yyyz, g_0_0_z_0_xxx_yyzz, g_0_0_z_0_xxx_yzzz, g_0_0_z_0_xxx_zzzz, g_x_0_z_0_xxx_xxxx, g_x_0_z_0_xxx_xxxxx, g_x_0_z_0_xxx_xxxxy, g_x_0_z_0_xxx_xxxxz, g_x_0_z_0_xxx_xxxy, g_x_0_z_0_xxx_xxxyy, g_x_0_z_0_xxx_xxxyz, g_x_0_z_0_xxx_xxxz, g_x_0_z_0_xxx_xxxzz, g_x_0_z_0_xxx_xxyy, g_x_0_z_0_xxx_xxyyy, g_x_0_z_0_xxx_xxyyz, g_x_0_z_0_xxx_xxyz, g_x_0_z_0_xxx_xxyzz, g_x_0_z_0_xxx_xxzz, g_x_0_z_0_xxx_xxzzz, g_x_0_z_0_xxx_xyyy, g_x_0_z_0_xxx_xyyyy, g_x_0_z_0_xxx_xyyyz, g_x_0_z_0_xxx_xyyz, g_x_0_z_0_xxx_xyyzz, g_x_0_z_0_xxx_xyzz, g_x_0_z_0_xxx_xyzzz, g_x_0_z_0_xxx_xzzz, g_x_0_z_0_xxx_xzzzz, g_x_0_z_0_xxx_yyyy, g_x_0_z_0_xxx_yyyz, g_x_0_z_0_xxx_yyzz, g_x_0_z_0_xxx_yzzz, g_x_0_z_0_xxx_zzzz, g_x_0_z_0_xxxx_xxxx, g_x_0_z_0_xxxx_xxxy, g_x_0_z_0_xxxx_xxxz, g_x_0_z_0_xxxx_xxyy, g_x_0_z_0_xxxx_xxyz, g_x_0_z_0_xxxx_xxzz, g_x_0_z_0_xxxx_xyyy, g_x_0_z_0_xxxx_xyyz, g_x_0_z_0_xxxx_xyzz, g_x_0_z_0_xxxx_xzzz, g_x_0_z_0_xxxx_yyyy, g_x_0_z_0_xxxx_yyyz, g_x_0_z_0_xxxx_yyzz, g_x_0_z_0_xxxx_yzzz, g_x_0_z_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxx_xxxx[k] = -g_0_0_z_0_xxx_xxxx[k] - g_x_0_z_0_xxx_xxxx[k] * ab_x + g_x_0_z_0_xxx_xxxxx[k];

                g_x_0_z_0_xxxx_xxxy[k] = -g_0_0_z_0_xxx_xxxy[k] - g_x_0_z_0_xxx_xxxy[k] * ab_x + g_x_0_z_0_xxx_xxxxy[k];

                g_x_0_z_0_xxxx_xxxz[k] = -g_0_0_z_0_xxx_xxxz[k] - g_x_0_z_0_xxx_xxxz[k] * ab_x + g_x_0_z_0_xxx_xxxxz[k];

                g_x_0_z_0_xxxx_xxyy[k] = -g_0_0_z_0_xxx_xxyy[k] - g_x_0_z_0_xxx_xxyy[k] * ab_x + g_x_0_z_0_xxx_xxxyy[k];

                g_x_0_z_0_xxxx_xxyz[k] = -g_0_0_z_0_xxx_xxyz[k] - g_x_0_z_0_xxx_xxyz[k] * ab_x + g_x_0_z_0_xxx_xxxyz[k];

                g_x_0_z_0_xxxx_xxzz[k] = -g_0_0_z_0_xxx_xxzz[k] - g_x_0_z_0_xxx_xxzz[k] * ab_x + g_x_0_z_0_xxx_xxxzz[k];

                g_x_0_z_0_xxxx_xyyy[k] = -g_0_0_z_0_xxx_xyyy[k] - g_x_0_z_0_xxx_xyyy[k] * ab_x + g_x_0_z_0_xxx_xxyyy[k];

                g_x_0_z_0_xxxx_xyyz[k] = -g_0_0_z_0_xxx_xyyz[k] - g_x_0_z_0_xxx_xyyz[k] * ab_x + g_x_0_z_0_xxx_xxyyz[k];

                g_x_0_z_0_xxxx_xyzz[k] = -g_0_0_z_0_xxx_xyzz[k] - g_x_0_z_0_xxx_xyzz[k] * ab_x + g_x_0_z_0_xxx_xxyzz[k];

                g_x_0_z_0_xxxx_xzzz[k] = -g_0_0_z_0_xxx_xzzz[k] - g_x_0_z_0_xxx_xzzz[k] * ab_x + g_x_0_z_0_xxx_xxzzz[k];

                g_x_0_z_0_xxxx_yyyy[k] = -g_0_0_z_0_xxx_yyyy[k] - g_x_0_z_0_xxx_yyyy[k] * ab_x + g_x_0_z_0_xxx_xyyyy[k];

                g_x_0_z_0_xxxx_yyyz[k] = -g_0_0_z_0_xxx_yyyz[k] - g_x_0_z_0_xxx_yyyz[k] * ab_x + g_x_0_z_0_xxx_xyyyz[k];

                g_x_0_z_0_xxxx_yyzz[k] = -g_0_0_z_0_xxx_yyzz[k] - g_x_0_z_0_xxx_yyzz[k] * ab_x + g_x_0_z_0_xxx_xyyzz[k];

                g_x_0_z_0_xxxx_yzzz[k] = -g_0_0_z_0_xxx_yzzz[k] - g_x_0_z_0_xxx_yzzz[k] * ab_x + g_x_0_z_0_xxx_xyzzz[k];

                g_x_0_z_0_xxxx_zzzz[k] = -g_0_0_z_0_xxx_zzzz[k] - g_x_0_z_0_xxx_zzzz[k] * ab_x + g_x_0_z_0_xxx_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxy_xxxx = cbuffer.data(gg_geom_1010_off + 465 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_xxxy = cbuffer.data(gg_geom_1010_off + 466 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_xxxz = cbuffer.data(gg_geom_1010_off + 467 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_xxyy = cbuffer.data(gg_geom_1010_off + 468 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_xxyz = cbuffer.data(gg_geom_1010_off + 469 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_xxzz = cbuffer.data(gg_geom_1010_off + 470 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_xyyy = cbuffer.data(gg_geom_1010_off + 471 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_xyyz = cbuffer.data(gg_geom_1010_off + 472 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_xyzz = cbuffer.data(gg_geom_1010_off + 473 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_xzzz = cbuffer.data(gg_geom_1010_off + 474 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_yyyy = cbuffer.data(gg_geom_1010_off + 475 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_yyyz = cbuffer.data(gg_geom_1010_off + 476 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_yyzz = cbuffer.data(gg_geom_1010_off + 477 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_yzzz = cbuffer.data(gg_geom_1010_off + 478 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_zzzz = cbuffer.data(gg_geom_1010_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxx_xxxx, g_x_0_z_0_xxx_xxxxy, g_x_0_z_0_xxx_xxxy, g_x_0_z_0_xxx_xxxyy, g_x_0_z_0_xxx_xxxyz, g_x_0_z_0_xxx_xxxz, g_x_0_z_0_xxx_xxyy, g_x_0_z_0_xxx_xxyyy, g_x_0_z_0_xxx_xxyyz, g_x_0_z_0_xxx_xxyz, g_x_0_z_0_xxx_xxyzz, g_x_0_z_0_xxx_xxzz, g_x_0_z_0_xxx_xyyy, g_x_0_z_0_xxx_xyyyy, g_x_0_z_0_xxx_xyyyz, g_x_0_z_0_xxx_xyyz, g_x_0_z_0_xxx_xyyzz, g_x_0_z_0_xxx_xyzz, g_x_0_z_0_xxx_xyzzz, g_x_0_z_0_xxx_xzzz, g_x_0_z_0_xxx_yyyy, g_x_0_z_0_xxx_yyyyy, g_x_0_z_0_xxx_yyyyz, g_x_0_z_0_xxx_yyyz, g_x_0_z_0_xxx_yyyzz, g_x_0_z_0_xxx_yyzz, g_x_0_z_0_xxx_yyzzz, g_x_0_z_0_xxx_yzzz, g_x_0_z_0_xxx_yzzzz, g_x_0_z_0_xxx_zzzz, g_x_0_z_0_xxxy_xxxx, g_x_0_z_0_xxxy_xxxy, g_x_0_z_0_xxxy_xxxz, g_x_0_z_0_xxxy_xxyy, g_x_0_z_0_xxxy_xxyz, g_x_0_z_0_xxxy_xxzz, g_x_0_z_0_xxxy_xyyy, g_x_0_z_0_xxxy_xyyz, g_x_0_z_0_xxxy_xyzz, g_x_0_z_0_xxxy_xzzz, g_x_0_z_0_xxxy_yyyy, g_x_0_z_0_xxxy_yyyz, g_x_0_z_0_xxxy_yyzz, g_x_0_z_0_xxxy_yzzz, g_x_0_z_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxy_xxxx[k] = -g_x_0_z_0_xxx_xxxx[k] * ab_y + g_x_0_z_0_xxx_xxxxy[k];

                g_x_0_z_0_xxxy_xxxy[k] = -g_x_0_z_0_xxx_xxxy[k] * ab_y + g_x_0_z_0_xxx_xxxyy[k];

                g_x_0_z_0_xxxy_xxxz[k] = -g_x_0_z_0_xxx_xxxz[k] * ab_y + g_x_0_z_0_xxx_xxxyz[k];

                g_x_0_z_0_xxxy_xxyy[k] = -g_x_0_z_0_xxx_xxyy[k] * ab_y + g_x_0_z_0_xxx_xxyyy[k];

                g_x_0_z_0_xxxy_xxyz[k] = -g_x_0_z_0_xxx_xxyz[k] * ab_y + g_x_0_z_0_xxx_xxyyz[k];

                g_x_0_z_0_xxxy_xxzz[k] = -g_x_0_z_0_xxx_xxzz[k] * ab_y + g_x_0_z_0_xxx_xxyzz[k];

                g_x_0_z_0_xxxy_xyyy[k] = -g_x_0_z_0_xxx_xyyy[k] * ab_y + g_x_0_z_0_xxx_xyyyy[k];

                g_x_0_z_0_xxxy_xyyz[k] = -g_x_0_z_0_xxx_xyyz[k] * ab_y + g_x_0_z_0_xxx_xyyyz[k];

                g_x_0_z_0_xxxy_xyzz[k] = -g_x_0_z_0_xxx_xyzz[k] * ab_y + g_x_0_z_0_xxx_xyyzz[k];

                g_x_0_z_0_xxxy_xzzz[k] = -g_x_0_z_0_xxx_xzzz[k] * ab_y + g_x_0_z_0_xxx_xyzzz[k];

                g_x_0_z_0_xxxy_yyyy[k] = -g_x_0_z_0_xxx_yyyy[k] * ab_y + g_x_0_z_0_xxx_yyyyy[k];

                g_x_0_z_0_xxxy_yyyz[k] = -g_x_0_z_0_xxx_yyyz[k] * ab_y + g_x_0_z_0_xxx_yyyyz[k];

                g_x_0_z_0_xxxy_yyzz[k] = -g_x_0_z_0_xxx_yyzz[k] * ab_y + g_x_0_z_0_xxx_yyyzz[k];

                g_x_0_z_0_xxxy_yzzz[k] = -g_x_0_z_0_xxx_yzzz[k] * ab_y + g_x_0_z_0_xxx_yyzzz[k];

                g_x_0_z_0_xxxy_zzzz[k] = -g_x_0_z_0_xxx_zzzz[k] * ab_y + g_x_0_z_0_xxx_yzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxz_xxxx = cbuffer.data(gg_geom_1010_off + 480 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_xxxy = cbuffer.data(gg_geom_1010_off + 481 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_xxxz = cbuffer.data(gg_geom_1010_off + 482 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_xxyy = cbuffer.data(gg_geom_1010_off + 483 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_xxyz = cbuffer.data(gg_geom_1010_off + 484 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_xxzz = cbuffer.data(gg_geom_1010_off + 485 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_xyyy = cbuffer.data(gg_geom_1010_off + 486 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_xyyz = cbuffer.data(gg_geom_1010_off + 487 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_xyzz = cbuffer.data(gg_geom_1010_off + 488 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_xzzz = cbuffer.data(gg_geom_1010_off + 489 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_yyyy = cbuffer.data(gg_geom_1010_off + 490 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_yyyz = cbuffer.data(gg_geom_1010_off + 491 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_yyzz = cbuffer.data(gg_geom_1010_off + 492 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_yzzz = cbuffer.data(gg_geom_1010_off + 493 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_zzzz = cbuffer.data(gg_geom_1010_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxx_xxxx, g_x_0_z_0_xxx_xxxxz, g_x_0_z_0_xxx_xxxy, g_x_0_z_0_xxx_xxxyz, g_x_0_z_0_xxx_xxxz, g_x_0_z_0_xxx_xxxzz, g_x_0_z_0_xxx_xxyy, g_x_0_z_0_xxx_xxyyz, g_x_0_z_0_xxx_xxyz, g_x_0_z_0_xxx_xxyzz, g_x_0_z_0_xxx_xxzz, g_x_0_z_0_xxx_xxzzz, g_x_0_z_0_xxx_xyyy, g_x_0_z_0_xxx_xyyyz, g_x_0_z_0_xxx_xyyz, g_x_0_z_0_xxx_xyyzz, g_x_0_z_0_xxx_xyzz, g_x_0_z_0_xxx_xyzzz, g_x_0_z_0_xxx_xzzz, g_x_0_z_0_xxx_xzzzz, g_x_0_z_0_xxx_yyyy, g_x_0_z_0_xxx_yyyyz, g_x_0_z_0_xxx_yyyz, g_x_0_z_0_xxx_yyyzz, g_x_0_z_0_xxx_yyzz, g_x_0_z_0_xxx_yyzzz, g_x_0_z_0_xxx_yzzz, g_x_0_z_0_xxx_yzzzz, g_x_0_z_0_xxx_zzzz, g_x_0_z_0_xxx_zzzzz, g_x_0_z_0_xxxz_xxxx, g_x_0_z_0_xxxz_xxxy, g_x_0_z_0_xxxz_xxxz, g_x_0_z_0_xxxz_xxyy, g_x_0_z_0_xxxz_xxyz, g_x_0_z_0_xxxz_xxzz, g_x_0_z_0_xxxz_xyyy, g_x_0_z_0_xxxz_xyyz, g_x_0_z_0_xxxz_xyzz, g_x_0_z_0_xxxz_xzzz, g_x_0_z_0_xxxz_yyyy, g_x_0_z_0_xxxz_yyyz, g_x_0_z_0_xxxz_yyzz, g_x_0_z_0_xxxz_yzzz, g_x_0_z_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxz_xxxx[k] = -g_x_0_z_0_xxx_xxxx[k] * ab_z + g_x_0_z_0_xxx_xxxxz[k];

                g_x_0_z_0_xxxz_xxxy[k] = -g_x_0_z_0_xxx_xxxy[k] * ab_z + g_x_0_z_0_xxx_xxxyz[k];

                g_x_0_z_0_xxxz_xxxz[k] = -g_x_0_z_0_xxx_xxxz[k] * ab_z + g_x_0_z_0_xxx_xxxzz[k];

                g_x_0_z_0_xxxz_xxyy[k] = -g_x_0_z_0_xxx_xxyy[k] * ab_z + g_x_0_z_0_xxx_xxyyz[k];

                g_x_0_z_0_xxxz_xxyz[k] = -g_x_0_z_0_xxx_xxyz[k] * ab_z + g_x_0_z_0_xxx_xxyzz[k];

                g_x_0_z_0_xxxz_xxzz[k] = -g_x_0_z_0_xxx_xxzz[k] * ab_z + g_x_0_z_0_xxx_xxzzz[k];

                g_x_0_z_0_xxxz_xyyy[k] = -g_x_0_z_0_xxx_xyyy[k] * ab_z + g_x_0_z_0_xxx_xyyyz[k];

                g_x_0_z_0_xxxz_xyyz[k] = -g_x_0_z_0_xxx_xyyz[k] * ab_z + g_x_0_z_0_xxx_xyyzz[k];

                g_x_0_z_0_xxxz_xyzz[k] = -g_x_0_z_0_xxx_xyzz[k] * ab_z + g_x_0_z_0_xxx_xyzzz[k];

                g_x_0_z_0_xxxz_xzzz[k] = -g_x_0_z_0_xxx_xzzz[k] * ab_z + g_x_0_z_0_xxx_xzzzz[k];

                g_x_0_z_0_xxxz_yyyy[k] = -g_x_0_z_0_xxx_yyyy[k] * ab_z + g_x_0_z_0_xxx_yyyyz[k];

                g_x_0_z_0_xxxz_yyyz[k] = -g_x_0_z_0_xxx_yyyz[k] * ab_z + g_x_0_z_0_xxx_yyyzz[k];

                g_x_0_z_0_xxxz_yyzz[k] = -g_x_0_z_0_xxx_yyzz[k] * ab_z + g_x_0_z_0_xxx_yyzzz[k];

                g_x_0_z_0_xxxz_yzzz[k] = -g_x_0_z_0_xxx_yzzz[k] * ab_z + g_x_0_z_0_xxx_yzzzz[k];

                g_x_0_z_0_xxxz_zzzz[k] = -g_x_0_z_0_xxx_zzzz[k] * ab_z + g_x_0_z_0_xxx_zzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyy_xxxx = cbuffer.data(gg_geom_1010_off + 495 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_xxxy = cbuffer.data(gg_geom_1010_off + 496 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_xxxz = cbuffer.data(gg_geom_1010_off + 497 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_xxyy = cbuffer.data(gg_geom_1010_off + 498 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_xxyz = cbuffer.data(gg_geom_1010_off + 499 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_xxzz = cbuffer.data(gg_geom_1010_off + 500 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_xyyy = cbuffer.data(gg_geom_1010_off + 501 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_xyyz = cbuffer.data(gg_geom_1010_off + 502 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_xyzz = cbuffer.data(gg_geom_1010_off + 503 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_xzzz = cbuffer.data(gg_geom_1010_off + 504 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_yyyy = cbuffer.data(gg_geom_1010_off + 505 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_yyyz = cbuffer.data(gg_geom_1010_off + 506 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_yyzz = cbuffer.data(gg_geom_1010_off + 507 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_yzzz = cbuffer.data(gg_geom_1010_off + 508 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_zzzz = cbuffer.data(gg_geom_1010_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxy_xxxx, g_x_0_z_0_xxy_xxxxy, g_x_0_z_0_xxy_xxxy, g_x_0_z_0_xxy_xxxyy, g_x_0_z_0_xxy_xxxyz, g_x_0_z_0_xxy_xxxz, g_x_0_z_0_xxy_xxyy, g_x_0_z_0_xxy_xxyyy, g_x_0_z_0_xxy_xxyyz, g_x_0_z_0_xxy_xxyz, g_x_0_z_0_xxy_xxyzz, g_x_0_z_0_xxy_xxzz, g_x_0_z_0_xxy_xyyy, g_x_0_z_0_xxy_xyyyy, g_x_0_z_0_xxy_xyyyz, g_x_0_z_0_xxy_xyyz, g_x_0_z_0_xxy_xyyzz, g_x_0_z_0_xxy_xyzz, g_x_0_z_0_xxy_xyzzz, g_x_0_z_0_xxy_xzzz, g_x_0_z_0_xxy_yyyy, g_x_0_z_0_xxy_yyyyy, g_x_0_z_0_xxy_yyyyz, g_x_0_z_0_xxy_yyyz, g_x_0_z_0_xxy_yyyzz, g_x_0_z_0_xxy_yyzz, g_x_0_z_0_xxy_yyzzz, g_x_0_z_0_xxy_yzzz, g_x_0_z_0_xxy_yzzzz, g_x_0_z_0_xxy_zzzz, g_x_0_z_0_xxyy_xxxx, g_x_0_z_0_xxyy_xxxy, g_x_0_z_0_xxyy_xxxz, g_x_0_z_0_xxyy_xxyy, g_x_0_z_0_xxyy_xxyz, g_x_0_z_0_xxyy_xxzz, g_x_0_z_0_xxyy_xyyy, g_x_0_z_0_xxyy_xyyz, g_x_0_z_0_xxyy_xyzz, g_x_0_z_0_xxyy_xzzz, g_x_0_z_0_xxyy_yyyy, g_x_0_z_0_xxyy_yyyz, g_x_0_z_0_xxyy_yyzz, g_x_0_z_0_xxyy_yzzz, g_x_0_z_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyy_xxxx[k] = -g_x_0_z_0_xxy_xxxx[k] * ab_y + g_x_0_z_0_xxy_xxxxy[k];

                g_x_0_z_0_xxyy_xxxy[k] = -g_x_0_z_0_xxy_xxxy[k] * ab_y + g_x_0_z_0_xxy_xxxyy[k];

                g_x_0_z_0_xxyy_xxxz[k] = -g_x_0_z_0_xxy_xxxz[k] * ab_y + g_x_0_z_0_xxy_xxxyz[k];

                g_x_0_z_0_xxyy_xxyy[k] = -g_x_0_z_0_xxy_xxyy[k] * ab_y + g_x_0_z_0_xxy_xxyyy[k];

                g_x_0_z_0_xxyy_xxyz[k] = -g_x_0_z_0_xxy_xxyz[k] * ab_y + g_x_0_z_0_xxy_xxyyz[k];

                g_x_0_z_0_xxyy_xxzz[k] = -g_x_0_z_0_xxy_xxzz[k] * ab_y + g_x_0_z_0_xxy_xxyzz[k];

                g_x_0_z_0_xxyy_xyyy[k] = -g_x_0_z_0_xxy_xyyy[k] * ab_y + g_x_0_z_0_xxy_xyyyy[k];

                g_x_0_z_0_xxyy_xyyz[k] = -g_x_0_z_0_xxy_xyyz[k] * ab_y + g_x_0_z_0_xxy_xyyyz[k];

                g_x_0_z_0_xxyy_xyzz[k] = -g_x_0_z_0_xxy_xyzz[k] * ab_y + g_x_0_z_0_xxy_xyyzz[k];

                g_x_0_z_0_xxyy_xzzz[k] = -g_x_0_z_0_xxy_xzzz[k] * ab_y + g_x_0_z_0_xxy_xyzzz[k];

                g_x_0_z_0_xxyy_yyyy[k] = -g_x_0_z_0_xxy_yyyy[k] * ab_y + g_x_0_z_0_xxy_yyyyy[k];

                g_x_0_z_0_xxyy_yyyz[k] = -g_x_0_z_0_xxy_yyyz[k] * ab_y + g_x_0_z_0_xxy_yyyyz[k];

                g_x_0_z_0_xxyy_yyzz[k] = -g_x_0_z_0_xxy_yyzz[k] * ab_y + g_x_0_z_0_xxy_yyyzz[k];

                g_x_0_z_0_xxyy_yzzz[k] = -g_x_0_z_0_xxy_yzzz[k] * ab_y + g_x_0_z_0_xxy_yyzzz[k];

                g_x_0_z_0_xxyy_zzzz[k] = -g_x_0_z_0_xxy_zzzz[k] * ab_y + g_x_0_z_0_xxy_yzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyz_xxxx = cbuffer.data(gg_geom_1010_off + 510 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_xxxy = cbuffer.data(gg_geom_1010_off + 511 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_xxxz = cbuffer.data(gg_geom_1010_off + 512 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_xxyy = cbuffer.data(gg_geom_1010_off + 513 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_xxyz = cbuffer.data(gg_geom_1010_off + 514 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_xxzz = cbuffer.data(gg_geom_1010_off + 515 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_xyyy = cbuffer.data(gg_geom_1010_off + 516 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_xyyz = cbuffer.data(gg_geom_1010_off + 517 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_xyzz = cbuffer.data(gg_geom_1010_off + 518 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_xzzz = cbuffer.data(gg_geom_1010_off + 519 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_yyyy = cbuffer.data(gg_geom_1010_off + 520 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_yyyz = cbuffer.data(gg_geom_1010_off + 521 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_yyzz = cbuffer.data(gg_geom_1010_off + 522 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_yzzz = cbuffer.data(gg_geom_1010_off + 523 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_zzzz = cbuffer.data(gg_geom_1010_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyz_xxxx, g_x_0_z_0_xxyz_xxxy, g_x_0_z_0_xxyz_xxxz, g_x_0_z_0_xxyz_xxyy, g_x_0_z_0_xxyz_xxyz, g_x_0_z_0_xxyz_xxzz, g_x_0_z_0_xxyz_xyyy, g_x_0_z_0_xxyz_xyyz, g_x_0_z_0_xxyz_xyzz, g_x_0_z_0_xxyz_xzzz, g_x_0_z_0_xxyz_yyyy, g_x_0_z_0_xxyz_yyyz, g_x_0_z_0_xxyz_yyzz, g_x_0_z_0_xxyz_yzzz, g_x_0_z_0_xxyz_zzzz, g_x_0_z_0_xxz_xxxx, g_x_0_z_0_xxz_xxxxy, g_x_0_z_0_xxz_xxxy, g_x_0_z_0_xxz_xxxyy, g_x_0_z_0_xxz_xxxyz, g_x_0_z_0_xxz_xxxz, g_x_0_z_0_xxz_xxyy, g_x_0_z_0_xxz_xxyyy, g_x_0_z_0_xxz_xxyyz, g_x_0_z_0_xxz_xxyz, g_x_0_z_0_xxz_xxyzz, g_x_0_z_0_xxz_xxzz, g_x_0_z_0_xxz_xyyy, g_x_0_z_0_xxz_xyyyy, g_x_0_z_0_xxz_xyyyz, g_x_0_z_0_xxz_xyyz, g_x_0_z_0_xxz_xyyzz, g_x_0_z_0_xxz_xyzz, g_x_0_z_0_xxz_xyzzz, g_x_0_z_0_xxz_xzzz, g_x_0_z_0_xxz_yyyy, g_x_0_z_0_xxz_yyyyy, g_x_0_z_0_xxz_yyyyz, g_x_0_z_0_xxz_yyyz, g_x_0_z_0_xxz_yyyzz, g_x_0_z_0_xxz_yyzz, g_x_0_z_0_xxz_yyzzz, g_x_0_z_0_xxz_yzzz, g_x_0_z_0_xxz_yzzzz, g_x_0_z_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyz_xxxx[k] = -g_x_0_z_0_xxz_xxxx[k] * ab_y + g_x_0_z_0_xxz_xxxxy[k];

                g_x_0_z_0_xxyz_xxxy[k] = -g_x_0_z_0_xxz_xxxy[k] * ab_y + g_x_0_z_0_xxz_xxxyy[k];

                g_x_0_z_0_xxyz_xxxz[k] = -g_x_0_z_0_xxz_xxxz[k] * ab_y + g_x_0_z_0_xxz_xxxyz[k];

                g_x_0_z_0_xxyz_xxyy[k] = -g_x_0_z_0_xxz_xxyy[k] * ab_y + g_x_0_z_0_xxz_xxyyy[k];

                g_x_0_z_0_xxyz_xxyz[k] = -g_x_0_z_0_xxz_xxyz[k] * ab_y + g_x_0_z_0_xxz_xxyyz[k];

                g_x_0_z_0_xxyz_xxzz[k] = -g_x_0_z_0_xxz_xxzz[k] * ab_y + g_x_0_z_0_xxz_xxyzz[k];

                g_x_0_z_0_xxyz_xyyy[k] = -g_x_0_z_0_xxz_xyyy[k] * ab_y + g_x_0_z_0_xxz_xyyyy[k];

                g_x_0_z_0_xxyz_xyyz[k] = -g_x_0_z_0_xxz_xyyz[k] * ab_y + g_x_0_z_0_xxz_xyyyz[k];

                g_x_0_z_0_xxyz_xyzz[k] = -g_x_0_z_0_xxz_xyzz[k] * ab_y + g_x_0_z_0_xxz_xyyzz[k];

                g_x_0_z_0_xxyz_xzzz[k] = -g_x_0_z_0_xxz_xzzz[k] * ab_y + g_x_0_z_0_xxz_xyzzz[k];

                g_x_0_z_0_xxyz_yyyy[k] = -g_x_0_z_0_xxz_yyyy[k] * ab_y + g_x_0_z_0_xxz_yyyyy[k];

                g_x_0_z_0_xxyz_yyyz[k] = -g_x_0_z_0_xxz_yyyz[k] * ab_y + g_x_0_z_0_xxz_yyyyz[k];

                g_x_0_z_0_xxyz_yyzz[k] = -g_x_0_z_0_xxz_yyzz[k] * ab_y + g_x_0_z_0_xxz_yyyzz[k];

                g_x_0_z_0_xxyz_yzzz[k] = -g_x_0_z_0_xxz_yzzz[k] * ab_y + g_x_0_z_0_xxz_yyzzz[k];

                g_x_0_z_0_xxyz_zzzz[k] = -g_x_0_z_0_xxz_zzzz[k] * ab_y + g_x_0_z_0_xxz_yzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxzz_xxxx = cbuffer.data(gg_geom_1010_off + 525 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_xxxy = cbuffer.data(gg_geom_1010_off + 526 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_xxxz = cbuffer.data(gg_geom_1010_off + 527 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_xxyy = cbuffer.data(gg_geom_1010_off + 528 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_xxyz = cbuffer.data(gg_geom_1010_off + 529 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_xxzz = cbuffer.data(gg_geom_1010_off + 530 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_xyyy = cbuffer.data(gg_geom_1010_off + 531 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_xyyz = cbuffer.data(gg_geom_1010_off + 532 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_xyzz = cbuffer.data(gg_geom_1010_off + 533 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_xzzz = cbuffer.data(gg_geom_1010_off + 534 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_yyyy = cbuffer.data(gg_geom_1010_off + 535 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_yyyz = cbuffer.data(gg_geom_1010_off + 536 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_yyzz = cbuffer.data(gg_geom_1010_off + 537 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_yzzz = cbuffer.data(gg_geom_1010_off + 538 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_zzzz = cbuffer.data(gg_geom_1010_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxz_xxxx, g_x_0_z_0_xxz_xxxxz, g_x_0_z_0_xxz_xxxy, g_x_0_z_0_xxz_xxxyz, g_x_0_z_0_xxz_xxxz, g_x_0_z_0_xxz_xxxzz, g_x_0_z_0_xxz_xxyy, g_x_0_z_0_xxz_xxyyz, g_x_0_z_0_xxz_xxyz, g_x_0_z_0_xxz_xxyzz, g_x_0_z_0_xxz_xxzz, g_x_0_z_0_xxz_xxzzz, g_x_0_z_0_xxz_xyyy, g_x_0_z_0_xxz_xyyyz, g_x_0_z_0_xxz_xyyz, g_x_0_z_0_xxz_xyyzz, g_x_0_z_0_xxz_xyzz, g_x_0_z_0_xxz_xyzzz, g_x_0_z_0_xxz_xzzz, g_x_0_z_0_xxz_xzzzz, g_x_0_z_0_xxz_yyyy, g_x_0_z_0_xxz_yyyyz, g_x_0_z_0_xxz_yyyz, g_x_0_z_0_xxz_yyyzz, g_x_0_z_0_xxz_yyzz, g_x_0_z_0_xxz_yyzzz, g_x_0_z_0_xxz_yzzz, g_x_0_z_0_xxz_yzzzz, g_x_0_z_0_xxz_zzzz, g_x_0_z_0_xxz_zzzzz, g_x_0_z_0_xxzz_xxxx, g_x_0_z_0_xxzz_xxxy, g_x_0_z_0_xxzz_xxxz, g_x_0_z_0_xxzz_xxyy, g_x_0_z_0_xxzz_xxyz, g_x_0_z_0_xxzz_xxzz, g_x_0_z_0_xxzz_xyyy, g_x_0_z_0_xxzz_xyyz, g_x_0_z_0_xxzz_xyzz, g_x_0_z_0_xxzz_xzzz, g_x_0_z_0_xxzz_yyyy, g_x_0_z_0_xxzz_yyyz, g_x_0_z_0_xxzz_yyzz, g_x_0_z_0_xxzz_yzzz, g_x_0_z_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxzz_xxxx[k] = -g_x_0_z_0_xxz_xxxx[k] * ab_z + g_x_0_z_0_xxz_xxxxz[k];

                g_x_0_z_0_xxzz_xxxy[k] = -g_x_0_z_0_xxz_xxxy[k] * ab_z + g_x_0_z_0_xxz_xxxyz[k];

                g_x_0_z_0_xxzz_xxxz[k] = -g_x_0_z_0_xxz_xxxz[k] * ab_z + g_x_0_z_0_xxz_xxxzz[k];

                g_x_0_z_0_xxzz_xxyy[k] = -g_x_0_z_0_xxz_xxyy[k] * ab_z + g_x_0_z_0_xxz_xxyyz[k];

                g_x_0_z_0_xxzz_xxyz[k] = -g_x_0_z_0_xxz_xxyz[k] * ab_z + g_x_0_z_0_xxz_xxyzz[k];

                g_x_0_z_0_xxzz_xxzz[k] = -g_x_0_z_0_xxz_xxzz[k] * ab_z + g_x_0_z_0_xxz_xxzzz[k];

                g_x_0_z_0_xxzz_xyyy[k] = -g_x_0_z_0_xxz_xyyy[k] * ab_z + g_x_0_z_0_xxz_xyyyz[k];

                g_x_0_z_0_xxzz_xyyz[k] = -g_x_0_z_0_xxz_xyyz[k] * ab_z + g_x_0_z_0_xxz_xyyzz[k];

                g_x_0_z_0_xxzz_xyzz[k] = -g_x_0_z_0_xxz_xyzz[k] * ab_z + g_x_0_z_0_xxz_xyzzz[k];

                g_x_0_z_0_xxzz_xzzz[k] = -g_x_0_z_0_xxz_xzzz[k] * ab_z + g_x_0_z_0_xxz_xzzzz[k];

                g_x_0_z_0_xxzz_yyyy[k] = -g_x_0_z_0_xxz_yyyy[k] * ab_z + g_x_0_z_0_xxz_yyyyz[k];

                g_x_0_z_0_xxzz_yyyz[k] = -g_x_0_z_0_xxz_yyyz[k] * ab_z + g_x_0_z_0_xxz_yyyzz[k];

                g_x_0_z_0_xxzz_yyzz[k] = -g_x_0_z_0_xxz_yyzz[k] * ab_z + g_x_0_z_0_xxz_yyzzz[k];

                g_x_0_z_0_xxzz_yzzz[k] = -g_x_0_z_0_xxz_yzzz[k] * ab_z + g_x_0_z_0_xxz_yzzzz[k];

                g_x_0_z_0_xxzz_zzzz[k] = -g_x_0_z_0_xxz_zzzz[k] * ab_z + g_x_0_z_0_xxz_zzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyy_xxxx = cbuffer.data(gg_geom_1010_off + 540 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_xxxy = cbuffer.data(gg_geom_1010_off + 541 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_xxxz = cbuffer.data(gg_geom_1010_off + 542 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_xxyy = cbuffer.data(gg_geom_1010_off + 543 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_xxyz = cbuffer.data(gg_geom_1010_off + 544 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_xxzz = cbuffer.data(gg_geom_1010_off + 545 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_xyyy = cbuffer.data(gg_geom_1010_off + 546 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_xyyz = cbuffer.data(gg_geom_1010_off + 547 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_xyzz = cbuffer.data(gg_geom_1010_off + 548 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_xzzz = cbuffer.data(gg_geom_1010_off + 549 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_yyyy = cbuffer.data(gg_geom_1010_off + 550 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_yyyz = cbuffer.data(gg_geom_1010_off + 551 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_yyzz = cbuffer.data(gg_geom_1010_off + 552 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_yzzz = cbuffer.data(gg_geom_1010_off + 553 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_zzzz = cbuffer.data(gg_geom_1010_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyy_xxxx, g_x_0_z_0_xyy_xxxxy, g_x_0_z_0_xyy_xxxy, g_x_0_z_0_xyy_xxxyy, g_x_0_z_0_xyy_xxxyz, g_x_0_z_0_xyy_xxxz, g_x_0_z_0_xyy_xxyy, g_x_0_z_0_xyy_xxyyy, g_x_0_z_0_xyy_xxyyz, g_x_0_z_0_xyy_xxyz, g_x_0_z_0_xyy_xxyzz, g_x_0_z_0_xyy_xxzz, g_x_0_z_0_xyy_xyyy, g_x_0_z_0_xyy_xyyyy, g_x_0_z_0_xyy_xyyyz, g_x_0_z_0_xyy_xyyz, g_x_0_z_0_xyy_xyyzz, g_x_0_z_0_xyy_xyzz, g_x_0_z_0_xyy_xyzzz, g_x_0_z_0_xyy_xzzz, g_x_0_z_0_xyy_yyyy, g_x_0_z_0_xyy_yyyyy, g_x_0_z_0_xyy_yyyyz, g_x_0_z_0_xyy_yyyz, g_x_0_z_0_xyy_yyyzz, g_x_0_z_0_xyy_yyzz, g_x_0_z_0_xyy_yyzzz, g_x_0_z_0_xyy_yzzz, g_x_0_z_0_xyy_yzzzz, g_x_0_z_0_xyy_zzzz, g_x_0_z_0_xyyy_xxxx, g_x_0_z_0_xyyy_xxxy, g_x_0_z_0_xyyy_xxxz, g_x_0_z_0_xyyy_xxyy, g_x_0_z_0_xyyy_xxyz, g_x_0_z_0_xyyy_xxzz, g_x_0_z_0_xyyy_xyyy, g_x_0_z_0_xyyy_xyyz, g_x_0_z_0_xyyy_xyzz, g_x_0_z_0_xyyy_xzzz, g_x_0_z_0_xyyy_yyyy, g_x_0_z_0_xyyy_yyyz, g_x_0_z_0_xyyy_yyzz, g_x_0_z_0_xyyy_yzzz, g_x_0_z_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyy_xxxx[k] = -g_x_0_z_0_xyy_xxxx[k] * ab_y + g_x_0_z_0_xyy_xxxxy[k];

                g_x_0_z_0_xyyy_xxxy[k] = -g_x_0_z_0_xyy_xxxy[k] * ab_y + g_x_0_z_0_xyy_xxxyy[k];

                g_x_0_z_0_xyyy_xxxz[k] = -g_x_0_z_0_xyy_xxxz[k] * ab_y + g_x_0_z_0_xyy_xxxyz[k];

                g_x_0_z_0_xyyy_xxyy[k] = -g_x_0_z_0_xyy_xxyy[k] * ab_y + g_x_0_z_0_xyy_xxyyy[k];

                g_x_0_z_0_xyyy_xxyz[k] = -g_x_0_z_0_xyy_xxyz[k] * ab_y + g_x_0_z_0_xyy_xxyyz[k];

                g_x_0_z_0_xyyy_xxzz[k] = -g_x_0_z_0_xyy_xxzz[k] * ab_y + g_x_0_z_0_xyy_xxyzz[k];

                g_x_0_z_0_xyyy_xyyy[k] = -g_x_0_z_0_xyy_xyyy[k] * ab_y + g_x_0_z_0_xyy_xyyyy[k];

                g_x_0_z_0_xyyy_xyyz[k] = -g_x_0_z_0_xyy_xyyz[k] * ab_y + g_x_0_z_0_xyy_xyyyz[k];

                g_x_0_z_0_xyyy_xyzz[k] = -g_x_0_z_0_xyy_xyzz[k] * ab_y + g_x_0_z_0_xyy_xyyzz[k];

                g_x_0_z_0_xyyy_xzzz[k] = -g_x_0_z_0_xyy_xzzz[k] * ab_y + g_x_0_z_0_xyy_xyzzz[k];

                g_x_0_z_0_xyyy_yyyy[k] = -g_x_0_z_0_xyy_yyyy[k] * ab_y + g_x_0_z_0_xyy_yyyyy[k];

                g_x_0_z_0_xyyy_yyyz[k] = -g_x_0_z_0_xyy_yyyz[k] * ab_y + g_x_0_z_0_xyy_yyyyz[k];

                g_x_0_z_0_xyyy_yyzz[k] = -g_x_0_z_0_xyy_yyzz[k] * ab_y + g_x_0_z_0_xyy_yyyzz[k];

                g_x_0_z_0_xyyy_yzzz[k] = -g_x_0_z_0_xyy_yzzz[k] * ab_y + g_x_0_z_0_xyy_yyzzz[k];

                g_x_0_z_0_xyyy_zzzz[k] = -g_x_0_z_0_xyy_zzzz[k] * ab_y + g_x_0_z_0_xyy_yzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyz_xxxx = cbuffer.data(gg_geom_1010_off + 555 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_xxxy = cbuffer.data(gg_geom_1010_off + 556 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_xxxz = cbuffer.data(gg_geom_1010_off + 557 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_xxyy = cbuffer.data(gg_geom_1010_off + 558 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_xxyz = cbuffer.data(gg_geom_1010_off + 559 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_xxzz = cbuffer.data(gg_geom_1010_off + 560 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_xyyy = cbuffer.data(gg_geom_1010_off + 561 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_xyyz = cbuffer.data(gg_geom_1010_off + 562 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_xyzz = cbuffer.data(gg_geom_1010_off + 563 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_xzzz = cbuffer.data(gg_geom_1010_off + 564 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_yyyy = cbuffer.data(gg_geom_1010_off + 565 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_yyyz = cbuffer.data(gg_geom_1010_off + 566 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_yyzz = cbuffer.data(gg_geom_1010_off + 567 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_yzzz = cbuffer.data(gg_geom_1010_off + 568 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_zzzz = cbuffer.data(gg_geom_1010_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyz_xxxx, g_x_0_z_0_xyyz_xxxy, g_x_0_z_0_xyyz_xxxz, g_x_0_z_0_xyyz_xxyy, g_x_0_z_0_xyyz_xxyz, g_x_0_z_0_xyyz_xxzz, g_x_0_z_0_xyyz_xyyy, g_x_0_z_0_xyyz_xyyz, g_x_0_z_0_xyyz_xyzz, g_x_0_z_0_xyyz_xzzz, g_x_0_z_0_xyyz_yyyy, g_x_0_z_0_xyyz_yyyz, g_x_0_z_0_xyyz_yyzz, g_x_0_z_0_xyyz_yzzz, g_x_0_z_0_xyyz_zzzz, g_x_0_z_0_xyz_xxxx, g_x_0_z_0_xyz_xxxxy, g_x_0_z_0_xyz_xxxy, g_x_0_z_0_xyz_xxxyy, g_x_0_z_0_xyz_xxxyz, g_x_0_z_0_xyz_xxxz, g_x_0_z_0_xyz_xxyy, g_x_0_z_0_xyz_xxyyy, g_x_0_z_0_xyz_xxyyz, g_x_0_z_0_xyz_xxyz, g_x_0_z_0_xyz_xxyzz, g_x_0_z_0_xyz_xxzz, g_x_0_z_0_xyz_xyyy, g_x_0_z_0_xyz_xyyyy, g_x_0_z_0_xyz_xyyyz, g_x_0_z_0_xyz_xyyz, g_x_0_z_0_xyz_xyyzz, g_x_0_z_0_xyz_xyzz, g_x_0_z_0_xyz_xyzzz, g_x_0_z_0_xyz_xzzz, g_x_0_z_0_xyz_yyyy, g_x_0_z_0_xyz_yyyyy, g_x_0_z_0_xyz_yyyyz, g_x_0_z_0_xyz_yyyz, g_x_0_z_0_xyz_yyyzz, g_x_0_z_0_xyz_yyzz, g_x_0_z_0_xyz_yyzzz, g_x_0_z_0_xyz_yzzz, g_x_0_z_0_xyz_yzzzz, g_x_0_z_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyz_xxxx[k] = -g_x_0_z_0_xyz_xxxx[k] * ab_y + g_x_0_z_0_xyz_xxxxy[k];

                g_x_0_z_0_xyyz_xxxy[k] = -g_x_0_z_0_xyz_xxxy[k] * ab_y + g_x_0_z_0_xyz_xxxyy[k];

                g_x_0_z_0_xyyz_xxxz[k] = -g_x_0_z_0_xyz_xxxz[k] * ab_y + g_x_0_z_0_xyz_xxxyz[k];

                g_x_0_z_0_xyyz_xxyy[k] = -g_x_0_z_0_xyz_xxyy[k] * ab_y + g_x_0_z_0_xyz_xxyyy[k];

                g_x_0_z_0_xyyz_xxyz[k] = -g_x_0_z_0_xyz_xxyz[k] * ab_y + g_x_0_z_0_xyz_xxyyz[k];

                g_x_0_z_0_xyyz_xxzz[k] = -g_x_0_z_0_xyz_xxzz[k] * ab_y + g_x_0_z_0_xyz_xxyzz[k];

                g_x_0_z_0_xyyz_xyyy[k] = -g_x_0_z_0_xyz_xyyy[k] * ab_y + g_x_0_z_0_xyz_xyyyy[k];

                g_x_0_z_0_xyyz_xyyz[k] = -g_x_0_z_0_xyz_xyyz[k] * ab_y + g_x_0_z_0_xyz_xyyyz[k];

                g_x_0_z_0_xyyz_xyzz[k] = -g_x_0_z_0_xyz_xyzz[k] * ab_y + g_x_0_z_0_xyz_xyyzz[k];

                g_x_0_z_0_xyyz_xzzz[k] = -g_x_0_z_0_xyz_xzzz[k] * ab_y + g_x_0_z_0_xyz_xyzzz[k];

                g_x_0_z_0_xyyz_yyyy[k] = -g_x_0_z_0_xyz_yyyy[k] * ab_y + g_x_0_z_0_xyz_yyyyy[k];

                g_x_0_z_0_xyyz_yyyz[k] = -g_x_0_z_0_xyz_yyyz[k] * ab_y + g_x_0_z_0_xyz_yyyyz[k];

                g_x_0_z_0_xyyz_yyzz[k] = -g_x_0_z_0_xyz_yyzz[k] * ab_y + g_x_0_z_0_xyz_yyyzz[k];

                g_x_0_z_0_xyyz_yzzz[k] = -g_x_0_z_0_xyz_yzzz[k] * ab_y + g_x_0_z_0_xyz_yyzzz[k];

                g_x_0_z_0_xyyz_zzzz[k] = -g_x_0_z_0_xyz_zzzz[k] * ab_y + g_x_0_z_0_xyz_yzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyzz_xxxx = cbuffer.data(gg_geom_1010_off + 570 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_xxxy = cbuffer.data(gg_geom_1010_off + 571 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_xxxz = cbuffer.data(gg_geom_1010_off + 572 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_xxyy = cbuffer.data(gg_geom_1010_off + 573 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_xxyz = cbuffer.data(gg_geom_1010_off + 574 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_xxzz = cbuffer.data(gg_geom_1010_off + 575 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_xyyy = cbuffer.data(gg_geom_1010_off + 576 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_xyyz = cbuffer.data(gg_geom_1010_off + 577 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_xyzz = cbuffer.data(gg_geom_1010_off + 578 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_xzzz = cbuffer.data(gg_geom_1010_off + 579 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_yyyy = cbuffer.data(gg_geom_1010_off + 580 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_yyyz = cbuffer.data(gg_geom_1010_off + 581 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_yyzz = cbuffer.data(gg_geom_1010_off + 582 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_yzzz = cbuffer.data(gg_geom_1010_off + 583 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_zzzz = cbuffer.data(gg_geom_1010_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyzz_xxxx, g_x_0_z_0_xyzz_xxxy, g_x_0_z_0_xyzz_xxxz, g_x_0_z_0_xyzz_xxyy, g_x_0_z_0_xyzz_xxyz, g_x_0_z_0_xyzz_xxzz, g_x_0_z_0_xyzz_xyyy, g_x_0_z_0_xyzz_xyyz, g_x_0_z_0_xyzz_xyzz, g_x_0_z_0_xyzz_xzzz, g_x_0_z_0_xyzz_yyyy, g_x_0_z_0_xyzz_yyyz, g_x_0_z_0_xyzz_yyzz, g_x_0_z_0_xyzz_yzzz, g_x_0_z_0_xyzz_zzzz, g_x_0_z_0_xzz_xxxx, g_x_0_z_0_xzz_xxxxy, g_x_0_z_0_xzz_xxxy, g_x_0_z_0_xzz_xxxyy, g_x_0_z_0_xzz_xxxyz, g_x_0_z_0_xzz_xxxz, g_x_0_z_0_xzz_xxyy, g_x_0_z_0_xzz_xxyyy, g_x_0_z_0_xzz_xxyyz, g_x_0_z_0_xzz_xxyz, g_x_0_z_0_xzz_xxyzz, g_x_0_z_0_xzz_xxzz, g_x_0_z_0_xzz_xyyy, g_x_0_z_0_xzz_xyyyy, g_x_0_z_0_xzz_xyyyz, g_x_0_z_0_xzz_xyyz, g_x_0_z_0_xzz_xyyzz, g_x_0_z_0_xzz_xyzz, g_x_0_z_0_xzz_xyzzz, g_x_0_z_0_xzz_xzzz, g_x_0_z_0_xzz_yyyy, g_x_0_z_0_xzz_yyyyy, g_x_0_z_0_xzz_yyyyz, g_x_0_z_0_xzz_yyyz, g_x_0_z_0_xzz_yyyzz, g_x_0_z_0_xzz_yyzz, g_x_0_z_0_xzz_yyzzz, g_x_0_z_0_xzz_yzzz, g_x_0_z_0_xzz_yzzzz, g_x_0_z_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyzz_xxxx[k] = -g_x_0_z_0_xzz_xxxx[k] * ab_y + g_x_0_z_0_xzz_xxxxy[k];

                g_x_0_z_0_xyzz_xxxy[k] = -g_x_0_z_0_xzz_xxxy[k] * ab_y + g_x_0_z_0_xzz_xxxyy[k];

                g_x_0_z_0_xyzz_xxxz[k] = -g_x_0_z_0_xzz_xxxz[k] * ab_y + g_x_0_z_0_xzz_xxxyz[k];

                g_x_0_z_0_xyzz_xxyy[k] = -g_x_0_z_0_xzz_xxyy[k] * ab_y + g_x_0_z_0_xzz_xxyyy[k];

                g_x_0_z_0_xyzz_xxyz[k] = -g_x_0_z_0_xzz_xxyz[k] * ab_y + g_x_0_z_0_xzz_xxyyz[k];

                g_x_0_z_0_xyzz_xxzz[k] = -g_x_0_z_0_xzz_xxzz[k] * ab_y + g_x_0_z_0_xzz_xxyzz[k];

                g_x_0_z_0_xyzz_xyyy[k] = -g_x_0_z_0_xzz_xyyy[k] * ab_y + g_x_0_z_0_xzz_xyyyy[k];

                g_x_0_z_0_xyzz_xyyz[k] = -g_x_0_z_0_xzz_xyyz[k] * ab_y + g_x_0_z_0_xzz_xyyyz[k];

                g_x_0_z_0_xyzz_xyzz[k] = -g_x_0_z_0_xzz_xyzz[k] * ab_y + g_x_0_z_0_xzz_xyyzz[k];

                g_x_0_z_0_xyzz_xzzz[k] = -g_x_0_z_0_xzz_xzzz[k] * ab_y + g_x_0_z_0_xzz_xyzzz[k];

                g_x_0_z_0_xyzz_yyyy[k] = -g_x_0_z_0_xzz_yyyy[k] * ab_y + g_x_0_z_0_xzz_yyyyy[k];

                g_x_0_z_0_xyzz_yyyz[k] = -g_x_0_z_0_xzz_yyyz[k] * ab_y + g_x_0_z_0_xzz_yyyyz[k];

                g_x_0_z_0_xyzz_yyzz[k] = -g_x_0_z_0_xzz_yyzz[k] * ab_y + g_x_0_z_0_xzz_yyyzz[k];

                g_x_0_z_0_xyzz_yzzz[k] = -g_x_0_z_0_xzz_yzzz[k] * ab_y + g_x_0_z_0_xzz_yyzzz[k];

                g_x_0_z_0_xyzz_zzzz[k] = -g_x_0_z_0_xzz_zzzz[k] * ab_y + g_x_0_z_0_xzz_yzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xzzz_xxxx = cbuffer.data(gg_geom_1010_off + 585 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_xxxy = cbuffer.data(gg_geom_1010_off + 586 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_xxxz = cbuffer.data(gg_geom_1010_off + 587 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_xxyy = cbuffer.data(gg_geom_1010_off + 588 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_xxyz = cbuffer.data(gg_geom_1010_off + 589 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_xxzz = cbuffer.data(gg_geom_1010_off + 590 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_xyyy = cbuffer.data(gg_geom_1010_off + 591 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_xyyz = cbuffer.data(gg_geom_1010_off + 592 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_xyzz = cbuffer.data(gg_geom_1010_off + 593 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_xzzz = cbuffer.data(gg_geom_1010_off + 594 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_yyyy = cbuffer.data(gg_geom_1010_off + 595 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_yyyz = cbuffer.data(gg_geom_1010_off + 596 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_yyzz = cbuffer.data(gg_geom_1010_off + 597 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_yzzz = cbuffer.data(gg_geom_1010_off + 598 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_zzzz = cbuffer.data(gg_geom_1010_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xzz_xxxx, g_x_0_z_0_xzz_xxxxz, g_x_0_z_0_xzz_xxxy, g_x_0_z_0_xzz_xxxyz, g_x_0_z_0_xzz_xxxz, g_x_0_z_0_xzz_xxxzz, g_x_0_z_0_xzz_xxyy, g_x_0_z_0_xzz_xxyyz, g_x_0_z_0_xzz_xxyz, g_x_0_z_0_xzz_xxyzz, g_x_0_z_0_xzz_xxzz, g_x_0_z_0_xzz_xxzzz, g_x_0_z_0_xzz_xyyy, g_x_0_z_0_xzz_xyyyz, g_x_0_z_0_xzz_xyyz, g_x_0_z_0_xzz_xyyzz, g_x_0_z_0_xzz_xyzz, g_x_0_z_0_xzz_xyzzz, g_x_0_z_0_xzz_xzzz, g_x_0_z_0_xzz_xzzzz, g_x_0_z_0_xzz_yyyy, g_x_0_z_0_xzz_yyyyz, g_x_0_z_0_xzz_yyyz, g_x_0_z_0_xzz_yyyzz, g_x_0_z_0_xzz_yyzz, g_x_0_z_0_xzz_yyzzz, g_x_0_z_0_xzz_yzzz, g_x_0_z_0_xzz_yzzzz, g_x_0_z_0_xzz_zzzz, g_x_0_z_0_xzz_zzzzz, g_x_0_z_0_xzzz_xxxx, g_x_0_z_0_xzzz_xxxy, g_x_0_z_0_xzzz_xxxz, g_x_0_z_0_xzzz_xxyy, g_x_0_z_0_xzzz_xxyz, g_x_0_z_0_xzzz_xxzz, g_x_0_z_0_xzzz_xyyy, g_x_0_z_0_xzzz_xyyz, g_x_0_z_0_xzzz_xyzz, g_x_0_z_0_xzzz_xzzz, g_x_0_z_0_xzzz_yyyy, g_x_0_z_0_xzzz_yyyz, g_x_0_z_0_xzzz_yyzz, g_x_0_z_0_xzzz_yzzz, g_x_0_z_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xzzz_xxxx[k] = -g_x_0_z_0_xzz_xxxx[k] * ab_z + g_x_0_z_0_xzz_xxxxz[k];

                g_x_0_z_0_xzzz_xxxy[k] = -g_x_0_z_0_xzz_xxxy[k] * ab_z + g_x_0_z_0_xzz_xxxyz[k];

                g_x_0_z_0_xzzz_xxxz[k] = -g_x_0_z_0_xzz_xxxz[k] * ab_z + g_x_0_z_0_xzz_xxxzz[k];

                g_x_0_z_0_xzzz_xxyy[k] = -g_x_0_z_0_xzz_xxyy[k] * ab_z + g_x_0_z_0_xzz_xxyyz[k];

                g_x_0_z_0_xzzz_xxyz[k] = -g_x_0_z_0_xzz_xxyz[k] * ab_z + g_x_0_z_0_xzz_xxyzz[k];

                g_x_0_z_0_xzzz_xxzz[k] = -g_x_0_z_0_xzz_xxzz[k] * ab_z + g_x_0_z_0_xzz_xxzzz[k];

                g_x_0_z_0_xzzz_xyyy[k] = -g_x_0_z_0_xzz_xyyy[k] * ab_z + g_x_0_z_0_xzz_xyyyz[k];

                g_x_0_z_0_xzzz_xyyz[k] = -g_x_0_z_0_xzz_xyyz[k] * ab_z + g_x_0_z_0_xzz_xyyzz[k];

                g_x_0_z_0_xzzz_xyzz[k] = -g_x_0_z_0_xzz_xyzz[k] * ab_z + g_x_0_z_0_xzz_xyzzz[k];

                g_x_0_z_0_xzzz_xzzz[k] = -g_x_0_z_0_xzz_xzzz[k] * ab_z + g_x_0_z_0_xzz_xzzzz[k];

                g_x_0_z_0_xzzz_yyyy[k] = -g_x_0_z_0_xzz_yyyy[k] * ab_z + g_x_0_z_0_xzz_yyyyz[k];

                g_x_0_z_0_xzzz_yyyz[k] = -g_x_0_z_0_xzz_yyyz[k] * ab_z + g_x_0_z_0_xzz_yyyzz[k];

                g_x_0_z_0_xzzz_yyzz[k] = -g_x_0_z_0_xzz_yyzz[k] * ab_z + g_x_0_z_0_xzz_yyzzz[k];

                g_x_0_z_0_xzzz_yzzz[k] = -g_x_0_z_0_xzz_yzzz[k] * ab_z + g_x_0_z_0_xzz_yzzzz[k];

                g_x_0_z_0_xzzz_zzzz[k] = -g_x_0_z_0_xzz_zzzz[k] * ab_z + g_x_0_z_0_xzz_zzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyy_xxxx = cbuffer.data(gg_geom_1010_off + 600 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_xxxy = cbuffer.data(gg_geom_1010_off + 601 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_xxxz = cbuffer.data(gg_geom_1010_off + 602 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_xxyy = cbuffer.data(gg_geom_1010_off + 603 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_xxyz = cbuffer.data(gg_geom_1010_off + 604 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_xxzz = cbuffer.data(gg_geom_1010_off + 605 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_xyyy = cbuffer.data(gg_geom_1010_off + 606 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_xyyz = cbuffer.data(gg_geom_1010_off + 607 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_xyzz = cbuffer.data(gg_geom_1010_off + 608 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_xzzz = cbuffer.data(gg_geom_1010_off + 609 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_yyyy = cbuffer.data(gg_geom_1010_off + 610 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_yyyz = cbuffer.data(gg_geom_1010_off + 611 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_yyzz = cbuffer.data(gg_geom_1010_off + 612 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_yzzz = cbuffer.data(gg_geom_1010_off + 613 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_zzzz = cbuffer.data(gg_geom_1010_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyy_xxxx, g_x_0_z_0_yyy_xxxxy, g_x_0_z_0_yyy_xxxy, g_x_0_z_0_yyy_xxxyy, g_x_0_z_0_yyy_xxxyz, g_x_0_z_0_yyy_xxxz, g_x_0_z_0_yyy_xxyy, g_x_0_z_0_yyy_xxyyy, g_x_0_z_0_yyy_xxyyz, g_x_0_z_0_yyy_xxyz, g_x_0_z_0_yyy_xxyzz, g_x_0_z_0_yyy_xxzz, g_x_0_z_0_yyy_xyyy, g_x_0_z_0_yyy_xyyyy, g_x_0_z_0_yyy_xyyyz, g_x_0_z_0_yyy_xyyz, g_x_0_z_0_yyy_xyyzz, g_x_0_z_0_yyy_xyzz, g_x_0_z_0_yyy_xyzzz, g_x_0_z_0_yyy_xzzz, g_x_0_z_0_yyy_yyyy, g_x_0_z_0_yyy_yyyyy, g_x_0_z_0_yyy_yyyyz, g_x_0_z_0_yyy_yyyz, g_x_0_z_0_yyy_yyyzz, g_x_0_z_0_yyy_yyzz, g_x_0_z_0_yyy_yyzzz, g_x_0_z_0_yyy_yzzz, g_x_0_z_0_yyy_yzzzz, g_x_0_z_0_yyy_zzzz, g_x_0_z_0_yyyy_xxxx, g_x_0_z_0_yyyy_xxxy, g_x_0_z_0_yyyy_xxxz, g_x_0_z_0_yyyy_xxyy, g_x_0_z_0_yyyy_xxyz, g_x_0_z_0_yyyy_xxzz, g_x_0_z_0_yyyy_xyyy, g_x_0_z_0_yyyy_xyyz, g_x_0_z_0_yyyy_xyzz, g_x_0_z_0_yyyy_xzzz, g_x_0_z_0_yyyy_yyyy, g_x_0_z_0_yyyy_yyyz, g_x_0_z_0_yyyy_yyzz, g_x_0_z_0_yyyy_yzzz, g_x_0_z_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyy_xxxx[k] = -g_x_0_z_0_yyy_xxxx[k] * ab_y + g_x_0_z_0_yyy_xxxxy[k];

                g_x_0_z_0_yyyy_xxxy[k] = -g_x_0_z_0_yyy_xxxy[k] * ab_y + g_x_0_z_0_yyy_xxxyy[k];

                g_x_0_z_0_yyyy_xxxz[k] = -g_x_0_z_0_yyy_xxxz[k] * ab_y + g_x_0_z_0_yyy_xxxyz[k];

                g_x_0_z_0_yyyy_xxyy[k] = -g_x_0_z_0_yyy_xxyy[k] * ab_y + g_x_0_z_0_yyy_xxyyy[k];

                g_x_0_z_0_yyyy_xxyz[k] = -g_x_0_z_0_yyy_xxyz[k] * ab_y + g_x_0_z_0_yyy_xxyyz[k];

                g_x_0_z_0_yyyy_xxzz[k] = -g_x_0_z_0_yyy_xxzz[k] * ab_y + g_x_0_z_0_yyy_xxyzz[k];

                g_x_0_z_0_yyyy_xyyy[k] = -g_x_0_z_0_yyy_xyyy[k] * ab_y + g_x_0_z_0_yyy_xyyyy[k];

                g_x_0_z_0_yyyy_xyyz[k] = -g_x_0_z_0_yyy_xyyz[k] * ab_y + g_x_0_z_0_yyy_xyyyz[k];

                g_x_0_z_0_yyyy_xyzz[k] = -g_x_0_z_0_yyy_xyzz[k] * ab_y + g_x_0_z_0_yyy_xyyzz[k];

                g_x_0_z_0_yyyy_xzzz[k] = -g_x_0_z_0_yyy_xzzz[k] * ab_y + g_x_0_z_0_yyy_xyzzz[k];

                g_x_0_z_0_yyyy_yyyy[k] = -g_x_0_z_0_yyy_yyyy[k] * ab_y + g_x_0_z_0_yyy_yyyyy[k];

                g_x_0_z_0_yyyy_yyyz[k] = -g_x_0_z_0_yyy_yyyz[k] * ab_y + g_x_0_z_0_yyy_yyyyz[k];

                g_x_0_z_0_yyyy_yyzz[k] = -g_x_0_z_0_yyy_yyzz[k] * ab_y + g_x_0_z_0_yyy_yyyzz[k];

                g_x_0_z_0_yyyy_yzzz[k] = -g_x_0_z_0_yyy_yzzz[k] * ab_y + g_x_0_z_0_yyy_yyzzz[k];

                g_x_0_z_0_yyyy_zzzz[k] = -g_x_0_z_0_yyy_zzzz[k] * ab_y + g_x_0_z_0_yyy_yzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyz_xxxx = cbuffer.data(gg_geom_1010_off + 615 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_xxxy = cbuffer.data(gg_geom_1010_off + 616 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_xxxz = cbuffer.data(gg_geom_1010_off + 617 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_xxyy = cbuffer.data(gg_geom_1010_off + 618 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_xxyz = cbuffer.data(gg_geom_1010_off + 619 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_xxzz = cbuffer.data(gg_geom_1010_off + 620 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_xyyy = cbuffer.data(gg_geom_1010_off + 621 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_xyyz = cbuffer.data(gg_geom_1010_off + 622 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_xyzz = cbuffer.data(gg_geom_1010_off + 623 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_xzzz = cbuffer.data(gg_geom_1010_off + 624 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_yyyy = cbuffer.data(gg_geom_1010_off + 625 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_yyyz = cbuffer.data(gg_geom_1010_off + 626 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_yyzz = cbuffer.data(gg_geom_1010_off + 627 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_yzzz = cbuffer.data(gg_geom_1010_off + 628 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_zzzz = cbuffer.data(gg_geom_1010_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyz_xxxx, g_x_0_z_0_yyyz_xxxy, g_x_0_z_0_yyyz_xxxz, g_x_0_z_0_yyyz_xxyy, g_x_0_z_0_yyyz_xxyz, g_x_0_z_0_yyyz_xxzz, g_x_0_z_0_yyyz_xyyy, g_x_0_z_0_yyyz_xyyz, g_x_0_z_0_yyyz_xyzz, g_x_0_z_0_yyyz_xzzz, g_x_0_z_0_yyyz_yyyy, g_x_0_z_0_yyyz_yyyz, g_x_0_z_0_yyyz_yyzz, g_x_0_z_0_yyyz_yzzz, g_x_0_z_0_yyyz_zzzz, g_x_0_z_0_yyz_xxxx, g_x_0_z_0_yyz_xxxxy, g_x_0_z_0_yyz_xxxy, g_x_0_z_0_yyz_xxxyy, g_x_0_z_0_yyz_xxxyz, g_x_0_z_0_yyz_xxxz, g_x_0_z_0_yyz_xxyy, g_x_0_z_0_yyz_xxyyy, g_x_0_z_0_yyz_xxyyz, g_x_0_z_0_yyz_xxyz, g_x_0_z_0_yyz_xxyzz, g_x_0_z_0_yyz_xxzz, g_x_0_z_0_yyz_xyyy, g_x_0_z_0_yyz_xyyyy, g_x_0_z_0_yyz_xyyyz, g_x_0_z_0_yyz_xyyz, g_x_0_z_0_yyz_xyyzz, g_x_0_z_0_yyz_xyzz, g_x_0_z_0_yyz_xyzzz, g_x_0_z_0_yyz_xzzz, g_x_0_z_0_yyz_yyyy, g_x_0_z_0_yyz_yyyyy, g_x_0_z_0_yyz_yyyyz, g_x_0_z_0_yyz_yyyz, g_x_0_z_0_yyz_yyyzz, g_x_0_z_0_yyz_yyzz, g_x_0_z_0_yyz_yyzzz, g_x_0_z_0_yyz_yzzz, g_x_0_z_0_yyz_yzzzz, g_x_0_z_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyz_xxxx[k] = -g_x_0_z_0_yyz_xxxx[k] * ab_y + g_x_0_z_0_yyz_xxxxy[k];

                g_x_0_z_0_yyyz_xxxy[k] = -g_x_0_z_0_yyz_xxxy[k] * ab_y + g_x_0_z_0_yyz_xxxyy[k];

                g_x_0_z_0_yyyz_xxxz[k] = -g_x_0_z_0_yyz_xxxz[k] * ab_y + g_x_0_z_0_yyz_xxxyz[k];

                g_x_0_z_0_yyyz_xxyy[k] = -g_x_0_z_0_yyz_xxyy[k] * ab_y + g_x_0_z_0_yyz_xxyyy[k];

                g_x_0_z_0_yyyz_xxyz[k] = -g_x_0_z_0_yyz_xxyz[k] * ab_y + g_x_0_z_0_yyz_xxyyz[k];

                g_x_0_z_0_yyyz_xxzz[k] = -g_x_0_z_0_yyz_xxzz[k] * ab_y + g_x_0_z_0_yyz_xxyzz[k];

                g_x_0_z_0_yyyz_xyyy[k] = -g_x_0_z_0_yyz_xyyy[k] * ab_y + g_x_0_z_0_yyz_xyyyy[k];

                g_x_0_z_0_yyyz_xyyz[k] = -g_x_0_z_0_yyz_xyyz[k] * ab_y + g_x_0_z_0_yyz_xyyyz[k];

                g_x_0_z_0_yyyz_xyzz[k] = -g_x_0_z_0_yyz_xyzz[k] * ab_y + g_x_0_z_0_yyz_xyyzz[k];

                g_x_0_z_0_yyyz_xzzz[k] = -g_x_0_z_0_yyz_xzzz[k] * ab_y + g_x_0_z_0_yyz_xyzzz[k];

                g_x_0_z_0_yyyz_yyyy[k] = -g_x_0_z_0_yyz_yyyy[k] * ab_y + g_x_0_z_0_yyz_yyyyy[k];

                g_x_0_z_0_yyyz_yyyz[k] = -g_x_0_z_0_yyz_yyyz[k] * ab_y + g_x_0_z_0_yyz_yyyyz[k];

                g_x_0_z_0_yyyz_yyzz[k] = -g_x_0_z_0_yyz_yyzz[k] * ab_y + g_x_0_z_0_yyz_yyyzz[k];

                g_x_0_z_0_yyyz_yzzz[k] = -g_x_0_z_0_yyz_yzzz[k] * ab_y + g_x_0_z_0_yyz_yyzzz[k];

                g_x_0_z_0_yyyz_zzzz[k] = -g_x_0_z_0_yyz_zzzz[k] * ab_y + g_x_0_z_0_yyz_yzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyzz_xxxx = cbuffer.data(gg_geom_1010_off + 630 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_xxxy = cbuffer.data(gg_geom_1010_off + 631 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_xxxz = cbuffer.data(gg_geom_1010_off + 632 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_xxyy = cbuffer.data(gg_geom_1010_off + 633 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_xxyz = cbuffer.data(gg_geom_1010_off + 634 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_xxzz = cbuffer.data(gg_geom_1010_off + 635 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_xyyy = cbuffer.data(gg_geom_1010_off + 636 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_xyyz = cbuffer.data(gg_geom_1010_off + 637 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_xyzz = cbuffer.data(gg_geom_1010_off + 638 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_xzzz = cbuffer.data(gg_geom_1010_off + 639 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_yyyy = cbuffer.data(gg_geom_1010_off + 640 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_yyyz = cbuffer.data(gg_geom_1010_off + 641 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_yyzz = cbuffer.data(gg_geom_1010_off + 642 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_yzzz = cbuffer.data(gg_geom_1010_off + 643 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_zzzz = cbuffer.data(gg_geom_1010_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyzz_xxxx, g_x_0_z_0_yyzz_xxxy, g_x_0_z_0_yyzz_xxxz, g_x_0_z_0_yyzz_xxyy, g_x_0_z_0_yyzz_xxyz, g_x_0_z_0_yyzz_xxzz, g_x_0_z_0_yyzz_xyyy, g_x_0_z_0_yyzz_xyyz, g_x_0_z_0_yyzz_xyzz, g_x_0_z_0_yyzz_xzzz, g_x_0_z_0_yyzz_yyyy, g_x_0_z_0_yyzz_yyyz, g_x_0_z_0_yyzz_yyzz, g_x_0_z_0_yyzz_yzzz, g_x_0_z_0_yyzz_zzzz, g_x_0_z_0_yzz_xxxx, g_x_0_z_0_yzz_xxxxy, g_x_0_z_0_yzz_xxxy, g_x_0_z_0_yzz_xxxyy, g_x_0_z_0_yzz_xxxyz, g_x_0_z_0_yzz_xxxz, g_x_0_z_0_yzz_xxyy, g_x_0_z_0_yzz_xxyyy, g_x_0_z_0_yzz_xxyyz, g_x_0_z_0_yzz_xxyz, g_x_0_z_0_yzz_xxyzz, g_x_0_z_0_yzz_xxzz, g_x_0_z_0_yzz_xyyy, g_x_0_z_0_yzz_xyyyy, g_x_0_z_0_yzz_xyyyz, g_x_0_z_0_yzz_xyyz, g_x_0_z_0_yzz_xyyzz, g_x_0_z_0_yzz_xyzz, g_x_0_z_0_yzz_xyzzz, g_x_0_z_0_yzz_xzzz, g_x_0_z_0_yzz_yyyy, g_x_0_z_0_yzz_yyyyy, g_x_0_z_0_yzz_yyyyz, g_x_0_z_0_yzz_yyyz, g_x_0_z_0_yzz_yyyzz, g_x_0_z_0_yzz_yyzz, g_x_0_z_0_yzz_yyzzz, g_x_0_z_0_yzz_yzzz, g_x_0_z_0_yzz_yzzzz, g_x_0_z_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyzz_xxxx[k] = -g_x_0_z_0_yzz_xxxx[k] * ab_y + g_x_0_z_0_yzz_xxxxy[k];

                g_x_0_z_0_yyzz_xxxy[k] = -g_x_0_z_0_yzz_xxxy[k] * ab_y + g_x_0_z_0_yzz_xxxyy[k];

                g_x_0_z_0_yyzz_xxxz[k] = -g_x_0_z_0_yzz_xxxz[k] * ab_y + g_x_0_z_0_yzz_xxxyz[k];

                g_x_0_z_0_yyzz_xxyy[k] = -g_x_0_z_0_yzz_xxyy[k] * ab_y + g_x_0_z_0_yzz_xxyyy[k];

                g_x_0_z_0_yyzz_xxyz[k] = -g_x_0_z_0_yzz_xxyz[k] * ab_y + g_x_0_z_0_yzz_xxyyz[k];

                g_x_0_z_0_yyzz_xxzz[k] = -g_x_0_z_0_yzz_xxzz[k] * ab_y + g_x_0_z_0_yzz_xxyzz[k];

                g_x_0_z_0_yyzz_xyyy[k] = -g_x_0_z_0_yzz_xyyy[k] * ab_y + g_x_0_z_0_yzz_xyyyy[k];

                g_x_0_z_0_yyzz_xyyz[k] = -g_x_0_z_0_yzz_xyyz[k] * ab_y + g_x_0_z_0_yzz_xyyyz[k];

                g_x_0_z_0_yyzz_xyzz[k] = -g_x_0_z_0_yzz_xyzz[k] * ab_y + g_x_0_z_0_yzz_xyyzz[k];

                g_x_0_z_0_yyzz_xzzz[k] = -g_x_0_z_0_yzz_xzzz[k] * ab_y + g_x_0_z_0_yzz_xyzzz[k];

                g_x_0_z_0_yyzz_yyyy[k] = -g_x_0_z_0_yzz_yyyy[k] * ab_y + g_x_0_z_0_yzz_yyyyy[k];

                g_x_0_z_0_yyzz_yyyz[k] = -g_x_0_z_0_yzz_yyyz[k] * ab_y + g_x_0_z_0_yzz_yyyyz[k];

                g_x_0_z_0_yyzz_yyzz[k] = -g_x_0_z_0_yzz_yyzz[k] * ab_y + g_x_0_z_0_yzz_yyyzz[k];

                g_x_0_z_0_yyzz_yzzz[k] = -g_x_0_z_0_yzz_yzzz[k] * ab_y + g_x_0_z_0_yzz_yyzzz[k];

                g_x_0_z_0_yyzz_zzzz[k] = -g_x_0_z_0_yzz_zzzz[k] * ab_y + g_x_0_z_0_yzz_yzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yzzz_xxxx = cbuffer.data(gg_geom_1010_off + 645 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_xxxy = cbuffer.data(gg_geom_1010_off + 646 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_xxxz = cbuffer.data(gg_geom_1010_off + 647 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_xxyy = cbuffer.data(gg_geom_1010_off + 648 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_xxyz = cbuffer.data(gg_geom_1010_off + 649 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_xxzz = cbuffer.data(gg_geom_1010_off + 650 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_xyyy = cbuffer.data(gg_geom_1010_off + 651 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_xyyz = cbuffer.data(gg_geom_1010_off + 652 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_xyzz = cbuffer.data(gg_geom_1010_off + 653 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_xzzz = cbuffer.data(gg_geom_1010_off + 654 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_yyyy = cbuffer.data(gg_geom_1010_off + 655 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_yyyz = cbuffer.data(gg_geom_1010_off + 656 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_yyzz = cbuffer.data(gg_geom_1010_off + 657 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_yzzz = cbuffer.data(gg_geom_1010_off + 658 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_zzzz = cbuffer.data(gg_geom_1010_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yzzz_xxxx, g_x_0_z_0_yzzz_xxxy, g_x_0_z_0_yzzz_xxxz, g_x_0_z_0_yzzz_xxyy, g_x_0_z_0_yzzz_xxyz, g_x_0_z_0_yzzz_xxzz, g_x_0_z_0_yzzz_xyyy, g_x_0_z_0_yzzz_xyyz, g_x_0_z_0_yzzz_xyzz, g_x_0_z_0_yzzz_xzzz, g_x_0_z_0_yzzz_yyyy, g_x_0_z_0_yzzz_yyyz, g_x_0_z_0_yzzz_yyzz, g_x_0_z_0_yzzz_yzzz, g_x_0_z_0_yzzz_zzzz, g_x_0_z_0_zzz_xxxx, g_x_0_z_0_zzz_xxxxy, g_x_0_z_0_zzz_xxxy, g_x_0_z_0_zzz_xxxyy, g_x_0_z_0_zzz_xxxyz, g_x_0_z_0_zzz_xxxz, g_x_0_z_0_zzz_xxyy, g_x_0_z_0_zzz_xxyyy, g_x_0_z_0_zzz_xxyyz, g_x_0_z_0_zzz_xxyz, g_x_0_z_0_zzz_xxyzz, g_x_0_z_0_zzz_xxzz, g_x_0_z_0_zzz_xyyy, g_x_0_z_0_zzz_xyyyy, g_x_0_z_0_zzz_xyyyz, g_x_0_z_0_zzz_xyyz, g_x_0_z_0_zzz_xyyzz, g_x_0_z_0_zzz_xyzz, g_x_0_z_0_zzz_xyzzz, g_x_0_z_0_zzz_xzzz, g_x_0_z_0_zzz_yyyy, g_x_0_z_0_zzz_yyyyy, g_x_0_z_0_zzz_yyyyz, g_x_0_z_0_zzz_yyyz, g_x_0_z_0_zzz_yyyzz, g_x_0_z_0_zzz_yyzz, g_x_0_z_0_zzz_yyzzz, g_x_0_z_0_zzz_yzzz, g_x_0_z_0_zzz_yzzzz, g_x_0_z_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yzzz_xxxx[k] = -g_x_0_z_0_zzz_xxxx[k] * ab_y + g_x_0_z_0_zzz_xxxxy[k];

                g_x_0_z_0_yzzz_xxxy[k] = -g_x_0_z_0_zzz_xxxy[k] * ab_y + g_x_0_z_0_zzz_xxxyy[k];

                g_x_0_z_0_yzzz_xxxz[k] = -g_x_0_z_0_zzz_xxxz[k] * ab_y + g_x_0_z_0_zzz_xxxyz[k];

                g_x_0_z_0_yzzz_xxyy[k] = -g_x_0_z_0_zzz_xxyy[k] * ab_y + g_x_0_z_0_zzz_xxyyy[k];

                g_x_0_z_0_yzzz_xxyz[k] = -g_x_0_z_0_zzz_xxyz[k] * ab_y + g_x_0_z_0_zzz_xxyyz[k];

                g_x_0_z_0_yzzz_xxzz[k] = -g_x_0_z_0_zzz_xxzz[k] * ab_y + g_x_0_z_0_zzz_xxyzz[k];

                g_x_0_z_0_yzzz_xyyy[k] = -g_x_0_z_0_zzz_xyyy[k] * ab_y + g_x_0_z_0_zzz_xyyyy[k];

                g_x_0_z_0_yzzz_xyyz[k] = -g_x_0_z_0_zzz_xyyz[k] * ab_y + g_x_0_z_0_zzz_xyyyz[k];

                g_x_0_z_0_yzzz_xyzz[k] = -g_x_0_z_0_zzz_xyzz[k] * ab_y + g_x_0_z_0_zzz_xyyzz[k];

                g_x_0_z_0_yzzz_xzzz[k] = -g_x_0_z_0_zzz_xzzz[k] * ab_y + g_x_0_z_0_zzz_xyzzz[k];

                g_x_0_z_0_yzzz_yyyy[k] = -g_x_0_z_0_zzz_yyyy[k] * ab_y + g_x_0_z_0_zzz_yyyyy[k];

                g_x_0_z_0_yzzz_yyyz[k] = -g_x_0_z_0_zzz_yyyz[k] * ab_y + g_x_0_z_0_zzz_yyyyz[k];

                g_x_0_z_0_yzzz_yyzz[k] = -g_x_0_z_0_zzz_yyzz[k] * ab_y + g_x_0_z_0_zzz_yyyzz[k];

                g_x_0_z_0_yzzz_yzzz[k] = -g_x_0_z_0_zzz_yzzz[k] * ab_y + g_x_0_z_0_zzz_yyzzz[k];

                g_x_0_z_0_yzzz_zzzz[k] = -g_x_0_z_0_zzz_zzzz[k] * ab_y + g_x_0_z_0_zzz_yzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_zzzz_xxxx = cbuffer.data(gg_geom_1010_off + 660 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_xxxy = cbuffer.data(gg_geom_1010_off + 661 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_xxxz = cbuffer.data(gg_geom_1010_off + 662 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_xxyy = cbuffer.data(gg_geom_1010_off + 663 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_xxyz = cbuffer.data(gg_geom_1010_off + 664 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_xxzz = cbuffer.data(gg_geom_1010_off + 665 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_xyyy = cbuffer.data(gg_geom_1010_off + 666 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_xyyz = cbuffer.data(gg_geom_1010_off + 667 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_xyzz = cbuffer.data(gg_geom_1010_off + 668 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_xzzz = cbuffer.data(gg_geom_1010_off + 669 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_yyyy = cbuffer.data(gg_geom_1010_off + 670 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_yyyz = cbuffer.data(gg_geom_1010_off + 671 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_yyzz = cbuffer.data(gg_geom_1010_off + 672 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_yzzz = cbuffer.data(gg_geom_1010_off + 673 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_zzzz = cbuffer.data(gg_geom_1010_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_zzz_xxxx, g_x_0_z_0_zzz_xxxxz, g_x_0_z_0_zzz_xxxy, g_x_0_z_0_zzz_xxxyz, g_x_0_z_0_zzz_xxxz, g_x_0_z_0_zzz_xxxzz, g_x_0_z_0_zzz_xxyy, g_x_0_z_0_zzz_xxyyz, g_x_0_z_0_zzz_xxyz, g_x_0_z_0_zzz_xxyzz, g_x_0_z_0_zzz_xxzz, g_x_0_z_0_zzz_xxzzz, g_x_0_z_0_zzz_xyyy, g_x_0_z_0_zzz_xyyyz, g_x_0_z_0_zzz_xyyz, g_x_0_z_0_zzz_xyyzz, g_x_0_z_0_zzz_xyzz, g_x_0_z_0_zzz_xyzzz, g_x_0_z_0_zzz_xzzz, g_x_0_z_0_zzz_xzzzz, g_x_0_z_0_zzz_yyyy, g_x_0_z_0_zzz_yyyyz, g_x_0_z_0_zzz_yyyz, g_x_0_z_0_zzz_yyyzz, g_x_0_z_0_zzz_yyzz, g_x_0_z_0_zzz_yyzzz, g_x_0_z_0_zzz_yzzz, g_x_0_z_0_zzz_yzzzz, g_x_0_z_0_zzz_zzzz, g_x_0_z_0_zzz_zzzzz, g_x_0_z_0_zzzz_xxxx, g_x_0_z_0_zzzz_xxxy, g_x_0_z_0_zzzz_xxxz, g_x_0_z_0_zzzz_xxyy, g_x_0_z_0_zzzz_xxyz, g_x_0_z_0_zzzz_xxzz, g_x_0_z_0_zzzz_xyyy, g_x_0_z_0_zzzz_xyyz, g_x_0_z_0_zzzz_xyzz, g_x_0_z_0_zzzz_xzzz, g_x_0_z_0_zzzz_yyyy, g_x_0_z_0_zzzz_yyyz, g_x_0_z_0_zzzz_yyzz, g_x_0_z_0_zzzz_yzzz, g_x_0_z_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_zzzz_xxxx[k] = -g_x_0_z_0_zzz_xxxx[k] * ab_z + g_x_0_z_0_zzz_xxxxz[k];

                g_x_0_z_0_zzzz_xxxy[k] = -g_x_0_z_0_zzz_xxxy[k] * ab_z + g_x_0_z_0_zzz_xxxyz[k];

                g_x_0_z_0_zzzz_xxxz[k] = -g_x_0_z_0_zzz_xxxz[k] * ab_z + g_x_0_z_0_zzz_xxxzz[k];

                g_x_0_z_0_zzzz_xxyy[k] = -g_x_0_z_0_zzz_xxyy[k] * ab_z + g_x_0_z_0_zzz_xxyyz[k];

                g_x_0_z_0_zzzz_xxyz[k] = -g_x_0_z_0_zzz_xxyz[k] * ab_z + g_x_0_z_0_zzz_xxyzz[k];

                g_x_0_z_0_zzzz_xxzz[k] = -g_x_0_z_0_zzz_xxzz[k] * ab_z + g_x_0_z_0_zzz_xxzzz[k];

                g_x_0_z_0_zzzz_xyyy[k] = -g_x_0_z_0_zzz_xyyy[k] * ab_z + g_x_0_z_0_zzz_xyyyz[k];

                g_x_0_z_0_zzzz_xyyz[k] = -g_x_0_z_0_zzz_xyyz[k] * ab_z + g_x_0_z_0_zzz_xyyzz[k];

                g_x_0_z_0_zzzz_xyzz[k] = -g_x_0_z_0_zzz_xyzz[k] * ab_z + g_x_0_z_0_zzz_xyzzz[k];

                g_x_0_z_0_zzzz_xzzz[k] = -g_x_0_z_0_zzz_xzzz[k] * ab_z + g_x_0_z_0_zzz_xzzzz[k];

                g_x_0_z_0_zzzz_yyyy[k] = -g_x_0_z_0_zzz_yyyy[k] * ab_z + g_x_0_z_0_zzz_yyyyz[k];

                g_x_0_z_0_zzzz_yyyz[k] = -g_x_0_z_0_zzz_yyyz[k] * ab_z + g_x_0_z_0_zzz_yyyzz[k];

                g_x_0_z_0_zzzz_yyzz[k] = -g_x_0_z_0_zzz_yyzz[k] * ab_z + g_x_0_z_0_zzz_yyzzz[k];

                g_x_0_z_0_zzzz_yzzz[k] = -g_x_0_z_0_zzz_yzzz[k] * ab_z + g_x_0_z_0_zzz_yzzzz[k];

                g_x_0_z_0_zzzz_zzzz[k] = -g_x_0_z_0_zzz_zzzz[k] * ab_z + g_x_0_z_0_zzz_zzzzz[k];
            }

            /// Set up 675-690 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxx_xxxx = cbuffer.data(gg_geom_1010_off + 675 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_xxxy = cbuffer.data(gg_geom_1010_off + 676 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_xxxz = cbuffer.data(gg_geom_1010_off + 677 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_xxyy = cbuffer.data(gg_geom_1010_off + 678 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_xxyz = cbuffer.data(gg_geom_1010_off + 679 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_xxzz = cbuffer.data(gg_geom_1010_off + 680 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_xyyy = cbuffer.data(gg_geom_1010_off + 681 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_xyyz = cbuffer.data(gg_geom_1010_off + 682 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_xyzz = cbuffer.data(gg_geom_1010_off + 683 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_xzzz = cbuffer.data(gg_geom_1010_off + 684 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_yyyy = cbuffer.data(gg_geom_1010_off + 685 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_yyyz = cbuffer.data(gg_geom_1010_off + 686 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_yyzz = cbuffer.data(gg_geom_1010_off + 687 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_yzzz = cbuffer.data(gg_geom_1010_off + 688 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_zzzz = cbuffer.data(gg_geom_1010_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxx_xxxx, g_y_0_x_0_xxx_xxxxx, g_y_0_x_0_xxx_xxxxy, g_y_0_x_0_xxx_xxxxz, g_y_0_x_0_xxx_xxxy, g_y_0_x_0_xxx_xxxyy, g_y_0_x_0_xxx_xxxyz, g_y_0_x_0_xxx_xxxz, g_y_0_x_0_xxx_xxxzz, g_y_0_x_0_xxx_xxyy, g_y_0_x_0_xxx_xxyyy, g_y_0_x_0_xxx_xxyyz, g_y_0_x_0_xxx_xxyz, g_y_0_x_0_xxx_xxyzz, g_y_0_x_0_xxx_xxzz, g_y_0_x_0_xxx_xxzzz, g_y_0_x_0_xxx_xyyy, g_y_0_x_0_xxx_xyyyy, g_y_0_x_0_xxx_xyyyz, g_y_0_x_0_xxx_xyyz, g_y_0_x_0_xxx_xyyzz, g_y_0_x_0_xxx_xyzz, g_y_0_x_0_xxx_xyzzz, g_y_0_x_0_xxx_xzzz, g_y_0_x_0_xxx_xzzzz, g_y_0_x_0_xxx_yyyy, g_y_0_x_0_xxx_yyyz, g_y_0_x_0_xxx_yyzz, g_y_0_x_0_xxx_yzzz, g_y_0_x_0_xxx_zzzz, g_y_0_x_0_xxxx_xxxx, g_y_0_x_0_xxxx_xxxy, g_y_0_x_0_xxxx_xxxz, g_y_0_x_0_xxxx_xxyy, g_y_0_x_0_xxxx_xxyz, g_y_0_x_0_xxxx_xxzz, g_y_0_x_0_xxxx_xyyy, g_y_0_x_0_xxxx_xyyz, g_y_0_x_0_xxxx_xyzz, g_y_0_x_0_xxxx_xzzz, g_y_0_x_0_xxxx_yyyy, g_y_0_x_0_xxxx_yyyz, g_y_0_x_0_xxxx_yyzz, g_y_0_x_0_xxxx_yzzz, g_y_0_x_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxx_xxxx[k] = -g_y_0_x_0_xxx_xxxx[k] * ab_x + g_y_0_x_0_xxx_xxxxx[k];

                g_y_0_x_0_xxxx_xxxy[k] = -g_y_0_x_0_xxx_xxxy[k] * ab_x + g_y_0_x_0_xxx_xxxxy[k];

                g_y_0_x_0_xxxx_xxxz[k] = -g_y_0_x_0_xxx_xxxz[k] * ab_x + g_y_0_x_0_xxx_xxxxz[k];

                g_y_0_x_0_xxxx_xxyy[k] = -g_y_0_x_0_xxx_xxyy[k] * ab_x + g_y_0_x_0_xxx_xxxyy[k];

                g_y_0_x_0_xxxx_xxyz[k] = -g_y_0_x_0_xxx_xxyz[k] * ab_x + g_y_0_x_0_xxx_xxxyz[k];

                g_y_0_x_0_xxxx_xxzz[k] = -g_y_0_x_0_xxx_xxzz[k] * ab_x + g_y_0_x_0_xxx_xxxzz[k];

                g_y_0_x_0_xxxx_xyyy[k] = -g_y_0_x_0_xxx_xyyy[k] * ab_x + g_y_0_x_0_xxx_xxyyy[k];

                g_y_0_x_0_xxxx_xyyz[k] = -g_y_0_x_0_xxx_xyyz[k] * ab_x + g_y_0_x_0_xxx_xxyyz[k];

                g_y_0_x_0_xxxx_xyzz[k] = -g_y_0_x_0_xxx_xyzz[k] * ab_x + g_y_0_x_0_xxx_xxyzz[k];

                g_y_0_x_0_xxxx_xzzz[k] = -g_y_0_x_0_xxx_xzzz[k] * ab_x + g_y_0_x_0_xxx_xxzzz[k];

                g_y_0_x_0_xxxx_yyyy[k] = -g_y_0_x_0_xxx_yyyy[k] * ab_x + g_y_0_x_0_xxx_xyyyy[k];

                g_y_0_x_0_xxxx_yyyz[k] = -g_y_0_x_0_xxx_yyyz[k] * ab_x + g_y_0_x_0_xxx_xyyyz[k];

                g_y_0_x_0_xxxx_yyzz[k] = -g_y_0_x_0_xxx_yyzz[k] * ab_x + g_y_0_x_0_xxx_xyyzz[k];

                g_y_0_x_0_xxxx_yzzz[k] = -g_y_0_x_0_xxx_yzzz[k] * ab_x + g_y_0_x_0_xxx_xyzzz[k];

                g_y_0_x_0_xxxx_zzzz[k] = -g_y_0_x_0_xxx_zzzz[k] * ab_x + g_y_0_x_0_xxx_xzzzz[k];
            }

            /// Set up 690-705 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxy_xxxx = cbuffer.data(gg_geom_1010_off + 690 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_xxxy = cbuffer.data(gg_geom_1010_off + 691 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_xxxz = cbuffer.data(gg_geom_1010_off + 692 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_xxyy = cbuffer.data(gg_geom_1010_off + 693 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_xxyz = cbuffer.data(gg_geom_1010_off + 694 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_xxzz = cbuffer.data(gg_geom_1010_off + 695 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_xyyy = cbuffer.data(gg_geom_1010_off + 696 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_xyyz = cbuffer.data(gg_geom_1010_off + 697 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_xyzz = cbuffer.data(gg_geom_1010_off + 698 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_xzzz = cbuffer.data(gg_geom_1010_off + 699 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_yyyy = cbuffer.data(gg_geom_1010_off + 700 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_yyyz = cbuffer.data(gg_geom_1010_off + 701 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_yyzz = cbuffer.data(gg_geom_1010_off + 702 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_yzzz = cbuffer.data(gg_geom_1010_off + 703 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_zzzz = cbuffer.data(gg_geom_1010_off + 704 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxy_xxxx, g_y_0_x_0_xxxy_xxxy, g_y_0_x_0_xxxy_xxxz, g_y_0_x_0_xxxy_xxyy, g_y_0_x_0_xxxy_xxyz, g_y_0_x_0_xxxy_xxzz, g_y_0_x_0_xxxy_xyyy, g_y_0_x_0_xxxy_xyyz, g_y_0_x_0_xxxy_xyzz, g_y_0_x_0_xxxy_xzzz, g_y_0_x_0_xxxy_yyyy, g_y_0_x_0_xxxy_yyyz, g_y_0_x_0_xxxy_yyzz, g_y_0_x_0_xxxy_yzzz, g_y_0_x_0_xxxy_zzzz, g_y_0_x_0_xxy_xxxx, g_y_0_x_0_xxy_xxxxx, g_y_0_x_0_xxy_xxxxy, g_y_0_x_0_xxy_xxxxz, g_y_0_x_0_xxy_xxxy, g_y_0_x_0_xxy_xxxyy, g_y_0_x_0_xxy_xxxyz, g_y_0_x_0_xxy_xxxz, g_y_0_x_0_xxy_xxxzz, g_y_0_x_0_xxy_xxyy, g_y_0_x_0_xxy_xxyyy, g_y_0_x_0_xxy_xxyyz, g_y_0_x_0_xxy_xxyz, g_y_0_x_0_xxy_xxyzz, g_y_0_x_0_xxy_xxzz, g_y_0_x_0_xxy_xxzzz, g_y_0_x_0_xxy_xyyy, g_y_0_x_0_xxy_xyyyy, g_y_0_x_0_xxy_xyyyz, g_y_0_x_0_xxy_xyyz, g_y_0_x_0_xxy_xyyzz, g_y_0_x_0_xxy_xyzz, g_y_0_x_0_xxy_xyzzz, g_y_0_x_0_xxy_xzzz, g_y_0_x_0_xxy_xzzzz, g_y_0_x_0_xxy_yyyy, g_y_0_x_0_xxy_yyyz, g_y_0_x_0_xxy_yyzz, g_y_0_x_0_xxy_yzzz, g_y_0_x_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxy_xxxx[k] = -g_y_0_x_0_xxy_xxxx[k] * ab_x + g_y_0_x_0_xxy_xxxxx[k];

                g_y_0_x_0_xxxy_xxxy[k] = -g_y_0_x_0_xxy_xxxy[k] * ab_x + g_y_0_x_0_xxy_xxxxy[k];

                g_y_0_x_0_xxxy_xxxz[k] = -g_y_0_x_0_xxy_xxxz[k] * ab_x + g_y_0_x_0_xxy_xxxxz[k];

                g_y_0_x_0_xxxy_xxyy[k] = -g_y_0_x_0_xxy_xxyy[k] * ab_x + g_y_0_x_0_xxy_xxxyy[k];

                g_y_0_x_0_xxxy_xxyz[k] = -g_y_0_x_0_xxy_xxyz[k] * ab_x + g_y_0_x_0_xxy_xxxyz[k];

                g_y_0_x_0_xxxy_xxzz[k] = -g_y_0_x_0_xxy_xxzz[k] * ab_x + g_y_0_x_0_xxy_xxxzz[k];

                g_y_0_x_0_xxxy_xyyy[k] = -g_y_0_x_0_xxy_xyyy[k] * ab_x + g_y_0_x_0_xxy_xxyyy[k];

                g_y_0_x_0_xxxy_xyyz[k] = -g_y_0_x_0_xxy_xyyz[k] * ab_x + g_y_0_x_0_xxy_xxyyz[k];

                g_y_0_x_0_xxxy_xyzz[k] = -g_y_0_x_0_xxy_xyzz[k] * ab_x + g_y_0_x_0_xxy_xxyzz[k];

                g_y_0_x_0_xxxy_xzzz[k] = -g_y_0_x_0_xxy_xzzz[k] * ab_x + g_y_0_x_0_xxy_xxzzz[k];

                g_y_0_x_0_xxxy_yyyy[k] = -g_y_0_x_0_xxy_yyyy[k] * ab_x + g_y_0_x_0_xxy_xyyyy[k];

                g_y_0_x_0_xxxy_yyyz[k] = -g_y_0_x_0_xxy_yyyz[k] * ab_x + g_y_0_x_0_xxy_xyyyz[k];

                g_y_0_x_0_xxxy_yyzz[k] = -g_y_0_x_0_xxy_yyzz[k] * ab_x + g_y_0_x_0_xxy_xyyzz[k];

                g_y_0_x_0_xxxy_yzzz[k] = -g_y_0_x_0_xxy_yzzz[k] * ab_x + g_y_0_x_0_xxy_xyzzz[k];

                g_y_0_x_0_xxxy_zzzz[k] = -g_y_0_x_0_xxy_zzzz[k] * ab_x + g_y_0_x_0_xxy_xzzzz[k];
            }

            /// Set up 705-720 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxz_xxxx = cbuffer.data(gg_geom_1010_off + 705 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_xxxy = cbuffer.data(gg_geom_1010_off + 706 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_xxxz = cbuffer.data(gg_geom_1010_off + 707 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_xxyy = cbuffer.data(gg_geom_1010_off + 708 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_xxyz = cbuffer.data(gg_geom_1010_off + 709 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_xxzz = cbuffer.data(gg_geom_1010_off + 710 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_xyyy = cbuffer.data(gg_geom_1010_off + 711 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_xyyz = cbuffer.data(gg_geom_1010_off + 712 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_xyzz = cbuffer.data(gg_geom_1010_off + 713 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_xzzz = cbuffer.data(gg_geom_1010_off + 714 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_yyyy = cbuffer.data(gg_geom_1010_off + 715 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_yyyz = cbuffer.data(gg_geom_1010_off + 716 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_yyzz = cbuffer.data(gg_geom_1010_off + 717 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_yzzz = cbuffer.data(gg_geom_1010_off + 718 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_zzzz = cbuffer.data(gg_geom_1010_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxz_xxxx, g_y_0_x_0_xxxz_xxxy, g_y_0_x_0_xxxz_xxxz, g_y_0_x_0_xxxz_xxyy, g_y_0_x_0_xxxz_xxyz, g_y_0_x_0_xxxz_xxzz, g_y_0_x_0_xxxz_xyyy, g_y_0_x_0_xxxz_xyyz, g_y_0_x_0_xxxz_xyzz, g_y_0_x_0_xxxz_xzzz, g_y_0_x_0_xxxz_yyyy, g_y_0_x_0_xxxz_yyyz, g_y_0_x_0_xxxz_yyzz, g_y_0_x_0_xxxz_yzzz, g_y_0_x_0_xxxz_zzzz, g_y_0_x_0_xxz_xxxx, g_y_0_x_0_xxz_xxxxx, g_y_0_x_0_xxz_xxxxy, g_y_0_x_0_xxz_xxxxz, g_y_0_x_0_xxz_xxxy, g_y_0_x_0_xxz_xxxyy, g_y_0_x_0_xxz_xxxyz, g_y_0_x_0_xxz_xxxz, g_y_0_x_0_xxz_xxxzz, g_y_0_x_0_xxz_xxyy, g_y_0_x_0_xxz_xxyyy, g_y_0_x_0_xxz_xxyyz, g_y_0_x_0_xxz_xxyz, g_y_0_x_0_xxz_xxyzz, g_y_0_x_0_xxz_xxzz, g_y_0_x_0_xxz_xxzzz, g_y_0_x_0_xxz_xyyy, g_y_0_x_0_xxz_xyyyy, g_y_0_x_0_xxz_xyyyz, g_y_0_x_0_xxz_xyyz, g_y_0_x_0_xxz_xyyzz, g_y_0_x_0_xxz_xyzz, g_y_0_x_0_xxz_xyzzz, g_y_0_x_0_xxz_xzzz, g_y_0_x_0_xxz_xzzzz, g_y_0_x_0_xxz_yyyy, g_y_0_x_0_xxz_yyyz, g_y_0_x_0_xxz_yyzz, g_y_0_x_0_xxz_yzzz, g_y_0_x_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxz_xxxx[k] = -g_y_0_x_0_xxz_xxxx[k] * ab_x + g_y_0_x_0_xxz_xxxxx[k];

                g_y_0_x_0_xxxz_xxxy[k] = -g_y_0_x_0_xxz_xxxy[k] * ab_x + g_y_0_x_0_xxz_xxxxy[k];

                g_y_0_x_0_xxxz_xxxz[k] = -g_y_0_x_0_xxz_xxxz[k] * ab_x + g_y_0_x_0_xxz_xxxxz[k];

                g_y_0_x_0_xxxz_xxyy[k] = -g_y_0_x_0_xxz_xxyy[k] * ab_x + g_y_0_x_0_xxz_xxxyy[k];

                g_y_0_x_0_xxxz_xxyz[k] = -g_y_0_x_0_xxz_xxyz[k] * ab_x + g_y_0_x_0_xxz_xxxyz[k];

                g_y_0_x_0_xxxz_xxzz[k] = -g_y_0_x_0_xxz_xxzz[k] * ab_x + g_y_0_x_0_xxz_xxxzz[k];

                g_y_0_x_0_xxxz_xyyy[k] = -g_y_0_x_0_xxz_xyyy[k] * ab_x + g_y_0_x_0_xxz_xxyyy[k];

                g_y_0_x_0_xxxz_xyyz[k] = -g_y_0_x_0_xxz_xyyz[k] * ab_x + g_y_0_x_0_xxz_xxyyz[k];

                g_y_0_x_0_xxxz_xyzz[k] = -g_y_0_x_0_xxz_xyzz[k] * ab_x + g_y_0_x_0_xxz_xxyzz[k];

                g_y_0_x_0_xxxz_xzzz[k] = -g_y_0_x_0_xxz_xzzz[k] * ab_x + g_y_0_x_0_xxz_xxzzz[k];

                g_y_0_x_0_xxxz_yyyy[k] = -g_y_0_x_0_xxz_yyyy[k] * ab_x + g_y_0_x_0_xxz_xyyyy[k];

                g_y_0_x_0_xxxz_yyyz[k] = -g_y_0_x_0_xxz_yyyz[k] * ab_x + g_y_0_x_0_xxz_xyyyz[k];

                g_y_0_x_0_xxxz_yyzz[k] = -g_y_0_x_0_xxz_yyzz[k] * ab_x + g_y_0_x_0_xxz_xyyzz[k];

                g_y_0_x_0_xxxz_yzzz[k] = -g_y_0_x_0_xxz_yzzz[k] * ab_x + g_y_0_x_0_xxz_xyzzz[k];

                g_y_0_x_0_xxxz_zzzz[k] = -g_y_0_x_0_xxz_zzzz[k] * ab_x + g_y_0_x_0_xxz_xzzzz[k];
            }

            /// Set up 720-735 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyy_xxxx = cbuffer.data(gg_geom_1010_off + 720 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_xxxy = cbuffer.data(gg_geom_1010_off + 721 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_xxxz = cbuffer.data(gg_geom_1010_off + 722 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_xxyy = cbuffer.data(gg_geom_1010_off + 723 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_xxyz = cbuffer.data(gg_geom_1010_off + 724 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_xxzz = cbuffer.data(gg_geom_1010_off + 725 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_xyyy = cbuffer.data(gg_geom_1010_off + 726 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_xyyz = cbuffer.data(gg_geom_1010_off + 727 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_xyzz = cbuffer.data(gg_geom_1010_off + 728 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_xzzz = cbuffer.data(gg_geom_1010_off + 729 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_yyyy = cbuffer.data(gg_geom_1010_off + 730 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_yyyz = cbuffer.data(gg_geom_1010_off + 731 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_yyzz = cbuffer.data(gg_geom_1010_off + 732 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_yzzz = cbuffer.data(gg_geom_1010_off + 733 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_zzzz = cbuffer.data(gg_geom_1010_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyy_xxxx, g_y_0_x_0_xxyy_xxxy, g_y_0_x_0_xxyy_xxxz, g_y_0_x_0_xxyy_xxyy, g_y_0_x_0_xxyy_xxyz, g_y_0_x_0_xxyy_xxzz, g_y_0_x_0_xxyy_xyyy, g_y_0_x_0_xxyy_xyyz, g_y_0_x_0_xxyy_xyzz, g_y_0_x_0_xxyy_xzzz, g_y_0_x_0_xxyy_yyyy, g_y_0_x_0_xxyy_yyyz, g_y_0_x_0_xxyy_yyzz, g_y_0_x_0_xxyy_yzzz, g_y_0_x_0_xxyy_zzzz, g_y_0_x_0_xyy_xxxx, g_y_0_x_0_xyy_xxxxx, g_y_0_x_0_xyy_xxxxy, g_y_0_x_0_xyy_xxxxz, g_y_0_x_0_xyy_xxxy, g_y_0_x_0_xyy_xxxyy, g_y_0_x_0_xyy_xxxyz, g_y_0_x_0_xyy_xxxz, g_y_0_x_0_xyy_xxxzz, g_y_0_x_0_xyy_xxyy, g_y_0_x_0_xyy_xxyyy, g_y_0_x_0_xyy_xxyyz, g_y_0_x_0_xyy_xxyz, g_y_0_x_0_xyy_xxyzz, g_y_0_x_0_xyy_xxzz, g_y_0_x_0_xyy_xxzzz, g_y_0_x_0_xyy_xyyy, g_y_0_x_0_xyy_xyyyy, g_y_0_x_0_xyy_xyyyz, g_y_0_x_0_xyy_xyyz, g_y_0_x_0_xyy_xyyzz, g_y_0_x_0_xyy_xyzz, g_y_0_x_0_xyy_xyzzz, g_y_0_x_0_xyy_xzzz, g_y_0_x_0_xyy_xzzzz, g_y_0_x_0_xyy_yyyy, g_y_0_x_0_xyy_yyyz, g_y_0_x_0_xyy_yyzz, g_y_0_x_0_xyy_yzzz, g_y_0_x_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyy_xxxx[k] = -g_y_0_x_0_xyy_xxxx[k] * ab_x + g_y_0_x_0_xyy_xxxxx[k];

                g_y_0_x_0_xxyy_xxxy[k] = -g_y_0_x_0_xyy_xxxy[k] * ab_x + g_y_0_x_0_xyy_xxxxy[k];

                g_y_0_x_0_xxyy_xxxz[k] = -g_y_0_x_0_xyy_xxxz[k] * ab_x + g_y_0_x_0_xyy_xxxxz[k];

                g_y_0_x_0_xxyy_xxyy[k] = -g_y_0_x_0_xyy_xxyy[k] * ab_x + g_y_0_x_0_xyy_xxxyy[k];

                g_y_0_x_0_xxyy_xxyz[k] = -g_y_0_x_0_xyy_xxyz[k] * ab_x + g_y_0_x_0_xyy_xxxyz[k];

                g_y_0_x_0_xxyy_xxzz[k] = -g_y_0_x_0_xyy_xxzz[k] * ab_x + g_y_0_x_0_xyy_xxxzz[k];

                g_y_0_x_0_xxyy_xyyy[k] = -g_y_0_x_0_xyy_xyyy[k] * ab_x + g_y_0_x_0_xyy_xxyyy[k];

                g_y_0_x_0_xxyy_xyyz[k] = -g_y_0_x_0_xyy_xyyz[k] * ab_x + g_y_0_x_0_xyy_xxyyz[k];

                g_y_0_x_0_xxyy_xyzz[k] = -g_y_0_x_0_xyy_xyzz[k] * ab_x + g_y_0_x_0_xyy_xxyzz[k];

                g_y_0_x_0_xxyy_xzzz[k] = -g_y_0_x_0_xyy_xzzz[k] * ab_x + g_y_0_x_0_xyy_xxzzz[k];

                g_y_0_x_0_xxyy_yyyy[k] = -g_y_0_x_0_xyy_yyyy[k] * ab_x + g_y_0_x_0_xyy_xyyyy[k];

                g_y_0_x_0_xxyy_yyyz[k] = -g_y_0_x_0_xyy_yyyz[k] * ab_x + g_y_0_x_0_xyy_xyyyz[k];

                g_y_0_x_0_xxyy_yyzz[k] = -g_y_0_x_0_xyy_yyzz[k] * ab_x + g_y_0_x_0_xyy_xyyzz[k];

                g_y_0_x_0_xxyy_yzzz[k] = -g_y_0_x_0_xyy_yzzz[k] * ab_x + g_y_0_x_0_xyy_xyzzz[k];

                g_y_0_x_0_xxyy_zzzz[k] = -g_y_0_x_0_xyy_zzzz[k] * ab_x + g_y_0_x_0_xyy_xzzzz[k];
            }

            /// Set up 735-750 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyz_xxxx = cbuffer.data(gg_geom_1010_off + 735 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_xxxy = cbuffer.data(gg_geom_1010_off + 736 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_xxxz = cbuffer.data(gg_geom_1010_off + 737 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_xxyy = cbuffer.data(gg_geom_1010_off + 738 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_xxyz = cbuffer.data(gg_geom_1010_off + 739 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_xxzz = cbuffer.data(gg_geom_1010_off + 740 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_xyyy = cbuffer.data(gg_geom_1010_off + 741 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_xyyz = cbuffer.data(gg_geom_1010_off + 742 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_xyzz = cbuffer.data(gg_geom_1010_off + 743 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_xzzz = cbuffer.data(gg_geom_1010_off + 744 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_yyyy = cbuffer.data(gg_geom_1010_off + 745 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_yyyz = cbuffer.data(gg_geom_1010_off + 746 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_yyzz = cbuffer.data(gg_geom_1010_off + 747 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_yzzz = cbuffer.data(gg_geom_1010_off + 748 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_zzzz = cbuffer.data(gg_geom_1010_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyz_xxxx, g_y_0_x_0_xxyz_xxxy, g_y_0_x_0_xxyz_xxxz, g_y_0_x_0_xxyz_xxyy, g_y_0_x_0_xxyz_xxyz, g_y_0_x_0_xxyz_xxzz, g_y_0_x_0_xxyz_xyyy, g_y_0_x_0_xxyz_xyyz, g_y_0_x_0_xxyz_xyzz, g_y_0_x_0_xxyz_xzzz, g_y_0_x_0_xxyz_yyyy, g_y_0_x_0_xxyz_yyyz, g_y_0_x_0_xxyz_yyzz, g_y_0_x_0_xxyz_yzzz, g_y_0_x_0_xxyz_zzzz, g_y_0_x_0_xyz_xxxx, g_y_0_x_0_xyz_xxxxx, g_y_0_x_0_xyz_xxxxy, g_y_0_x_0_xyz_xxxxz, g_y_0_x_0_xyz_xxxy, g_y_0_x_0_xyz_xxxyy, g_y_0_x_0_xyz_xxxyz, g_y_0_x_0_xyz_xxxz, g_y_0_x_0_xyz_xxxzz, g_y_0_x_0_xyz_xxyy, g_y_0_x_0_xyz_xxyyy, g_y_0_x_0_xyz_xxyyz, g_y_0_x_0_xyz_xxyz, g_y_0_x_0_xyz_xxyzz, g_y_0_x_0_xyz_xxzz, g_y_0_x_0_xyz_xxzzz, g_y_0_x_0_xyz_xyyy, g_y_0_x_0_xyz_xyyyy, g_y_0_x_0_xyz_xyyyz, g_y_0_x_0_xyz_xyyz, g_y_0_x_0_xyz_xyyzz, g_y_0_x_0_xyz_xyzz, g_y_0_x_0_xyz_xyzzz, g_y_0_x_0_xyz_xzzz, g_y_0_x_0_xyz_xzzzz, g_y_0_x_0_xyz_yyyy, g_y_0_x_0_xyz_yyyz, g_y_0_x_0_xyz_yyzz, g_y_0_x_0_xyz_yzzz, g_y_0_x_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyz_xxxx[k] = -g_y_0_x_0_xyz_xxxx[k] * ab_x + g_y_0_x_0_xyz_xxxxx[k];

                g_y_0_x_0_xxyz_xxxy[k] = -g_y_0_x_0_xyz_xxxy[k] * ab_x + g_y_0_x_0_xyz_xxxxy[k];

                g_y_0_x_0_xxyz_xxxz[k] = -g_y_0_x_0_xyz_xxxz[k] * ab_x + g_y_0_x_0_xyz_xxxxz[k];

                g_y_0_x_0_xxyz_xxyy[k] = -g_y_0_x_0_xyz_xxyy[k] * ab_x + g_y_0_x_0_xyz_xxxyy[k];

                g_y_0_x_0_xxyz_xxyz[k] = -g_y_0_x_0_xyz_xxyz[k] * ab_x + g_y_0_x_0_xyz_xxxyz[k];

                g_y_0_x_0_xxyz_xxzz[k] = -g_y_0_x_0_xyz_xxzz[k] * ab_x + g_y_0_x_0_xyz_xxxzz[k];

                g_y_0_x_0_xxyz_xyyy[k] = -g_y_0_x_0_xyz_xyyy[k] * ab_x + g_y_0_x_0_xyz_xxyyy[k];

                g_y_0_x_0_xxyz_xyyz[k] = -g_y_0_x_0_xyz_xyyz[k] * ab_x + g_y_0_x_0_xyz_xxyyz[k];

                g_y_0_x_0_xxyz_xyzz[k] = -g_y_0_x_0_xyz_xyzz[k] * ab_x + g_y_0_x_0_xyz_xxyzz[k];

                g_y_0_x_0_xxyz_xzzz[k] = -g_y_0_x_0_xyz_xzzz[k] * ab_x + g_y_0_x_0_xyz_xxzzz[k];

                g_y_0_x_0_xxyz_yyyy[k] = -g_y_0_x_0_xyz_yyyy[k] * ab_x + g_y_0_x_0_xyz_xyyyy[k];

                g_y_0_x_0_xxyz_yyyz[k] = -g_y_0_x_0_xyz_yyyz[k] * ab_x + g_y_0_x_0_xyz_xyyyz[k];

                g_y_0_x_0_xxyz_yyzz[k] = -g_y_0_x_0_xyz_yyzz[k] * ab_x + g_y_0_x_0_xyz_xyyzz[k];

                g_y_0_x_0_xxyz_yzzz[k] = -g_y_0_x_0_xyz_yzzz[k] * ab_x + g_y_0_x_0_xyz_xyzzz[k];

                g_y_0_x_0_xxyz_zzzz[k] = -g_y_0_x_0_xyz_zzzz[k] * ab_x + g_y_0_x_0_xyz_xzzzz[k];
            }

            /// Set up 750-765 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxzz_xxxx = cbuffer.data(gg_geom_1010_off + 750 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_xxxy = cbuffer.data(gg_geom_1010_off + 751 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_xxxz = cbuffer.data(gg_geom_1010_off + 752 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_xxyy = cbuffer.data(gg_geom_1010_off + 753 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_xxyz = cbuffer.data(gg_geom_1010_off + 754 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_xxzz = cbuffer.data(gg_geom_1010_off + 755 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_xyyy = cbuffer.data(gg_geom_1010_off + 756 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_xyyz = cbuffer.data(gg_geom_1010_off + 757 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_xyzz = cbuffer.data(gg_geom_1010_off + 758 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_xzzz = cbuffer.data(gg_geom_1010_off + 759 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_yyyy = cbuffer.data(gg_geom_1010_off + 760 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_yyyz = cbuffer.data(gg_geom_1010_off + 761 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_yyzz = cbuffer.data(gg_geom_1010_off + 762 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_yzzz = cbuffer.data(gg_geom_1010_off + 763 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_zzzz = cbuffer.data(gg_geom_1010_off + 764 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxzz_xxxx, g_y_0_x_0_xxzz_xxxy, g_y_0_x_0_xxzz_xxxz, g_y_0_x_0_xxzz_xxyy, g_y_0_x_0_xxzz_xxyz, g_y_0_x_0_xxzz_xxzz, g_y_0_x_0_xxzz_xyyy, g_y_0_x_0_xxzz_xyyz, g_y_0_x_0_xxzz_xyzz, g_y_0_x_0_xxzz_xzzz, g_y_0_x_0_xxzz_yyyy, g_y_0_x_0_xxzz_yyyz, g_y_0_x_0_xxzz_yyzz, g_y_0_x_0_xxzz_yzzz, g_y_0_x_0_xxzz_zzzz, g_y_0_x_0_xzz_xxxx, g_y_0_x_0_xzz_xxxxx, g_y_0_x_0_xzz_xxxxy, g_y_0_x_0_xzz_xxxxz, g_y_0_x_0_xzz_xxxy, g_y_0_x_0_xzz_xxxyy, g_y_0_x_0_xzz_xxxyz, g_y_0_x_0_xzz_xxxz, g_y_0_x_0_xzz_xxxzz, g_y_0_x_0_xzz_xxyy, g_y_0_x_0_xzz_xxyyy, g_y_0_x_0_xzz_xxyyz, g_y_0_x_0_xzz_xxyz, g_y_0_x_0_xzz_xxyzz, g_y_0_x_0_xzz_xxzz, g_y_0_x_0_xzz_xxzzz, g_y_0_x_0_xzz_xyyy, g_y_0_x_0_xzz_xyyyy, g_y_0_x_0_xzz_xyyyz, g_y_0_x_0_xzz_xyyz, g_y_0_x_0_xzz_xyyzz, g_y_0_x_0_xzz_xyzz, g_y_0_x_0_xzz_xyzzz, g_y_0_x_0_xzz_xzzz, g_y_0_x_0_xzz_xzzzz, g_y_0_x_0_xzz_yyyy, g_y_0_x_0_xzz_yyyz, g_y_0_x_0_xzz_yyzz, g_y_0_x_0_xzz_yzzz, g_y_0_x_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxzz_xxxx[k] = -g_y_0_x_0_xzz_xxxx[k] * ab_x + g_y_0_x_0_xzz_xxxxx[k];

                g_y_0_x_0_xxzz_xxxy[k] = -g_y_0_x_0_xzz_xxxy[k] * ab_x + g_y_0_x_0_xzz_xxxxy[k];

                g_y_0_x_0_xxzz_xxxz[k] = -g_y_0_x_0_xzz_xxxz[k] * ab_x + g_y_0_x_0_xzz_xxxxz[k];

                g_y_0_x_0_xxzz_xxyy[k] = -g_y_0_x_0_xzz_xxyy[k] * ab_x + g_y_0_x_0_xzz_xxxyy[k];

                g_y_0_x_0_xxzz_xxyz[k] = -g_y_0_x_0_xzz_xxyz[k] * ab_x + g_y_0_x_0_xzz_xxxyz[k];

                g_y_0_x_0_xxzz_xxzz[k] = -g_y_0_x_0_xzz_xxzz[k] * ab_x + g_y_0_x_0_xzz_xxxzz[k];

                g_y_0_x_0_xxzz_xyyy[k] = -g_y_0_x_0_xzz_xyyy[k] * ab_x + g_y_0_x_0_xzz_xxyyy[k];

                g_y_0_x_0_xxzz_xyyz[k] = -g_y_0_x_0_xzz_xyyz[k] * ab_x + g_y_0_x_0_xzz_xxyyz[k];

                g_y_0_x_0_xxzz_xyzz[k] = -g_y_0_x_0_xzz_xyzz[k] * ab_x + g_y_0_x_0_xzz_xxyzz[k];

                g_y_0_x_0_xxzz_xzzz[k] = -g_y_0_x_0_xzz_xzzz[k] * ab_x + g_y_0_x_0_xzz_xxzzz[k];

                g_y_0_x_0_xxzz_yyyy[k] = -g_y_0_x_0_xzz_yyyy[k] * ab_x + g_y_0_x_0_xzz_xyyyy[k];

                g_y_0_x_0_xxzz_yyyz[k] = -g_y_0_x_0_xzz_yyyz[k] * ab_x + g_y_0_x_0_xzz_xyyyz[k];

                g_y_0_x_0_xxzz_yyzz[k] = -g_y_0_x_0_xzz_yyzz[k] * ab_x + g_y_0_x_0_xzz_xyyzz[k];

                g_y_0_x_0_xxzz_yzzz[k] = -g_y_0_x_0_xzz_yzzz[k] * ab_x + g_y_0_x_0_xzz_xyzzz[k];

                g_y_0_x_0_xxzz_zzzz[k] = -g_y_0_x_0_xzz_zzzz[k] * ab_x + g_y_0_x_0_xzz_xzzzz[k];
            }

            /// Set up 765-780 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyy_xxxx = cbuffer.data(gg_geom_1010_off + 765 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_xxxy = cbuffer.data(gg_geom_1010_off + 766 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_xxxz = cbuffer.data(gg_geom_1010_off + 767 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_xxyy = cbuffer.data(gg_geom_1010_off + 768 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_xxyz = cbuffer.data(gg_geom_1010_off + 769 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_xxzz = cbuffer.data(gg_geom_1010_off + 770 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_xyyy = cbuffer.data(gg_geom_1010_off + 771 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_xyyz = cbuffer.data(gg_geom_1010_off + 772 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_xyzz = cbuffer.data(gg_geom_1010_off + 773 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_xzzz = cbuffer.data(gg_geom_1010_off + 774 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_yyyy = cbuffer.data(gg_geom_1010_off + 775 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_yyyz = cbuffer.data(gg_geom_1010_off + 776 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_yyzz = cbuffer.data(gg_geom_1010_off + 777 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_yzzz = cbuffer.data(gg_geom_1010_off + 778 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_zzzz = cbuffer.data(gg_geom_1010_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyy_xxxx, g_y_0_x_0_xyyy_xxxy, g_y_0_x_0_xyyy_xxxz, g_y_0_x_0_xyyy_xxyy, g_y_0_x_0_xyyy_xxyz, g_y_0_x_0_xyyy_xxzz, g_y_0_x_0_xyyy_xyyy, g_y_0_x_0_xyyy_xyyz, g_y_0_x_0_xyyy_xyzz, g_y_0_x_0_xyyy_xzzz, g_y_0_x_0_xyyy_yyyy, g_y_0_x_0_xyyy_yyyz, g_y_0_x_0_xyyy_yyzz, g_y_0_x_0_xyyy_yzzz, g_y_0_x_0_xyyy_zzzz, g_y_0_x_0_yyy_xxxx, g_y_0_x_0_yyy_xxxxx, g_y_0_x_0_yyy_xxxxy, g_y_0_x_0_yyy_xxxxz, g_y_0_x_0_yyy_xxxy, g_y_0_x_0_yyy_xxxyy, g_y_0_x_0_yyy_xxxyz, g_y_0_x_0_yyy_xxxz, g_y_0_x_0_yyy_xxxzz, g_y_0_x_0_yyy_xxyy, g_y_0_x_0_yyy_xxyyy, g_y_0_x_0_yyy_xxyyz, g_y_0_x_0_yyy_xxyz, g_y_0_x_0_yyy_xxyzz, g_y_0_x_0_yyy_xxzz, g_y_0_x_0_yyy_xxzzz, g_y_0_x_0_yyy_xyyy, g_y_0_x_0_yyy_xyyyy, g_y_0_x_0_yyy_xyyyz, g_y_0_x_0_yyy_xyyz, g_y_0_x_0_yyy_xyyzz, g_y_0_x_0_yyy_xyzz, g_y_0_x_0_yyy_xyzzz, g_y_0_x_0_yyy_xzzz, g_y_0_x_0_yyy_xzzzz, g_y_0_x_0_yyy_yyyy, g_y_0_x_0_yyy_yyyz, g_y_0_x_0_yyy_yyzz, g_y_0_x_0_yyy_yzzz, g_y_0_x_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyy_xxxx[k] = -g_y_0_x_0_yyy_xxxx[k] * ab_x + g_y_0_x_0_yyy_xxxxx[k];

                g_y_0_x_0_xyyy_xxxy[k] = -g_y_0_x_0_yyy_xxxy[k] * ab_x + g_y_0_x_0_yyy_xxxxy[k];

                g_y_0_x_0_xyyy_xxxz[k] = -g_y_0_x_0_yyy_xxxz[k] * ab_x + g_y_0_x_0_yyy_xxxxz[k];

                g_y_0_x_0_xyyy_xxyy[k] = -g_y_0_x_0_yyy_xxyy[k] * ab_x + g_y_0_x_0_yyy_xxxyy[k];

                g_y_0_x_0_xyyy_xxyz[k] = -g_y_0_x_0_yyy_xxyz[k] * ab_x + g_y_0_x_0_yyy_xxxyz[k];

                g_y_0_x_0_xyyy_xxzz[k] = -g_y_0_x_0_yyy_xxzz[k] * ab_x + g_y_0_x_0_yyy_xxxzz[k];

                g_y_0_x_0_xyyy_xyyy[k] = -g_y_0_x_0_yyy_xyyy[k] * ab_x + g_y_0_x_0_yyy_xxyyy[k];

                g_y_0_x_0_xyyy_xyyz[k] = -g_y_0_x_0_yyy_xyyz[k] * ab_x + g_y_0_x_0_yyy_xxyyz[k];

                g_y_0_x_0_xyyy_xyzz[k] = -g_y_0_x_0_yyy_xyzz[k] * ab_x + g_y_0_x_0_yyy_xxyzz[k];

                g_y_0_x_0_xyyy_xzzz[k] = -g_y_0_x_0_yyy_xzzz[k] * ab_x + g_y_0_x_0_yyy_xxzzz[k];

                g_y_0_x_0_xyyy_yyyy[k] = -g_y_0_x_0_yyy_yyyy[k] * ab_x + g_y_0_x_0_yyy_xyyyy[k];

                g_y_0_x_0_xyyy_yyyz[k] = -g_y_0_x_0_yyy_yyyz[k] * ab_x + g_y_0_x_0_yyy_xyyyz[k];

                g_y_0_x_0_xyyy_yyzz[k] = -g_y_0_x_0_yyy_yyzz[k] * ab_x + g_y_0_x_0_yyy_xyyzz[k];

                g_y_0_x_0_xyyy_yzzz[k] = -g_y_0_x_0_yyy_yzzz[k] * ab_x + g_y_0_x_0_yyy_xyzzz[k];

                g_y_0_x_0_xyyy_zzzz[k] = -g_y_0_x_0_yyy_zzzz[k] * ab_x + g_y_0_x_0_yyy_xzzzz[k];
            }

            /// Set up 780-795 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyz_xxxx = cbuffer.data(gg_geom_1010_off + 780 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_xxxy = cbuffer.data(gg_geom_1010_off + 781 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_xxxz = cbuffer.data(gg_geom_1010_off + 782 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_xxyy = cbuffer.data(gg_geom_1010_off + 783 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_xxyz = cbuffer.data(gg_geom_1010_off + 784 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_xxzz = cbuffer.data(gg_geom_1010_off + 785 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_xyyy = cbuffer.data(gg_geom_1010_off + 786 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_xyyz = cbuffer.data(gg_geom_1010_off + 787 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_xyzz = cbuffer.data(gg_geom_1010_off + 788 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_xzzz = cbuffer.data(gg_geom_1010_off + 789 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_yyyy = cbuffer.data(gg_geom_1010_off + 790 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_yyyz = cbuffer.data(gg_geom_1010_off + 791 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_yyzz = cbuffer.data(gg_geom_1010_off + 792 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_yzzz = cbuffer.data(gg_geom_1010_off + 793 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_zzzz = cbuffer.data(gg_geom_1010_off + 794 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyz_xxxx, g_y_0_x_0_xyyz_xxxy, g_y_0_x_0_xyyz_xxxz, g_y_0_x_0_xyyz_xxyy, g_y_0_x_0_xyyz_xxyz, g_y_0_x_0_xyyz_xxzz, g_y_0_x_0_xyyz_xyyy, g_y_0_x_0_xyyz_xyyz, g_y_0_x_0_xyyz_xyzz, g_y_0_x_0_xyyz_xzzz, g_y_0_x_0_xyyz_yyyy, g_y_0_x_0_xyyz_yyyz, g_y_0_x_0_xyyz_yyzz, g_y_0_x_0_xyyz_yzzz, g_y_0_x_0_xyyz_zzzz, g_y_0_x_0_yyz_xxxx, g_y_0_x_0_yyz_xxxxx, g_y_0_x_0_yyz_xxxxy, g_y_0_x_0_yyz_xxxxz, g_y_0_x_0_yyz_xxxy, g_y_0_x_0_yyz_xxxyy, g_y_0_x_0_yyz_xxxyz, g_y_0_x_0_yyz_xxxz, g_y_0_x_0_yyz_xxxzz, g_y_0_x_0_yyz_xxyy, g_y_0_x_0_yyz_xxyyy, g_y_0_x_0_yyz_xxyyz, g_y_0_x_0_yyz_xxyz, g_y_0_x_0_yyz_xxyzz, g_y_0_x_0_yyz_xxzz, g_y_0_x_0_yyz_xxzzz, g_y_0_x_0_yyz_xyyy, g_y_0_x_0_yyz_xyyyy, g_y_0_x_0_yyz_xyyyz, g_y_0_x_0_yyz_xyyz, g_y_0_x_0_yyz_xyyzz, g_y_0_x_0_yyz_xyzz, g_y_0_x_0_yyz_xyzzz, g_y_0_x_0_yyz_xzzz, g_y_0_x_0_yyz_xzzzz, g_y_0_x_0_yyz_yyyy, g_y_0_x_0_yyz_yyyz, g_y_0_x_0_yyz_yyzz, g_y_0_x_0_yyz_yzzz, g_y_0_x_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyz_xxxx[k] = -g_y_0_x_0_yyz_xxxx[k] * ab_x + g_y_0_x_0_yyz_xxxxx[k];

                g_y_0_x_0_xyyz_xxxy[k] = -g_y_0_x_0_yyz_xxxy[k] * ab_x + g_y_0_x_0_yyz_xxxxy[k];

                g_y_0_x_0_xyyz_xxxz[k] = -g_y_0_x_0_yyz_xxxz[k] * ab_x + g_y_0_x_0_yyz_xxxxz[k];

                g_y_0_x_0_xyyz_xxyy[k] = -g_y_0_x_0_yyz_xxyy[k] * ab_x + g_y_0_x_0_yyz_xxxyy[k];

                g_y_0_x_0_xyyz_xxyz[k] = -g_y_0_x_0_yyz_xxyz[k] * ab_x + g_y_0_x_0_yyz_xxxyz[k];

                g_y_0_x_0_xyyz_xxzz[k] = -g_y_0_x_0_yyz_xxzz[k] * ab_x + g_y_0_x_0_yyz_xxxzz[k];

                g_y_0_x_0_xyyz_xyyy[k] = -g_y_0_x_0_yyz_xyyy[k] * ab_x + g_y_0_x_0_yyz_xxyyy[k];

                g_y_0_x_0_xyyz_xyyz[k] = -g_y_0_x_0_yyz_xyyz[k] * ab_x + g_y_0_x_0_yyz_xxyyz[k];

                g_y_0_x_0_xyyz_xyzz[k] = -g_y_0_x_0_yyz_xyzz[k] * ab_x + g_y_0_x_0_yyz_xxyzz[k];

                g_y_0_x_0_xyyz_xzzz[k] = -g_y_0_x_0_yyz_xzzz[k] * ab_x + g_y_0_x_0_yyz_xxzzz[k];

                g_y_0_x_0_xyyz_yyyy[k] = -g_y_0_x_0_yyz_yyyy[k] * ab_x + g_y_0_x_0_yyz_xyyyy[k];

                g_y_0_x_0_xyyz_yyyz[k] = -g_y_0_x_0_yyz_yyyz[k] * ab_x + g_y_0_x_0_yyz_xyyyz[k];

                g_y_0_x_0_xyyz_yyzz[k] = -g_y_0_x_0_yyz_yyzz[k] * ab_x + g_y_0_x_0_yyz_xyyzz[k];

                g_y_0_x_0_xyyz_yzzz[k] = -g_y_0_x_0_yyz_yzzz[k] * ab_x + g_y_0_x_0_yyz_xyzzz[k];

                g_y_0_x_0_xyyz_zzzz[k] = -g_y_0_x_0_yyz_zzzz[k] * ab_x + g_y_0_x_0_yyz_xzzzz[k];
            }

            /// Set up 795-810 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyzz_xxxx = cbuffer.data(gg_geom_1010_off + 795 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_xxxy = cbuffer.data(gg_geom_1010_off + 796 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_xxxz = cbuffer.data(gg_geom_1010_off + 797 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_xxyy = cbuffer.data(gg_geom_1010_off + 798 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_xxyz = cbuffer.data(gg_geom_1010_off + 799 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_xxzz = cbuffer.data(gg_geom_1010_off + 800 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_xyyy = cbuffer.data(gg_geom_1010_off + 801 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_xyyz = cbuffer.data(gg_geom_1010_off + 802 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_xyzz = cbuffer.data(gg_geom_1010_off + 803 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_xzzz = cbuffer.data(gg_geom_1010_off + 804 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_yyyy = cbuffer.data(gg_geom_1010_off + 805 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_yyyz = cbuffer.data(gg_geom_1010_off + 806 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_yyzz = cbuffer.data(gg_geom_1010_off + 807 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_yzzz = cbuffer.data(gg_geom_1010_off + 808 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_zzzz = cbuffer.data(gg_geom_1010_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyzz_xxxx, g_y_0_x_0_xyzz_xxxy, g_y_0_x_0_xyzz_xxxz, g_y_0_x_0_xyzz_xxyy, g_y_0_x_0_xyzz_xxyz, g_y_0_x_0_xyzz_xxzz, g_y_0_x_0_xyzz_xyyy, g_y_0_x_0_xyzz_xyyz, g_y_0_x_0_xyzz_xyzz, g_y_0_x_0_xyzz_xzzz, g_y_0_x_0_xyzz_yyyy, g_y_0_x_0_xyzz_yyyz, g_y_0_x_0_xyzz_yyzz, g_y_0_x_0_xyzz_yzzz, g_y_0_x_0_xyzz_zzzz, g_y_0_x_0_yzz_xxxx, g_y_0_x_0_yzz_xxxxx, g_y_0_x_0_yzz_xxxxy, g_y_0_x_0_yzz_xxxxz, g_y_0_x_0_yzz_xxxy, g_y_0_x_0_yzz_xxxyy, g_y_0_x_0_yzz_xxxyz, g_y_0_x_0_yzz_xxxz, g_y_0_x_0_yzz_xxxzz, g_y_0_x_0_yzz_xxyy, g_y_0_x_0_yzz_xxyyy, g_y_0_x_0_yzz_xxyyz, g_y_0_x_0_yzz_xxyz, g_y_0_x_0_yzz_xxyzz, g_y_0_x_0_yzz_xxzz, g_y_0_x_0_yzz_xxzzz, g_y_0_x_0_yzz_xyyy, g_y_0_x_0_yzz_xyyyy, g_y_0_x_0_yzz_xyyyz, g_y_0_x_0_yzz_xyyz, g_y_0_x_0_yzz_xyyzz, g_y_0_x_0_yzz_xyzz, g_y_0_x_0_yzz_xyzzz, g_y_0_x_0_yzz_xzzz, g_y_0_x_0_yzz_xzzzz, g_y_0_x_0_yzz_yyyy, g_y_0_x_0_yzz_yyyz, g_y_0_x_0_yzz_yyzz, g_y_0_x_0_yzz_yzzz, g_y_0_x_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyzz_xxxx[k] = -g_y_0_x_0_yzz_xxxx[k] * ab_x + g_y_0_x_0_yzz_xxxxx[k];

                g_y_0_x_0_xyzz_xxxy[k] = -g_y_0_x_0_yzz_xxxy[k] * ab_x + g_y_0_x_0_yzz_xxxxy[k];

                g_y_0_x_0_xyzz_xxxz[k] = -g_y_0_x_0_yzz_xxxz[k] * ab_x + g_y_0_x_0_yzz_xxxxz[k];

                g_y_0_x_0_xyzz_xxyy[k] = -g_y_0_x_0_yzz_xxyy[k] * ab_x + g_y_0_x_0_yzz_xxxyy[k];

                g_y_0_x_0_xyzz_xxyz[k] = -g_y_0_x_0_yzz_xxyz[k] * ab_x + g_y_0_x_0_yzz_xxxyz[k];

                g_y_0_x_0_xyzz_xxzz[k] = -g_y_0_x_0_yzz_xxzz[k] * ab_x + g_y_0_x_0_yzz_xxxzz[k];

                g_y_0_x_0_xyzz_xyyy[k] = -g_y_0_x_0_yzz_xyyy[k] * ab_x + g_y_0_x_0_yzz_xxyyy[k];

                g_y_0_x_0_xyzz_xyyz[k] = -g_y_0_x_0_yzz_xyyz[k] * ab_x + g_y_0_x_0_yzz_xxyyz[k];

                g_y_0_x_0_xyzz_xyzz[k] = -g_y_0_x_0_yzz_xyzz[k] * ab_x + g_y_0_x_0_yzz_xxyzz[k];

                g_y_0_x_0_xyzz_xzzz[k] = -g_y_0_x_0_yzz_xzzz[k] * ab_x + g_y_0_x_0_yzz_xxzzz[k];

                g_y_0_x_0_xyzz_yyyy[k] = -g_y_0_x_0_yzz_yyyy[k] * ab_x + g_y_0_x_0_yzz_xyyyy[k];

                g_y_0_x_0_xyzz_yyyz[k] = -g_y_0_x_0_yzz_yyyz[k] * ab_x + g_y_0_x_0_yzz_xyyyz[k];

                g_y_0_x_0_xyzz_yyzz[k] = -g_y_0_x_0_yzz_yyzz[k] * ab_x + g_y_0_x_0_yzz_xyyzz[k];

                g_y_0_x_0_xyzz_yzzz[k] = -g_y_0_x_0_yzz_yzzz[k] * ab_x + g_y_0_x_0_yzz_xyzzz[k];

                g_y_0_x_0_xyzz_zzzz[k] = -g_y_0_x_0_yzz_zzzz[k] * ab_x + g_y_0_x_0_yzz_xzzzz[k];
            }

            /// Set up 810-825 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xzzz_xxxx = cbuffer.data(gg_geom_1010_off + 810 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_xxxy = cbuffer.data(gg_geom_1010_off + 811 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_xxxz = cbuffer.data(gg_geom_1010_off + 812 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_xxyy = cbuffer.data(gg_geom_1010_off + 813 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_xxyz = cbuffer.data(gg_geom_1010_off + 814 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_xxzz = cbuffer.data(gg_geom_1010_off + 815 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_xyyy = cbuffer.data(gg_geom_1010_off + 816 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_xyyz = cbuffer.data(gg_geom_1010_off + 817 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_xyzz = cbuffer.data(gg_geom_1010_off + 818 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_xzzz = cbuffer.data(gg_geom_1010_off + 819 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_yyyy = cbuffer.data(gg_geom_1010_off + 820 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_yyyz = cbuffer.data(gg_geom_1010_off + 821 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_yyzz = cbuffer.data(gg_geom_1010_off + 822 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_yzzz = cbuffer.data(gg_geom_1010_off + 823 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_zzzz = cbuffer.data(gg_geom_1010_off + 824 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xzzz_xxxx, g_y_0_x_0_xzzz_xxxy, g_y_0_x_0_xzzz_xxxz, g_y_0_x_0_xzzz_xxyy, g_y_0_x_0_xzzz_xxyz, g_y_0_x_0_xzzz_xxzz, g_y_0_x_0_xzzz_xyyy, g_y_0_x_0_xzzz_xyyz, g_y_0_x_0_xzzz_xyzz, g_y_0_x_0_xzzz_xzzz, g_y_0_x_0_xzzz_yyyy, g_y_0_x_0_xzzz_yyyz, g_y_0_x_0_xzzz_yyzz, g_y_0_x_0_xzzz_yzzz, g_y_0_x_0_xzzz_zzzz, g_y_0_x_0_zzz_xxxx, g_y_0_x_0_zzz_xxxxx, g_y_0_x_0_zzz_xxxxy, g_y_0_x_0_zzz_xxxxz, g_y_0_x_0_zzz_xxxy, g_y_0_x_0_zzz_xxxyy, g_y_0_x_0_zzz_xxxyz, g_y_0_x_0_zzz_xxxz, g_y_0_x_0_zzz_xxxzz, g_y_0_x_0_zzz_xxyy, g_y_0_x_0_zzz_xxyyy, g_y_0_x_0_zzz_xxyyz, g_y_0_x_0_zzz_xxyz, g_y_0_x_0_zzz_xxyzz, g_y_0_x_0_zzz_xxzz, g_y_0_x_0_zzz_xxzzz, g_y_0_x_0_zzz_xyyy, g_y_0_x_0_zzz_xyyyy, g_y_0_x_0_zzz_xyyyz, g_y_0_x_0_zzz_xyyz, g_y_0_x_0_zzz_xyyzz, g_y_0_x_0_zzz_xyzz, g_y_0_x_0_zzz_xyzzz, g_y_0_x_0_zzz_xzzz, g_y_0_x_0_zzz_xzzzz, g_y_0_x_0_zzz_yyyy, g_y_0_x_0_zzz_yyyz, g_y_0_x_0_zzz_yyzz, g_y_0_x_0_zzz_yzzz, g_y_0_x_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xzzz_xxxx[k] = -g_y_0_x_0_zzz_xxxx[k] * ab_x + g_y_0_x_0_zzz_xxxxx[k];

                g_y_0_x_0_xzzz_xxxy[k] = -g_y_0_x_0_zzz_xxxy[k] * ab_x + g_y_0_x_0_zzz_xxxxy[k];

                g_y_0_x_0_xzzz_xxxz[k] = -g_y_0_x_0_zzz_xxxz[k] * ab_x + g_y_0_x_0_zzz_xxxxz[k];

                g_y_0_x_0_xzzz_xxyy[k] = -g_y_0_x_0_zzz_xxyy[k] * ab_x + g_y_0_x_0_zzz_xxxyy[k];

                g_y_0_x_0_xzzz_xxyz[k] = -g_y_0_x_0_zzz_xxyz[k] * ab_x + g_y_0_x_0_zzz_xxxyz[k];

                g_y_0_x_0_xzzz_xxzz[k] = -g_y_0_x_0_zzz_xxzz[k] * ab_x + g_y_0_x_0_zzz_xxxzz[k];

                g_y_0_x_0_xzzz_xyyy[k] = -g_y_0_x_0_zzz_xyyy[k] * ab_x + g_y_0_x_0_zzz_xxyyy[k];

                g_y_0_x_0_xzzz_xyyz[k] = -g_y_0_x_0_zzz_xyyz[k] * ab_x + g_y_0_x_0_zzz_xxyyz[k];

                g_y_0_x_0_xzzz_xyzz[k] = -g_y_0_x_0_zzz_xyzz[k] * ab_x + g_y_0_x_0_zzz_xxyzz[k];

                g_y_0_x_0_xzzz_xzzz[k] = -g_y_0_x_0_zzz_xzzz[k] * ab_x + g_y_0_x_0_zzz_xxzzz[k];

                g_y_0_x_0_xzzz_yyyy[k] = -g_y_0_x_0_zzz_yyyy[k] * ab_x + g_y_0_x_0_zzz_xyyyy[k];

                g_y_0_x_0_xzzz_yyyz[k] = -g_y_0_x_0_zzz_yyyz[k] * ab_x + g_y_0_x_0_zzz_xyyyz[k];

                g_y_0_x_0_xzzz_yyzz[k] = -g_y_0_x_0_zzz_yyzz[k] * ab_x + g_y_0_x_0_zzz_xyyzz[k];

                g_y_0_x_0_xzzz_yzzz[k] = -g_y_0_x_0_zzz_yzzz[k] * ab_x + g_y_0_x_0_zzz_xyzzz[k];

                g_y_0_x_0_xzzz_zzzz[k] = -g_y_0_x_0_zzz_zzzz[k] * ab_x + g_y_0_x_0_zzz_xzzzz[k];
            }

            /// Set up 825-840 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyy_xxxx = cbuffer.data(gg_geom_1010_off + 825 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_xxxy = cbuffer.data(gg_geom_1010_off + 826 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_xxxz = cbuffer.data(gg_geom_1010_off + 827 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_xxyy = cbuffer.data(gg_geom_1010_off + 828 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_xxyz = cbuffer.data(gg_geom_1010_off + 829 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_xxzz = cbuffer.data(gg_geom_1010_off + 830 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_xyyy = cbuffer.data(gg_geom_1010_off + 831 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_xyyz = cbuffer.data(gg_geom_1010_off + 832 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_xyzz = cbuffer.data(gg_geom_1010_off + 833 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_xzzz = cbuffer.data(gg_geom_1010_off + 834 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_yyyy = cbuffer.data(gg_geom_1010_off + 835 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_yyyz = cbuffer.data(gg_geom_1010_off + 836 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_yyzz = cbuffer.data(gg_geom_1010_off + 837 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_yzzz = cbuffer.data(gg_geom_1010_off + 838 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_zzzz = cbuffer.data(gg_geom_1010_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_yyy_xxxx, g_0_0_x_0_yyy_xxxy, g_0_0_x_0_yyy_xxxz, g_0_0_x_0_yyy_xxyy, g_0_0_x_0_yyy_xxyz, g_0_0_x_0_yyy_xxzz, g_0_0_x_0_yyy_xyyy, g_0_0_x_0_yyy_xyyz, g_0_0_x_0_yyy_xyzz, g_0_0_x_0_yyy_xzzz, g_0_0_x_0_yyy_yyyy, g_0_0_x_0_yyy_yyyz, g_0_0_x_0_yyy_yyzz, g_0_0_x_0_yyy_yzzz, g_0_0_x_0_yyy_zzzz, g_y_0_x_0_yyy_xxxx, g_y_0_x_0_yyy_xxxxy, g_y_0_x_0_yyy_xxxy, g_y_0_x_0_yyy_xxxyy, g_y_0_x_0_yyy_xxxyz, g_y_0_x_0_yyy_xxxz, g_y_0_x_0_yyy_xxyy, g_y_0_x_0_yyy_xxyyy, g_y_0_x_0_yyy_xxyyz, g_y_0_x_0_yyy_xxyz, g_y_0_x_0_yyy_xxyzz, g_y_0_x_0_yyy_xxzz, g_y_0_x_0_yyy_xyyy, g_y_0_x_0_yyy_xyyyy, g_y_0_x_0_yyy_xyyyz, g_y_0_x_0_yyy_xyyz, g_y_0_x_0_yyy_xyyzz, g_y_0_x_0_yyy_xyzz, g_y_0_x_0_yyy_xyzzz, g_y_0_x_0_yyy_xzzz, g_y_0_x_0_yyy_yyyy, g_y_0_x_0_yyy_yyyyy, g_y_0_x_0_yyy_yyyyz, g_y_0_x_0_yyy_yyyz, g_y_0_x_0_yyy_yyyzz, g_y_0_x_0_yyy_yyzz, g_y_0_x_0_yyy_yyzzz, g_y_0_x_0_yyy_yzzz, g_y_0_x_0_yyy_yzzzz, g_y_0_x_0_yyy_zzzz, g_y_0_x_0_yyyy_xxxx, g_y_0_x_0_yyyy_xxxy, g_y_0_x_0_yyyy_xxxz, g_y_0_x_0_yyyy_xxyy, g_y_0_x_0_yyyy_xxyz, g_y_0_x_0_yyyy_xxzz, g_y_0_x_0_yyyy_xyyy, g_y_0_x_0_yyyy_xyyz, g_y_0_x_0_yyyy_xyzz, g_y_0_x_0_yyyy_xzzz, g_y_0_x_0_yyyy_yyyy, g_y_0_x_0_yyyy_yyyz, g_y_0_x_0_yyyy_yyzz, g_y_0_x_0_yyyy_yzzz, g_y_0_x_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyy_xxxx[k] = -g_0_0_x_0_yyy_xxxx[k] - g_y_0_x_0_yyy_xxxx[k] * ab_y + g_y_0_x_0_yyy_xxxxy[k];

                g_y_0_x_0_yyyy_xxxy[k] = -g_0_0_x_0_yyy_xxxy[k] - g_y_0_x_0_yyy_xxxy[k] * ab_y + g_y_0_x_0_yyy_xxxyy[k];

                g_y_0_x_0_yyyy_xxxz[k] = -g_0_0_x_0_yyy_xxxz[k] - g_y_0_x_0_yyy_xxxz[k] * ab_y + g_y_0_x_0_yyy_xxxyz[k];

                g_y_0_x_0_yyyy_xxyy[k] = -g_0_0_x_0_yyy_xxyy[k] - g_y_0_x_0_yyy_xxyy[k] * ab_y + g_y_0_x_0_yyy_xxyyy[k];

                g_y_0_x_0_yyyy_xxyz[k] = -g_0_0_x_0_yyy_xxyz[k] - g_y_0_x_0_yyy_xxyz[k] * ab_y + g_y_0_x_0_yyy_xxyyz[k];

                g_y_0_x_0_yyyy_xxzz[k] = -g_0_0_x_0_yyy_xxzz[k] - g_y_0_x_0_yyy_xxzz[k] * ab_y + g_y_0_x_0_yyy_xxyzz[k];

                g_y_0_x_0_yyyy_xyyy[k] = -g_0_0_x_0_yyy_xyyy[k] - g_y_0_x_0_yyy_xyyy[k] * ab_y + g_y_0_x_0_yyy_xyyyy[k];

                g_y_0_x_0_yyyy_xyyz[k] = -g_0_0_x_0_yyy_xyyz[k] - g_y_0_x_0_yyy_xyyz[k] * ab_y + g_y_0_x_0_yyy_xyyyz[k];

                g_y_0_x_0_yyyy_xyzz[k] = -g_0_0_x_0_yyy_xyzz[k] - g_y_0_x_0_yyy_xyzz[k] * ab_y + g_y_0_x_0_yyy_xyyzz[k];

                g_y_0_x_0_yyyy_xzzz[k] = -g_0_0_x_0_yyy_xzzz[k] - g_y_0_x_0_yyy_xzzz[k] * ab_y + g_y_0_x_0_yyy_xyzzz[k];

                g_y_0_x_0_yyyy_yyyy[k] = -g_0_0_x_0_yyy_yyyy[k] - g_y_0_x_0_yyy_yyyy[k] * ab_y + g_y_0_x_0_yyy_yyyyy[k];

                g_y_0_x_0_yyyy_yyyz[k] = -g_0_0_x_0_yyy_yyyz[k] - g_y_0_x_0_yyy_yyyz[k] * ab_y + g_y_0_x_0_yyy_yyyyz[k];

                g_y_0_x_0_yyyy_yyzz[k] = -g_0_0_x_0_yyy_yyzz[k] - g_y_0_x_0_yyy_yyzz[k] * ab_y + g_y_0_x_0_yyy_yyyzz[k];

                g_y_0_x_0_yyyy_yzzz[k] = -g_0_0_x_0_yyy_yzzz[k] - g_y_0_x_0_yyy_yzzz[k] * ab_y + g_y_0_x_0_yyy_yyzzz[k];

                g_y_0_x_0_yyyy_zzzz[k] = -g_0_0_x_0_yyy_zzzz[k] - g_y_0_x_0_yyy_zzzz[k] * ab_y + g_y_0_x_0_yyy_yzzzz[k];
            }

            /// Set up 840-855 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyz_xxxx = cbuffer.data(gg_geom_1010_off + 840 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_xxxy = cbuffer.data(gg_geom_1010_off + 841 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_xxxz = cbuffer.data(gg_geom_1010_off + 842 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_xxyy = cbuffer.data(gg_geom_1010_off + 843 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_xxyz = cbuffer.data(gg_geom_1010_off + 844 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_xxzz = cbuffer.data(gg_geom_1010_off + 845 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_xyyy = cbuffer.data(gg_geom_1010_off + 846 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_xyyz = cbuffer.data(gg_geom_1010_off + 847 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_xyzz = cbuffer.data(gg_geom_1010_off + 848 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_xzzz = cbuffer.data(gg_geom_1010_off + 849 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_yyyy = cbuffer.data(gg_geom_1010_off + 850 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_yyyz = cbuffer.data(gg_geom_1010_off + 851 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_yyzz = cbuffer.data(gg_geom_1010_off + 852 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_yzzz = cbuffer.data(gg_geom_1010_off + 853 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_zzzz = cbuffer.data(gg_geom_1010_off + 854 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyy_xxxx, g_y_0_x_0_yyy_xxxxz, g_y_0_x_0_yyy_xxxy, g_y_0_x_0_yyy_xxxyz, g_y_0_x_0_yyy_xxxz, g_y_0_x_0_yyy_xxxzz, g_y_0_x_0_yyy_xxyy, g_y_0_x_0_yyy_xxyyz, g_y_0_x_0_yyy_xxyz, g_y_0_x_0_yyy_xxyzz, g_y_0_x_0_yyy_xxzz, g_y_0_x_0_yyy_xxzzz, g_y_0_x_0_yyy_xyyy, g_y_0_x_0_yyy_xyyyz, g_y_0_x_0_yyy_xyyz, g_y_0_x_0_yyy_xyyzz, g_y_0_x_0_yyy_xyzz, g_y_0_x_0_yyy_xyzzz, g_y_0_x_0_yyy_xzzz, g_y_0_x_0_yyy_xzzzz, g_y_0_x_0_yyy_yyyy, g_y_0_x_0_yyy_yyyyz, g_y_0_x_0_yyy_yyyz, g_y_0_x_0_yyy_yyyzz, g_y_0_x_0_yyy_yyzz, g_y_0_x_0_yyy_yyzzz, g_y_0_x_0_yyy_yzzz, g_y_0_x_0_yyy_yzzzz, g_y_0_x_0_yyy_zzzz, g_y_0_x_0_yyy_zzzzz, g_y_0_x_0_yyyz_xxxx, g_y_0_x_0_yyyz_xxxy, g_y_0_x_0_yyyz_xxxz, g_y_0_x_0_yyyz_xxyy, g_y_0_x_0_yyyz_xxyz, g_y_0_x_0_yyyz_xxzz, g_y_0_x_0_yyyz_xyyy, g_y_0_x_0_yyyz_xyyz, g_y_0_x_0_yyyz_xyzz, g_y_0_x_0_yyyz_xzzz, g_y_0_x_0_yyyz_yyyy, g_y_0_x_0_yyyz_yyyz, g_y_0_x_0_yyyz_yyzz, g_y_0_x_0_yyyz_yzzz, g_y_0_x_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyz_xxxx[k] = -g_y_0_x_0_yyy_xxxx[k] * ab_z + g_y_0_x_0_yyy_xxxxz[k];

                g_y_0_x_0_yyyz_xxxy[k] = -g_y_0_x_0_yyy_xxxy[k] * ab_z + g_y_0_x_0_yyy_xxxyz[k];

                g_y_0_x_0_yyyz_xxxz[k] = -g_y_0_x_0_yyy_xxxz[k] * ab_z + g_y_0_x_0_yyy_xxxzz[k];

                g_y_0_x_0_yyyz_xxyy[k] = -g_y_0_x_0_yyy_xxyy[k] * ab_z + g_y_0_x_0_yyy_xxyyz[k];

                g_y_0_x_0_yyyz_xxyz[k] = -g_y_0_x_0_yyy_xxyz[k] * ab_z + g_y_0_x_0_yyy_xxyzz[k];

                g_y_0_x_0_yyyz_xxzz[k] = -g_y_0_x_0_yyy_xxzz[k] * ab_z + g_y_0_x_0_yyy_xxzzz[k];

                g_y_0_x_0_yyyz_xyyy[k] = -g_y_0_x_0_yyy_xyyy[k] * ab_z + g_y_0_x_0_yyy_xyyyz[k];

                g_y_0_x_0_yyyz_xyyz[k] = -g_y_0_x_0_yyy_xyyz[k] * ab_z + g_y_0_x_0_yyy_xyyzz[k];

                g_y_0_x_0_yyyz_xyzz[k] = -g_y_0_x_0_yyy_xyzz[k] * ab_z + g_y_0_x_0_yyy_xyzzz[k];

                g_y_0_x_0_yyyz_xzzz[k] = -g_y_0_x_0_yyy_xzzz[k] * ab_z + g_y_0_x_0_yyy_xzzzz[k];

                g_y_0_x_0_yyyz_yyyy[k] = -g_y_0_x_0_yyy_yyyy[k] * ab_z + g_y_0_x_0_yyy_yyyyz[k];

                g_y_0_x_0_yyyz_yyyz[k] = -g_y_0_x_0_yyy_yyyz[k] * ab_z + g_y_0_x_0_yyy_yyyzz[k];

                g_y_0_x_0_yyyz_yyzz[k] = -g_y_0_x_0_yyy_yyzz[k] * ab_z + g_y_0_x_0_yyy_yyzzz[k];

                g_y_0_x_0_yyyz_yzzz[k] = -g_y_0_x_0_yyy_yzzz[k] * ab_z + g_y_0_x_0_yyy_yzzzz[k];

                g_y_0_x_0_yyyz_zzzz[k] = -g_y_0_x_0_yyy_zzzz[k] * ab_z + g_y_0_x_0_yyy_zzzzz[k];
            }

            /// Set up 855-870 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyzz_xxxx = cbuffer.data(gg_geom_1010_off + 855 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_xxxy = cbuffer.data(gg_geom_1010_off + 856 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_xxxz = cbuffer.data(gg_geom_1010_off + 857 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_xxyy = cbuffer.data(gg_geom_1010_off + 858 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_xxyz = cbuffer.data(gg_geom_1010_off + 859 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_xxzz = cbuffer.data(gg_geom_1010_off + 860 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_xyyy = cbuffer.data(gg_geom_1010_off + 861 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_xyyz = cbuffer.data(gg_geom_1010_off + 862 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_xyzz = cbuffer.data(gg_geom_1010_off + 863 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_xzzz = cbuffer.data(gg_geom_1010_off + 864 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_yyyy = cbuffer.data(gg_geom_1010_off + 865 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_yyyz = cbuffer.data(gg_geom_1010_off + 866 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_yyzz = cbuffer.data(gg_geom_1010_off + 867 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_yzzz = cbuffer.data(gg_geom_1010_off + 868 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_zzzz = cbuffer.data(gg_geom_1010_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyz_xxxx, g_y_0_x_0_yyz_xxxxz, g_y_0_x_0_yyz_xxxy, g_y_0_x_0_yyz_xxxyz, g_y_0_x_0_yyz_xxxz, g_y_0_x_0_yyz_xxxzz, g_y_0_x_0_yyz_xxyy, g_y_0_x_0_yyz_xxyyz, g_y_0_x_0_yyz_xxyz, g_y_0_x_0_yyz_xxyzz, g_y_0_x_0_yyz_xxzz, g_y_0_x_0_yyz_xxzzz, g_y_0_x_0_yyz_xyyy, g_y_0_x_0_yyz_xyyyz, g_y_0_x_0_yyz_xyyz, g_y_0_x_0_yyz_xyyzz, g_y_0_x_0_yyz_xyzz, g_y_0_x_0_yyz_xyzzz, g_y_0_x_0_yyz_xzzz, g_y_0_x_0_yyz_xzzzz, g_y_0_x_0_yyz_yyyy, g_y_0_x_0_yyz_yyyyz, g_y_0_x_0_yyz_yyyz, g_y_0_x_0_yyz_yyyzz, g_y_0_x_0_yyz_yyzz, g_y_0_x_0_yyz_yyzzz, g_y_0_x_0_yyz_yzzz, g_y_0_x_0_yyz_yzzzz, g_y_0_x_0_yyz_zzzz, g_y_0_x_0_yyz_zzzzz, g_y_0_x_0_yyzz_xxxx, g_y_0_x_0_yyzz_xxxy, g_y_0_x_0_yyzz_xxxz, g_y_0_x_0_yyzz_xxyy, g_y_0_x_0_yyzz_xxyz, g_y_0_x_0_yyzz_xxzz, g_y_0_x_0_yyzz_xyyy, g_y_0_x_0_yyzz_xyyz, g_y_0_x_0_yyzz_xyzz, g_y_0_x_0_yyzz_xzzz, g_y_0_x_0_yyzz_yyyy, g_y_0_x_0_yyzz_yyyz, g_y_0_x_0_yyzz_yyzz, g_y_0_x_0_yyzz_yzzz, g_y_0_x_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyzz_xxxx[k] = -g_y_0_x_0_yyz_xxxx[k] * ab_z + g_y_0_x_0_yyz_xxxxz[k];

                g_y_0_x_0_yyzz_xxxy[k] = -g_y_0_x_0_yyz_xxxy[k] * ab_z + g_y_0_x_0_yyz_xxxyz[k];

                g_y_0_x_0_yyzz_xxxz[k] = -g_y_0_x_0_yyz_xxxz[k] * ab_z + g_y_0_x_0_yyz_xxxzz[k];

                g_y_0_x_0_yyzz_xxyy[k] = -g_y_0_x_0_yyz_xxyy[k] * ab_z + g_y_0_x_0_yyz_xxyyz[k];

                g_y_0_x_0_yyzz_xxyz[k] = -g_y_0_x_0_yyz_xxyz[k] * ab_z + g_y_0_x_0_yyz_xxyzz[k];

                g_y_0_x_0_yyzz_xxzz[k] = -g_y_0_x_0_yyz_xxzz[k] * ab_z + g_y_0_x_0_yyz_xxzzz[k];

                g_y_0_x_0_yyzz_xyyy[k] = -g_y_0_x_0_yyz_xyyy[k] * ab_z + g_y_0_x_0_yyz_xyyyz[k];

                g_y_0_x_0_yyzz_xyyz[k] = -g_y_0_x_0_yyz_xyyz[k] * ab_z + g_y_0_x_0_yyz_xyyzz[k];

                g_y_0_x_0_yyzz_xyzz[k] = -g_y_0_x_0_yyz_xyzz[k] * ab_z + g_y_0_x_0_yyz_xyzzz[k];

                g_y_0_x_0_yyzz_xzzz[k] = -g_y_0_x_0_yyz_xzzz[k] * ab_z + g_y_0_x_0_yyz_xzzzz[k];

                g_y_0_x_0_yyzz_yyyy[k] = -g_y_0_x_0_yyz_yyyy[k] * ab_z + g_y_0_x_0_yyz_yyyyz[k];

                g_y_0_x_0_yyzz_yyyz[k] = -g_y_0_x_0_yyz_yyyz[k] * ab_z + g_y_0_x_0_yyz_yyyzz[k];

                g_y_0_x_0_yyzz_yyzz[k] = -g_y_0_x_0_yyz_yyzz[k] * ab_z + g_y_0_x_0_yyz_yyzzz[k];

                g_y_0_x_0_yyzz_yzzz[k] = -g_y_0_x_0_yyz_yzzz[k] * ab_z + g_y_0_x_0_yyz_yzzzz[k];

                g_y_0_x_0_yyzz_zzzz[k] = -g_y_0_x_0_yyz_zzzz[k] * ab_z + g_y_0_x_0_yyz_zzzzz[k];
            }

            /// Set up 870-885 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yzzz_xxxx = cbuffer.data(gg_geom_1010_off + 870 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_xxxy = cbuffer.data(gg_geom_1010_off + 871 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_xxxz = cbuffer.data(gg_geom_1010_off + 872 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_xxyy = cbuffer.data(gg_geom_1010_off + 873 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_xxyz = cbuffer.data(gg_geom_1010_off + 874 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_xxzz = cbuffer.data(gg_geom_1010_off + 875 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_xyyy = cbuffer.data(gg_geom_1010_off + 876 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_xyyz = cbuffer.data(gg_geom_1010_off + 877 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_xyzz = cbuffer.data(gg_geom_1010_off + 878 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_xzzz = cbuffer.data(gg_geom_1010_off + 879 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_yyyy = cbuffer.data(gg_geom_1010_off + 880 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_yyyz = cbuffer.data(gg_geom_1010_off + 881 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_yyzz = cbuffer.data(gg_geom_1010_off + 882 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_yzzz = cbuffer.data(gg_geom_1010_off + 883 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_zzzz = cbuffer.data(gg_geom_1010_off + 884 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yzz_xxxx, g_y_0_x_0_yzz_xxxxz, g_y_0_x_0_yzz_xxxy, g_y_0_x_0_yzz_xxxyz, g_y_0_x_0_yzz_xxxz, g_y_0_x_0_yzz_xxxzz, g_y_0_x_0_yzz_xxyy, g_y_0_x_0_yzz_xxyyz, g_y_0_x_0_yzz_xxyz, g_y_0_x_0_yzz_xxyzz, g_y_0_x_0_yzz_xxzz, g_y_0_x_0_yzz_xxzzz, g_y_0_x_0_yzz_xyyy, g_y_0_x_0_yzz_xyyyz, g_y_0_x_0_yzz_xyyz, g_y_0_x_0_yzz_xyyzz, g_y_0_x_0_yzz_xyzz, g_y_0_x_0_yzz_xyzzz, g_y_0_x_0_yzz_xzzz, g_y_0_x_0_yzz_xzzzz, g_y_0_x_0_yzz_yyyy, g_y_0_x_0_yzz_yyyyz, g_y_0_x_0_yzz_yyyz, g_y_0_x_0_yzz_yyyzz, g_y_0_x_0_yzz_yyzz, g_y_0_x_0_yzz_yyzzz, g_y_0_x_0_yzz_yzzz, g_y_0_x_0_yzz_yzzzz, g_y_0_x_0_yzz_zzzz, g_y_0_x_0_yzz_zzzzz, g_y_0_x_0_yzzz_xxxx, g_y_0_x_0_yzzz_xxxy, g_y_0_x_0_yzzz_xxxz, g_y_0_x_0_yzzz_xxyy, g_y_0_x_0_yzzz_xxyz, g_y_0_x_0_yzzz_xxzz, g_y_0_x_0_yzzz_xyyy, g_y_0_x_0_yzzz_xyyz, g_y_0_x_0_yzzz_xyzz, g_y_0_x_0_yzzz_xzzz, g_y_0_x_0_yzzz_yyyy, g_y_0_x_0_yzzz_yyyz, g_y_0_x_0_yzzz_yyzz, g_y_0_x_0_yzzz_yzzz, g_y_0_x_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yzzz_xxxx[k] = -g_y_0_x_0_yzz_xxxx[k] * ab_z + g_y_0_x_0_yzz_xxxxz[k];

                g_y_0_x_0_yzzz_xxxy[k] = -g_y_0_x_0_yzz_xxxy[k] * ab_z + g_y_0_x_0_yzz_xxxyz[k];

                g_y_0_x_0_yzzz_xxxz[k] = -g_y_0_x_0_yzz_xxxz[k] * ab_z + g_y_0_x_0_yzz_xxxzz[k];

                g_y_0_x_0_yzzz_xxyy[k] = -g_y_0_x_0_yzz_xxyy[k] * ab_z + g_y_0_x_0_yzz_xxyyz[k];

                g_y_0_x_0_yzzz_xxyz[k] = -g_y_0_x_0_yzz_xxyz[k] * ab_z + g_y_0_x_0_yzz_xxyzz[k];

                g_y_0_x_0_yzzz_xxzz[k] = -g_y_0_x_0_yzz_xxzz[k] * ab_z + g_y_0_x_0_yzz_xxzzz[k];

                g_y_0_x_0_yzzz_xyyy[k] = -g_y_0_x_0_yzz_xyyy[k] * ab_z + g_y_0_x_0_yzz_xyyyz[k];

                g_y_0_x_0_yzzz_xyyz[k] = -g_y_0_x_0_yzz_xyyz[k] * ab_z + g_y_0_x_0_yzz_xyyzz[k];

                g_y_0_x_0_yzzz_xyzz[k] = -g_y_0_x_0_yzz_xyzz[k] * ab_z + g_y_0_x_0_yzz_xyzzz[k];

                g_y_0_x_0_yzzz_xzzz[k] = -g_y_0_x_0_yzz_xzzz[k] * ab_z + g_y_0_x_0_yzz_xzzzz[k];

                g_y_0_x_0_yzzz_yyyy[k] = -g_y_0_x_0_yzz_yyyy[k] * ab_z + g_y_0_x_0_yzz_yyyyz[k];

                g_y_0_x_0_yzzz_yyyz[k] = -g_y_0_x_0_yzz_yyyz[k] * ab_z + g_y_0_x_0_yzz_yyyzz[k];

                g_y_0_x_0_yzzz_yyzz[k] = -g_y_0_x_0_yzz_yyzz[k] * ab_z + g_y_0_x_0_yzz_yyzzz[k];

                g_y_0_x_0_yzzz_yzzz[k] = -g_y_0_x_0_yzz_yzzz[k] * ab_z + g_y_0_x_0_yzz_yzzzz[k];

                g_y_0_x_0_yzzz_zzzz[k] = -g_y_0_x_0_yzz_zzzz[k] * ab_z + g_y_0_x_0_yzz_zzzzz[k];
            }

            /// Set up 885-900 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_zzzz_xxxx = cbuffer.data(gg_geom_1010_off + 885 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_xxxy = cbuffer.data(gg_geom_1010_off + 886 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_xxxz = cbuffer.data(gg_geom_1010_off + 887 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_xxyy = cbuffer.data(gg_geom_1010_off + 888 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_xxyz = cbuffer.data(gg_geom_1010_off + 889 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_xxzz = cbuffer.data(gg_geom_1010_off + 890 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_xyyy = cbuffer.data(gg_geom_1010_off + 891 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_xyyz = cbuffer.data(gg_geom_1010_off + 892 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_xyzz = cbuffer.data(gg_geom_1010_off + 893 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_xzzz = cbuffer.data(gg_geom_1010_off + 894 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_yyyy = cbuffer.data(gg_geom_1010_off + 895 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_yyyz = cbuffer.data(gg_geom_1010_off + 896 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_yyzz = cbuffer.data(gg_geom_1010_off + 897 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_yzzz = cbuffer.data(gg_geom_1010_off + 898 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_zzzz = cbuffer.data(gg_geom_1010_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_zzz_xxxx, g_y_0_x_0_zzz_xxxxz, g_y_0_x_0_zzz_xxxy, g_y_0_x_0_zzz_xxxyz, g_y_0_x_0_zzz_xxxz, g_y_0_x_0_zzz_xxxzz, g_y_0_x_0_zzz_xxyy, g_y_0_x_0_zzz_xxyyz, g_y_0_x_0_zzz_xxyz, g_y_0_x_0_zzz_xxyzz, g_y_0_x_0_zzz_xxzz, g_y_0_x_0_zzz_xxzzz, g_y_0_x_0_zzz_xyyy, g_y_0_x_0_zzz_xyyyz, g_y_0_x_0_zzz_xyyz, g_y_0_x_0_zzz_xyyzz, g_y_0_x_0_zzz_xyzz, g_y_0_x_0_zzz_xyzzz, g_y_0_x_0_zzz_xzzz, g_y_0_x_0_zzz_xzzzz, g_y_0_x_0_zzz_yyyy, g_y_0_x_0_zzz_yyyyz, g_y_0_x_0_zzz_yyyz, g_y_0_x_0_zzz_yyyzz, g_y_0_x_0_zzz_yyzz, g_y_0_x_0_zzz_yyzzz, g_y_0_x_0_zzz_yzzz, g_y_0_x_0_zzz_yzzzz, g_y_0_x_0_zzz_zzzz, g_y_0_x_0_zzz_zzzzz, g_y_0_x_0_zzzz_xxxx, g_y_0_x_0_zzzz_xxxy, g_y_0_x_0_zzzz_xxxz, g_y_0_x_0_zzzz_xxyy, g_y_0_x_0_zzzz_xxyz, g_y_0_x_0_zzzz_xxzz, g_y_0_x_0_zzzz_xyyy, g_y_0_x_0_zzzz_xyyz, g_y_0_x_0_zzzz_xyzz, g_y_0_x_0_zzzz_xzzz, g_y_0_x_0_zzzz_yyyy, g_y_0_x_0_zzzz_yyyz, g_y_0_x_0_zzzz_yyzz, g_y_0_x_0_zzzz_yzzz, g_y_0_x_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_zzzz_xxxx[k] = -g_y_0_x_0_zzz_xxxx[k] * ab_z + g_y_0_x_0_zzz_xxxxz[k];

                g_y_0_x_0_zzzz_xxxy[k] = -g_y_0_x_0_zzz_xxxy[k] * ab_z + g_y_0_x_0_zzz_xxxyz[k];

                g_y_0_x_0_zzzz_xxxz[k] = -g_y_0_x_0_zzz_xxxz[k] * ab_z + g_y_0_x_0_zzz_xxxzz[k];

                g_y_0_x_0_zzzz_xxyy[k] = -g_y_0_x_0_zzz_xxyy[k] * ab_z + g_y_0_x_0_zzz_xxyyz[k];

                g_y_0_x_0_zzzz_xxyz[k] = -g_y_0_x_0_zzz_xxyz[k] * ab_z + g_y_0_x_0_zzz_xxyzz[k];

                g_y_0_x_0_zzzz_xxzz[k] = -g_y_0_x_0_zzz_xxzz[k] * ab_z + g_y_0_x_0_zzz_xxzzz[k];

                g_y_0_x_0_zzzz_xyyy[k] = -g_y_0_x_0_zzz_xyyy[k] * ab_z + g_y_0_x_0_zzz_xyyyz[k];

                g_y_0_x_0_zzzz_xyyz[k] = -g_y_0_x_0_zzz_xyyz[k] * ab_z + g_y_0_x_0_zzz_xyyzz[k];

                g_y_0_x_0_zzzz_xyzz[k] = -g_y_0_x_0_zzz_xyzz[k] * ab_z + g_y_0_x_0_zzz_xyzzz[k];

                g_y_0_x_0_zzzz_xzzz[k] = -g_y_0_x_0_zzz_xzzz[k] * ab_z + g_y_0_x_0_zzz_xzzzz[k];

                g_y_0_x_0_zzzz_yyyy[k] = -g_y_0_x_0_zzz_yyyy[k] * ab_z + g_y_0_x_0_zzz_yyyyz[k];

                g_y_0_x_0_zzzz_yyyz[k] = -g_y_0_x_0_zzz_yyyz[k] * ab_z + g_y_0_x_0_zzz_yyyzz[k];

                g_y_0_x_0_zzzz_yyzz[k] = -g_y_0_x_0_zzz_yyzz[k] * ab_z + g_y_0_x_0_zzz_yyzzz[k];

                g_y_0_x_0_zzzz_yzzz[k] = -g_y_0_x_0_zzz_yzzz[k] * ab_z + g_y_0_x_0_zzz_yzzzz[k];

                g_y_0_x_0_zzzz_zzzz[k] = -g_y_0_x_0_zzz_zzzz[k] * ab_z + g_y_0_x_0_zzz_zzzzz[k];
            }

            /// Set up 900-915 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxx_xxxx = cbuffer.data(gg_geom_1010_off + 900 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_xxxy = cbuffer.data(gg_geom_1010_off + 901 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_xxxz = cbuffer.data(gg_geom_1010_off + 902 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_xxyy = cbuffer.data(gg_geom_1010_off + 903 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_xxyz = cbuffer.data(gg_geom_1010_off + 904 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_xxzz = cbuffer.data(gg_geom_1010_off + 905 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_xyyy = cbuffer.data(gg_geom_1010_off + 906 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_xyyz = cbuffer.data(gg_geom_1010_off + 907 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_xyzz = cbuffer.data(gg_geom_1010_off + 908 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_xzzz = cbuffer.data(gg_geom_1010_off + 909 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_yyyy = cbuffer.data(gg_geom_1010_off + 910 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_yyyz = cbuffer.data(gg_geom_1010_off + 911 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_yyzz = cbuffer.data(gg_geom_1010_off + 912 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_yzzz = cbuffer.data(gg_geom_1010_off + 913 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_zzzz = cbuffer.data(gg_geom_1010_off + 914 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxx_xxxx, g_y_0_y_0_xxx_xxxxx, g_y_0_y_0_xxx_xxxxy, g_y_0_y_0_xxx_xxxxz, g_y_0_y_0_xxx_xxxy, g_y_0_y_0_xxx_xxxyy, g_y_0_y_0_xxx_xxxyz, g_y_0_y_0_xxx_xxxz, g_y_0_y_0_xxx_xxxzz, g_y_0_y_0_xxx_xxyy, g_y_0_y_0_xxx_xxyyy, g_y_0_y_0_xxx_xxyyz, g_y_0_y_0_xxx_xxyz, g_y_0_y_0_xxx_xxyzz, g_y_0_y_0_xxx_xxzz, g_y_0_y_0_xxx_xxzzz, g_y_0_y_0_xxx_xyyy, g_y_0_y_0_xxx_xyyyy, g_y_0_y_0_xxx_xyyyz, g_y_0_y_0_xxx_xyyz, g_y_0_y_0_xxx_xyyzz, g_y_0_y_0_xxx_xyzz, g_y_0_y_0_xxx_xyzzz, g_y_0_y_0_xxx_xzzz, g_y_0_y_0_xxx_xzzzz, g_y_0_y_0_xxx_yyyy, g_y_0_y_0_xxx_yyyz, g_y_0_y_0_xxx_yyzz, g_y_0_y_0_xxx_yzzz, g_y_0_y_0_xxx_zzzz, g_y_0_y_0_xxxx_xxxx, g_y_0_y_0_xxxx_xxxy, g_y_0_y_0_xxxx_xxxz, g_y_0_y_0_xxxx_xxyy, g_y_0_y_0_xxxx_xxyz, g_y_0_y_0_xxxx_xxzz, g_y_0_y_0_xxxx_xyyy, g_y_0_y_0_xxxx_xyyz, g_y_0_y_0_xxxx_xyzz, g_y_0_y_0_xxxx_xzzz, g_y_0_y_0_xxxx_yyyy, g_y_0_y_0_xxxx_yyyz, g_y_0_y_0_xxxx_yyzz, g_y_0_y_0_xxxx_yzzz, g_y_0_y_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxx_xxxx[k] = -g_y_0_y_0_xxx_xxxx[k] * ab_x + g_y_0_y_0_xxx_xxxxx[k];

                g_y_0_y_0_xxxx_xxxy[k] = -g_y_0_y_0_xxx_xxxy[k] * ab_x + g_y_0_y_0_xxx_xxxxy[k];

                g_y_0_y_0_xxxx_xxxz[k] = -g_y_0_y_0_xxx_xxxz[k] * ab_x + g_y_0_y_0_xxx_xxxxz[k];

                g_y_0_y_0_xxxx_xxyy[k] = -g_y_0_y_0_xxx_xxyy[k] * ab_x + g_y_0_y_0_xxx_xxxyy[k];

                g_y_0_y_0_xxxx_xxyz[k] = -g_y_0_y_0_xxx_xxyz[k] * ab_x + g_y_0_y_0_xxx_xxxyz[k];

                g_y_0_y_0_xxxx_xxzz[k] = -g_y_0_y_0_xxx_xxzz[k] * ab_x + g_y_0_y_0_xxx_xxxzz[k];

                g_y_0_y_0_xxxx_xyyy[k] = -g_y_0_y_0_xxx_xyyy[k] * ab_x + g_y_0_y_0_xxx_xxyyy[k];

                g_y_0_y_0_xxxx_xyyz[k] = -g_y_0_y_0_xxx_xyyz[k] * ab_x + g_y_0_y_0_xxx_xxyyz[k];

                g_y_0_y_0_xxxx_xyzz[k] = -g_y_0_y_0_xxx_xyzz[k] * ab_x + g_y_0_y_0_xxx_xxyzz[k];

                g_y_0_y_0_xxxx_xzzz[k] = -g_y_0_y_0_xxx_xzzz[k] * ab_x + g_y_0_y_0_xxx_xxzzz[k];

                g_y_0_y_0_xxxx_yyyy[k] = -g_y_0_y_0_xxx_yyyy[k] * ab_x + g_y_0_y_0_xxx_xyyyy[k];

                g_y_0_y_0_xxxx_yyyz[k] = -g_y_0_y_0_xxx_yyyz[k] * ab_x + g_y_0_y_0_xxx_xyyyz[k];

                g_y_0_y_0_xxxx_yyzz[k] = -g_y_0_y_0_xxx_yyzz[k] * ab_x + g_y_0_y_0_xxx_xyyzz[k];

                g_y_0_y_0_xxxx_yzzz[k] = -g_y_0_y_0_xxx_yzzz[k] * ab_x + g_y_0_y_0_xxx_xyzzz[k];

                g_y_0_y_0_xxxx_zzzz[k] = -g_y_0_y_0_xxx_zzzz[k] * ab_x + g_y_0_y_0_xxx_xzzzz[k];
            }

            /// Set up 915-930 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxy_xxxx = cbuffer.data(gg_geom_1010_off + 915 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_xxxy = cbuffer.data(gg_geom_1010_off + 916 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_xxxz = cbuffer.data(gg_geom_1010_off + 917 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_xxyy = cbuffer.data(gg_geom_1010_off + 918 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_xxyz = cbuffer.data(gg_geom_1010_off + 919 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_xxzz = cbuffer.data(gg_geom_1010_off + 920 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_xyyy = cbuffer.data(gg_geom_1010_off + 921 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_xyyz = cbuffer.data(gg_geom_1010_off + 922 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_xyzz = cbuffer.data(gg_geom_1010_off + 923 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_xzzz = cbuffer.data(gg_geom_1010_off + 924 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_yyyy = cbuffer.data(gg_geom_1010_off + 925 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_yyyz = cbuffer.data(gg_geom_1010_off + 926 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_yyzz = cbuffer.data(gg_geom_1010_off + 927 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_yzzz = cbuffer.data(gg_geom_1010_off + 928 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_zzzz = cbuffer.data(gg_geom_1010_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxy_xxxx, g_y_0_y_0_xxxy_xxxy, g_y_0_y_0_xxxy_xxxz, g_y_0_y_0_xxxy_xxyy, g_y_0_y_0_xxxy_xxyz, g_y_0_y_0_xxxy_xxzz, g_y_0_y_0_xxxy_xyyy, g_y_0_y_0_xxxy_xyyz, g_y_0_y_0_xxxy_xyzz, g_y_0_y_0_xxxy_xzzz, g_y_0_y_0_xxxy_yyyy, g_y_0_y_0_xxxy_yyyz, g_y_0_y_0_xxxy_yyzz, g_y_0_y_0_xxxy_yzzz, g_y_0_y_0_xxxy_zzzz, g_y_0_y_0_xxy_xxxx, g_y_0_y_0_xxy_xxxxx, g_y_0_y_0_xxy_xxxxy, g_y_0_y_0_xxy_xxxxz, g_y_0_y_0_xxy_xxxy, g_y_0_y_0_xxy_xxxyy, g_y_0_y_0_xxy_xxxyz, g_y_0_y_0_xxy_xxxz, g_y_0_y_0_xxy_xxxzz, g_y_0_y_0_xxy_xxyy, g_y_0_y_0_xxy_xxyyy, g_y_0_y_0_xxy_xxyyz, g_y_0_y_0_xxy_xxyz, g_y_0_y_0_xxy_xxyzz, g_y_0_y_0_xxy_xxzz, g_y_0_y_0_xxy_xxzzz, g_y_0_y_0_xxy_xyyy, g_y_0_y_0_xxy_xyyyy, g_y_0_y_0_xxy_xyyyz, g_y_0_y_0_xxy_xyyz, g_y_0_y_0_xxy_xyyzz, g_y_0_y_0_xxy_xyzz, g_y_0_y_0_xxy_xyzzz, g_y_0_y_0_xxy_xzzz, g_y_0_y_0_xxy_xzzzz, g_y_0_y_0_xxy_yyyy, g_y_0_y_0_xxy_yyyz, g_y_0_y_0_xxy_yyzz, g_y_0_y_0_xxy_yzzz, g_y_0_y_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxy_xxxx[k] = -g_y_0_y_0_xxy_xxxx[k] * ab_x + g_y_0_y_0_xxy_xxxxx[k];

                g_y_0_y_0_xxxy_xxxy[k] = -g_y_0_y_0_xxy_xxxy[k] * ab_x + g_y_0_y_0_xxy_xxxxy[k];

                g_y_0_y_0_xxxy_xxxz[k] = -g_y_0_y_0_xxy_xxxz[k] * ab_x + g_y_0_y_0_xxy_xxxxz[k];

                g_y_0_y_0_xxxy_xxyy[k] = -g_y_0_y_0_xxy_xxyy[k] * ab_x + g_y_0_y_0_xxy_xxxyy[k];

                g_y_0_y_0_xxxy_xxyz[k] = -g_y_0_y_0_xxy_xxyz[k] * ab_x + g_y_0_y_0_xxy_xxxyz[k];

                g_y_0_y_0_xxxy_xxzz[k] = -g_y_0_y_0_xxy_xxzz[k] * ab_x + g_y_0_y_0_xxy_xxxzz[k];

                g_y_0_y_0_xxxy_xyyy[k] = -g_y_0_y_0_xxy_xyyy[k] * ab_x + g_y_0_y_0_xxy_xxyyy[k];

                g_y_0_y_0_xxxy_xyyz[k] = -g_y_0_y_0_xxy_xyyz[k] * ab_x + g_y_0_y_0_xxy_xxyyz[k];

                g_y_0_y_0_xxxy_xyzz[k] = -g_y_0_y_0_xxy_xyzz[k] * ab_x + g_y_0_y_0_xxy_xxyzz[k];

                g_y_0_y_0_xxxy_xzzz[k] = -g_y_0_y_0_xxy_xzzz[k] * ab_x + g_y_0_y_0_xxy_xxzzz[k];

                g_y_0_y_0_xxxy_yyyy[k] = -g_y_0_y_0_xxy_yyyy[k] * ab_x + g_y_0_y_0_xxy_xyyyy[k];

                g_y_0_y_0_xxxy_yyyz[k] = -g_y_0_y_0_xxy_yyyz[k] * ab_x + g_y_0_y_0_xxy_xyyyz[k];

                g_y_0_y_0_xxxy_yyzz[k] = -g_y_0_y_0_xxy_yyzz[k] * ab_x + g_y_0_y_0_xxy_xyyzz[k];

                g_y_0_y_0_xxxy_yzzz[k] = -g_y_0_y_0_xxy_yzzz[k] * ab_x + g_y_0_y_0_xxy_xyzzz[k];

                g_y_0_y_0_xxxy_zzzz[k] = -g_y_0_y_0_xxy_zzzz[k] * ab_x + g_y_0_y_0_xxy_xzzzz[k];
            }

            /// Set up 930-945 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxz_xxxx = cbuffer.data(gg_geom_1010_off + 930 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_xxxy = cbuffer.data(gg_geom_1010_off + 931 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_xxxz = cbuffer.data(gg_geom_1010_off + 932 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_xxyy = cbuffer.data(gg_geom_1010_off + 933 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_xxyz = cbuffer.data(gg_geom_1010_off + 934 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_xxzz = cbuffer.data(gg_geom_1010_off + 935 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_xyyy = cbuffer.data(gg_geom_1010_off + 936 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_xyyz = cbuffer.data(gg_geom_1010_off + 937 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_xyzz = cbuffer.data(gg_geom_1010_off + 938 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_xzzz = cbuffer.data(gg_geom_1010_off + 939 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_yyyy = cbuffer.data(gg_geom_1010_off + 940 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_yyyz = cbuffer.data(gg_geom_1010_off + 941 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_yyzz = cbuffer.data(gg_geom_1010_off + 942 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_yzzz = cbuffer.data(gg_geom_1010_off + 943 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_zzzz = cbuffer.data(gg_geom_1010_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxz_xxxx, g_y_0_y_0_xxxz_xxxy, g_y_0_y_0_xxxz_xxxz, g_y_0_y_0_xxxz_xxyy, g_y_0_y_0_xxxz_xxyz, g_y_0_y_0_xxxz_xxzz, g_y_0_y_0_xxxz_xyyy, g_y_0_y_0_xxxz_xyyz, g_y_0_y_0_xxxz_xyzz, g_y_0_y_0_xxxz_xzzz, g_y_0_y_0_xxxz_yyyy, g_y_0_y_0_xxxz_yyyz, g_y_0_y_0_xxxz_yyzz, g_y_0_y_0_xxxz_yzzz, g_y_0_y_0_xxxz_zzzz, g_y_0_y_0_xxz_xxxx, g_y_0_y_0_xxz_xxxxx, g_y_0_y_0_xxz_xxxxy, g_y_0_y_0_xxz_xxxxz, g_y_0_y_0_xxz_xxxy, g_y_0_y_0_xxz_xxxyy, g_y_0_y_0_xxz_xxxyz, g_y_0_y_0_xxz_xxxz, g_y_0_y_0_xxz_xxxzz, g_y_0_y_0_xxz_xxyy, g_y_0_y_0_xxz_xxyyy, g_y_0_y_0_xxz_xxyyz, g_y_0_y_0_xxz_xxyz, g_y_0_y_0_xxz_xxyzz, g_y_0_y_0_xxz_xxzz, g_y_0_y_0_xxz_xxzzz, g_y_0_y_0_xxz_xyyy, g_y_0_y_0_xxz_xyyyy, g_y_0_y_0_xxz_xyyyz, g_y_0_y_0_xxz_xyyz, g_y_0_y_0_xxz_xyyzz, g_y_0_y_0_xxz_xyzz, g_y_0_y_0_xxz_xyzzz, g_y_0_y_0_xxz_xzzz, g_y_0_y_0_xxz_xzzzz, g_y_0_y_0_xxz_yyyy, g_y_0_y_0_xxz_yyyz, g_y_0_y_0_xxz_yyzz, g_y_0_y_0_xxz_yzzz, g_y_0_y_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxz_xxxx[k] = -g_y_0_y_0_xxz_xxxx[k] * ab_x + g_y_0_y_0_xxz_xxxxx[k];

                g_y_0_y_0_xxxz_xxxy[k] = -g_y_0_y_0_xxz_xxxy[k] * ab_x + g_y_0_y_0_xxz_xxxxy[k];

                g_y_0_y_0_xxxz_xxxz[k] = -g_y_0_y_0_xxz_xxxz[k] * ab_x + g_y_0_y_0_xxz_xxxxz[k];

                g_y_0_y_0_xxxz_xxyy[k] = -g_y_0_y_0_xxz_xxyy[k] * ab_x + g_y_0_y_0_xxz_xxxyy[k];

                g_y_0_y_0_xxxz_xxyz[k] = -g_y_0_y_0_xxz_xxyz[k] * ab_x + g_y_0_y_0_xxz_xxxyz[k];

                g_y_0_y_0_xxxz_xxzz[k] = -g_y_0_y_0_xxz_xxzz[k] * ab_x + g_y_0_y_0_xxz_xxxzz[k];

                g_y_0_y_0_xxxz_xyyy[k] = -g_y_0_y_0_xxz_xyyy[k] * ab_x + g_y_0_y_0_xxz_xxyyy[k];

                g_y_0_y_0_xxxz_xyyz[k] = -g_y_0_y_0_xxz_xyyz[k] * ab_x + g_y_0_y_0_xxz_xxyyz[k];

                g_y_0_y_0_xxxz_xyzz[k] = -g_y_0_y_0_xxz_xyzz[k] * ab_x + g_y_0_y_0_xxz_xxyzz[k];

                g_y_0_y_0_xxxz_xzzz[k] = -g_y_0_y_0_xxz_xzzz[k] * ab_x + g_y_0_y_0_xxz_xxzzz[k];

                g_y_0_y_0_xxxz_yyyy[k] = -g_y_0_y_0_xxz_yyyy[k] * ab_x + g_y_0_y_0_xxz_xyyyy[k];

                g_y_0_y_0_xxxz_yyyz[k] = -g_y_0_y_0_xxz_yyyz[k] * ab_x + g_y_0_y_0_xxz_xyyyz[k];

                g_y_0_y_0_xxxz_yyzz[k] = -g_y_0_y_0_xxz_yyzz[k] * ab_x + g_y_0_y_0_xxz_xyyzz[k];

                g_y_0_y_0_xxxz_yzzz[k] = -g_y_0_y_0_xxz_yzzz[k] * ab_x + g_y_0_y_0_xxz_xyzzz[k];

                g_y_0_y_0_xxxz_zzzz[k] = -g_y_0_y_0_xxz_zzzz[k] * ab_x + g_y_0_y_0_xxz_xzzzz[k];
            }

            /// Set up 945-960 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyy_xxxx = cbuffer.data(gg_geom_1010_off + 945 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_xxxy = cbuffer.data(gg_geom_1010_off + 946 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_xxxz = cbuffer.data(gg_geom_1010_off + 947 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_xxyy = cbuffer.data(gg_geom_1010_off + 948 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_xxyz = cbuffer.data(gg_geom_1010_off + 949 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_xxzz = cbuffer.data(gg_geom_1010_off + 950 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_xyyy = cbuffer.data(gg_geom_1010_off + 951 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_xyyz = cbuffer.data(gg_geom_1010_off + 952 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_xyzz = cbuffer.data(gg_geom_1010_off + 953 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_xzzz = cbuffer.data(gg_geom_1010_off + 954 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_yyyy = cbuffer.data(gg_geom_1010_off + 955 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_yyyz = cbuffer.data(gg_geom_1010_off + 956 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_yyzz = cbuffer.data(gg_geom_1010_off + 957 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_yzzz = cbuffer.data(gg_geom_1010_off + 958 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_zzzz = cbuffer.data(gg_geom_1010_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyy_xxxx, g_y_0_y_0_xxyy_xxxy, g_y_0_y_0_xxyy_xxxz, g_y_0_y_0_xxyy_xxyy, g_y_0_y_0_xxyy_xxyz, g_y_0_y_0_xxyy_xxzz, g_y_0_y_0_xxyy_xyyy, g_y_0_y_0_xxyy_xyyz, g_y_0_y_0_xxyy_xyzz, g_y_0_y_0_xxyy_xzzz, g_y_0_y_0_xxyy_yyyy, g_y_0_y_0_xxyy_yyyz, g_y_0_y_0_xxyy_yyzz, g_y_0_y_0_xxyy_yzzz, g_y_0_y_0_xxyy_zzzz, g_y_0_y_0_xyy_xxxx, g_y_0_y_0_xyy_xxxxx, g_y_0_y_0_xyy_xxxxy, g_y_0_y_0_xyy_xxxxz, g_y_0_y_0_xyy_xxxy, g_y_0_y_0_xyy_xxxyy, g_y_0_y_0_xyy_xxxyz, g_y_0_y_0_xyy_xxxz, g_y_0_y_0_xyy_xxxzz, g_y_0_y_0_xyy_xxyy, g_y_0_y_0_xyy_xxyyy, g_y_0_y_0_xyy_xxyyz, g_y_0_y_0_xyy_xxyz, g_y_0_y_0_xyy_xxyzz, g_y_0_y_0_xyy_xxzz, g_y_0_y_0_xyy_xxzzz, g_y_0_y_0_xyy_xyyy, g_y_0_y_0_xyy_xyyyy, g_y_0_y_0_xyy_xyyyz, g_y_0_y_0_xyy_xyyz, g_y_0_y_0_xyy_xyyzz, g_y_0_y_0_xyy_xyzz, g_y_0_y_0_xyy_xyzzz, g_y_0_y_0_xyy_xzzz, g_y_0_y_0_xyy_xzzzz, g_y_0_y_0_xyy_yyyy, g_y_0_y_0_xyy_yyyz, g_y_0_y_0_xyy_yyzz, g_y_0_y_0_xyy_yzzz, g_y_0_y_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyy_xxxx[k] = -g_y_0_y_0_xyy_xxxx[k] * ab_x + g_y_0_y_0_xyy_xxxxx[k];

                g_y_0_y_0_xxyy_xxxy[k] = -g_y_0_y_0_xyy_xxxy[k] * ab_x + g_y_0_y_0_xyy_xxxxy[k];

                g_y_0_y_0_xxyy_xxxz[k] = -g_y_0_y_0_xyy_xxxz[k] * ab_x + g_y_0_y_0_xyy_xxxxz[k];

                g_y_0_y_0_xxyy_xxyy[k] = -g_y_0_y_0_xyy_xxyy[k] * ab_x + g_y_0_y_0_xyy_xxxyy[k];

                g_y_0_y_0_xxyy_xxyz[k] = -g_y_0_y_0_xyy_xxyz[k] * ab_x + g_y_0_y_0_xyy_xxxyz[k];

                g_y_0_y_0_xxyy_xxzz[k] = -g_y_0_y_0_xyy_xxzz[k] * ab_x + g_y_0_y_0_xyy_xxxzz[k];

                g_y_0_y_0_xxyy_xyyy[k] = -g_y_0_y_0_xyy_xyyy[k] * ab_x + g_y_0_y_0_xyy_xxyyy[k];

                g_y_0_y_0_xxyy_xyyz[k] = -g_y_0_y_0_xyy_xyyz[k] * ab_x + g_y_0_y_0_xyy_xxyyz[k];

                g_y_0_y_0_xxyy_xyzz[k] = -g_y_0_y_0_xyy_xyzz[k] * ab_x + g_y_0_y_0_xyy_xxyzz[k];

                g_y_0_y_0_xxyy_xzzz[k] = -g_y_0_y_0_xyy_xzzz[k] * ab_x + g_y_0_y_0_xyy_xxzzz[k];

                g_y_0_y_0_xxyy_yyyy[k] = -g_y_0_y_0_xyy_yyyy[k] * ab_x + g_y_0_y_0_xyy_xyyyy[k];

                g_y_0_y_0_xxyy_yyyz[k] = -g_y_0_y_0_xyy_yyyz[k] * ab_x + g_y_0_y_0_xyy_xyyyz[k];

                g_y_0_y_0_xxyy_yyzz[k] = -g_y_0_y_0_xyy_yyzz[k] * ab_x + g_y_0_y_0_xyy_xyyzz[k];

                g_y_0_y_0_xxyy_yzzz[k] = -g_y_0_y_0_xyy_yzzz[k] * ab_x + g_y_0_y_0_xyy_xyzzz[k];

                g_y_0_y_0_xxyy_zzzz[k] = -g_y_0_y_0_xyy_zzzz[k] * ab_x + g_y_0_y_0_xyy_xzzzz[k];
            }

            /// Set up 960-975 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyz_xxxx = cbuffer.data(gg_geom_1010_off + 960 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_xxxy = cbuffer.data(gg_geom_1010_off + 961 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_xxxz = cbuffer.data(gg_geom_1010_off + 962 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_xxyy = cbuffer.data(gg_geom_1010_off + 963 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_xxyz = cbuffer.data(gg_geom_1010_off + 964 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_xxzz = cbuffer.data(gg_geom_1010_off + 965 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_xyyy = cbuffer.data(gg_geom_1010_off + 966 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_xyyz = cbuffer.data(gg_geom_1010_off + 967 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_xyzz = cbuffer.data(gg_geom_1010_off + 968 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_xzzz = cbuffer.data(gg_geom_1010_off + 969 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_yyyy = cbuffer.data(gg_geom_1010_off + 970 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_yyyz = cbuffer.data(gg_geom_1010_off + 971 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_yyzz = cbuffer.data(gg_geom_1010_off + 972 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_yzzz = cbuffer.data(gg_geom_1010_off + 973 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_zzzz = cbuffer.data(gg_geom_1010_off + 974 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyz_xxxx, g_y_0_y_0_xxyz_xxxy, g_y_0_y_0_xxyz_xxxz, g_y_0_y_0_xxyz_xxyy, g_y_0_y_0_xxyz_xxyz, g_y_0_y_0_xxyz_xxzz, g_y_0_y_0_xxyz_xyyy, g_y_0_y_0_xxyz_xyyz, g_y_0_y_0_xxyz_xyzz, g_y_0_y_0_xxyz_xzzz, g_y_0_y_0_xxyz_yyyy, g_y_0_y_0_xxyz_yyyz, g_y_0_y_0_xxyz_yyzz, g_y_0_y_0_xxyz_yzzz, g_y_0_y_0_xxyz_zzzz, g_y_0_y_0_xyz_xxxx, g_y_0_y_0_xyz_xxxxx, g_y_0_y_0_xyz_xxxxy, g_y_0_y_0_xyz_xxxxz, g_y_0_y_0_xyz_xxxy, g_y_0_y_0_xyz_xxxyy, g_y_0_y_0_xyz_xxxyz, g_y_0_y_0_xyz_xxxz, g_y_0_y_0_xyz_xxxzz, g_y_0_y_0_xyz_xxyy, g_y_0_y_0_xyz_xxyyy, g_y_0_y_0_xyz_xxyyz, g_y_0_y_0_xyz_xxyz, g_y_0_y_0_xyz_xxyzz, g_y_0_y_0_xyz_xxzz, g_y_0_y_0_xyz_xxzzz, g_y_0_y_0_xyz_xyyy, g_y_0_y_0_xyz_xyyyy, g_y_0_y_0_xyz_xyyyz, g_y_0_y_0_xyz_xyyz, g_y_0_y_0_xyz_xyyzz, g_y_0_y_0_xyz_xyzz, g_y_0_y_0_xyz_xyzzz, g_y_0_y_0_xyz_xzzz, g_y_0_y_0_xyz_xzzzz, g_y_0_y_0_xyz_yyyy, g_y_0_y_0_xyz_yyyz, g_y_0_y_0_xyz_yyzz, g_y_0_y_0_xyz_yzzz, g_y_0_y_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyz_xxxx[k] = -g_y_0_y_0_xyz_xxxx[k] * ab_x + g_y_0_y_0_xyz_xxxxx[k];

                g_y_0_y_0_xxyz_xxxy[k] = -g_y_0_y_0_xyz_xxxy[k] * ab_x + g_y_0_y_0_xyz_xxxxy[k];

                g_y_0_y_0_xxyz_xxxz[k] = -g_y_0_y_0_xyz_xxxz[k] * ab_x + g_y_0_y_0_xyz_xxxxz[k];

                g_y_0_y_0_xxyz_xxyy[k] = -g_y_0_y_0_xyz_xxyy[k] * ab_x + g_y_0_y_0_xyz_xxxyy[k];

                g_y_0_y_0_xxyz_xxyz[k] = -g_y_0_y_0_xyz_xxyz[k] * ab_x + g_y_0_y_0_xyz_xxxyz[k];

                g_y_0_y_0_xxyz_xxzz[k] = -g_y_0_y_0_xyz_xxzz[k] * ab_x + g_y_0_y_0_xyz_xxxzz[k];

                g_y_0_y_0_xxyz_xyyy[k] = -g_y_0_y_0_xyz_xyyy[k] * ab_x + g_y_0_y_0_xyz_xxyyy[k];

                g_y_0_y_0_xxyz_xyyz[k] = -g_y_0_y_0_xyz_xyyz[k] * ab_x + g_y_0_y_0_xyz_xxyyz[k];

                g_y_0_y_0_xxyz_xyzz[k] = -g_y_0_y_0_xyz_xyzz[k] * ab_x + g_y_0_y_0_xyz_xxyzz[k];

                g_y_0_y_0_xxyz_xzzz[k] = -g_y_0_y_0_xyz_xzzz[k] * ab_x + g_y_0_y_0_xyz_xxzzz[k];

                g_y_0_y_0_xxyz_yyyy[k] = -g_y_0_y_0_xyz_yyyy[k] * ab_x + g_y_0_y_0_xyz_xyyyy[k];

                g_y_0_y_0_xxyz_yyyz[k] = -g_y_0_y_0_xyz_yyyz[k] * ab_x + g_y_0_y_0_xyz_xyyyz[k];

                g_y_0_y_0_xxyz_yyzz[k] = -g_y_0_y_0_xyz_yyzz[k] * ab_x + g_y_0_y_0_xyz_xyyzz[k];

                g_y_0_y_0_xxyz_yzzz[k] = -g_y_0_y_0_xyz_yzzz[k] * ab_x + g_y_0_y_0_xyz_xyzzz[k];

                g_y_0_y_0_xxyz_zzzz[k] = -g_y_0_y_0_xyz_zzzz[k] * ab_x + g_y_0_y_0_xyz_xzzzz[k];
            }

            /// Set up 975-990 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxzz_xxxx = cbuffer.data(gg_geom_1010_off + 975 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_xxxy = cbuffer.data(gg_geom_1010_off + 976 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_xxxz = cbuffer.data(gg_geom_1010_off + 977 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_xxyy = cbuffer.data(gg_geom_1010_off + 978 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_xxyz = cbuffer.data(gg_geom_1010_off + 979 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_xxzz = cbuffer.data(gg_geom_1010_off + 980 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_xyyy = cbuffer.data(gg_geom_1010_off + 981 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_xyyz = cbuffer.data(gg_geom_1010_off + 982 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_xyzz = cbuffer.data(gg_geom_1010_off + 983 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_xzzz = cbuffer.data(gg_geom_1010_off + 984 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_yyyy = cbuffer.data(gg_geom_1010_off + 985 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_yyyz = cbuffer.data(gg_geom_1010_off + 986 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_yyzz = cbuffer.data(gg_geom_1010_off + 987 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_yzzz = cbuffer.data(gg_geom_1010_off + 988 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_zzzz = cbuffer.data(gg_geom_1010_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxzz_xxxx, g_y_0_y_0_xxzz_xxxy, g_y_0_y_0_xxzz_xxxz, g_y_0_y_0_xxzz_xxyy, g_y_0_y_0_xxzz_xxyz, g_y_0_y_0_xxzz_xxzz, g_y_0_y_0_xxzz_xyyy, g_y_0_y_0_xxzz_xyyz, g_y_0_y_0_xxzz_xyzz, g_y_0_y_0_xxzz_xzzz, g_y_0_y_0_xxzz_yyyy, g_y_0_y_0_xxzz_yyyz, g_y_0_y_0_xxzz_yyzz, g_y_0_y_0_xxzz_yzzz, g_y_0_y_0_xxzz_zzzz, g_y_0_y_0_xzz_xxxx, g_y_0_y_0_xzz_xxxxx, g_y_0_y_0_xzz_xxxxy, g_y_0_y_0_xzz_xxxxz, g_y_0_y_0_xzz_xxxy, g_y_0_y_0_xzz_xxxyy, g_y_0_y_0_xzz_xxxyz, g_y_0_y_0_xzz_xxxz, g_y_0_y_0_xzz_xxxzz, g_y_0_y_0_xzz_xxyy, g_y_0_y_0_xzz_xxyyy, g_y_0_y_0_xzz_xxyyz, g_y_0_y_0_xzz_xxyz, g_y_0_y_0_xzz_xxyzz, g_y_0_y_0_xzz_xxzz, g_y_0_y_0_xzz_xxzzz, g_y_0_y_0_xzz_xyyy, g_y_0_y_0_xzz_xyyyy, g_y_0_y_0_xzz_xyyyz, g_y_0_y_0_xzz_xyyz, g_y_0_y_0_xzz_xyyzz, g_y_0_y_0_xzz_xyzz, g_y_0_y_0_xzz_xyzzz, g_y_0_y_0_xzz_xzzz, g_y_0_y_0_xzz_xzzzz, g_y_0_y_0_xzz_yyyy, g_y_0_y_0_xzz_yyyz, g_y_0_y_0_xzz_yyzz, g_y_0_y_0_xzz_yzzz, g_y_0_y_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxzz_xxxx[k] = -g_y_0_y_0_xzz_xxxx[k] * ab_x + g_y_0_y_0_xzz_xxxxx[k];

                g_y_0_y_0_xxzz_xxxy[k] = -g_y_0_y_0_xzz_xxxy[k] * ab_x + g_y_0_y_0_xzz_xxxxy[k];

                g_y_0_y_0_xxzz_xxxz[k] = -g_y_0_y_0_xzz_xxxz[k] * ab_x + g_y_0_y_0_xzz_xxxxz[k];

                g_y_0_y_0_xxzz_xxyy[k] = -g_y_0_y_0_xzz_xxyy[k] * ab_x + g_y_0_y_0_xzz_xxxyy[k];

                g_y_0_y_0_xxzz_xxyz[k] = -g_y_0_y_0_xzz_xxyz[k] * ab_x + g_y_0_y_0_xzz_xxxyz[k];

                g_y_0_y_0_xxzz_xxzz[k] = -g_y_0_y_0_xzz_xxzz[k] * ab_x + g_y_0_y_0_xzz_xxxzz[k];

                g_y_0_y_0_xxzz_xyyy[k] = -g_y_0_y_0_xzz_xyyy[k] * ab_x + g_y_0_y_0_xzz_xxyyy[k];

                g_y_0_y_0_xxzz_xyyz[k] = -g_y_0_y_0_xzz_xyyz[k] * ab_x + g_y_0_y_0_xzz_xxyyz[k];

                g_y_0_y_0_xxzz_xyzz[k] = -g_y_0_y_0_xzz_xyzz[k] * ab_x + g_y_0_y_0_xzz_xxyzz[k];

                g_y_0_y_0_xxzz_xzzz[k] = -g_y_0_y_0_xzz_xzzz[k] * ab_x + g_y_0_y_0_xzz_xxzzz[k];

                g_y_0_y_0_xxzz_yyyy[k] = -g_y_0_y_0_xzz_yyyy[k] * ab_x + g_y_0_y_0_xzz_xyyyy[k];

                g_y_0_y_0_xxzz_yyyz[k] = -g_y_0_y_0_xzz_yyyz[k] * ab_x + g_y_0_y_0_xzz_xyyyz[k];

                g_y_0_y_0_xxzz_yyzz[k] = -g_y_0_y_0_xzz_yyzz[k] * ab_x + g_y_0_y_0_xzz_xyyzz[k];

                g_y_0_y_0_xxzz_yzzz[k] = -g_y_0_y_0_xzz_yzzz[k] * ab_x + g_y_0_y_0_xzz_xyzzz[k];

                g_y_0_y_0_xxzz_zzzz[k] = -g_y_0_y_0_xzz_zzzz[k] * ab_x + g_y_0_y_0_xzz_xzzzz[k];
            }

            /// Set up 990-1005 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyy_xxxx = cbuffer.data(gg_geom_1010_off + 990 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_xxxy = cbuffer.data(gg_geom_1010_off + 991 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_xxxz = cbuffer.data(gg_geom_1010_off + 992 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_xxyy = cbuffer.data(gg_geom_1010_off + 993 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_xxyz = cbuffer.data(gg_geom_1010_off + 994 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_xxzz = cbuffer.data(gg_geom_1010_off + 995 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_xyyy = cbuffer.data(gg_geom_1010_off + 996 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_xyyz = cbuffer.data(gg_geom_1010_off + 997 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_xyzz = cbuffer.data(gg_geom_1010_off + 998 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_xzzz = cbuffer.data(gg_geom_1010_off + 999 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1000 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1001 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1002 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1003 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1004 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyy_xxxx, g_y_0_y_0_xyyy_xxxy, g_y_0_y_0_xyyy_xxxz, g_y_0_y_0_xyyy_xxyy, g_y_0_y_0_xyyy_xxyz, g_y_0_y_0_xyyy_xxzz, g_y_0_y_0_xyyy_xyyy, g_y_0_y_0_xyyy_xyyz, g_y_0_y_0_xyyy_xyzz, g_y_0_y_0_xyyy_xzzz, g_y_0_y_0_xyyy_yyyy, g_y_0_y_0_xyyy_yyyz, g_y_0_y_0_xyyy_yyzz, g_y_0_y_0_xyyy_yzzz, g_y_0_y_0_xyyy_zzzz, g_y_0_y_0_yyy_xxxx, g_y_0_y_0_yyy_xxxxx, g_y_0_y_0_yyy_xxxxy, g_y_0_y_0_yyy_xxxxz, g_y_0_y_0_yyy_xxxy, g_y_0_y_0_yyy_xxxyy, g_y_0_y_0_yyy_xxxyz, g_y_0_y_0_yyy_xxxz, g_y_0_y_0_yyy_xxxzz, g_y_0_y_0_yyy_xxyy, g_y_0_y_0_yyy_xxyyy, g_y_0_y_0_yyy_xxyyz, g_y_0_y_0_yyy_xxyz, g_y_0_y_0_yyy_xxyzz, g_y_0_y_0_yyy_xxzz, g_y_0_y_0_yyy_xxzzz, g_y_0_y_0_yyy_xyyy, g_y_0_y_0_yyy_xyyyy, g_y_0_y_0_yyy_xyyyz, g_y_0_y_0_yyy_xyyz, g_y_0_y_0_yyy_xyyzz, g_y_0_y_0_yyy_xyzz, g_y_0_y_0_yyy_xyzzz, g_y_0_y_0_yyy_xzzz, g_y_0_y_0_yyy_xzzzz, g_y_0_y_0_yyy_yyyy, g_y_0_y_0_yyy_yyyz, g_y_0_y_0_yyy_yyzz, g_y_0_y_0_yyy_yzzz, g_y_0_y_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyy_xxxx[k] = -g_y_0_y_0_yyy_xxxx[k] * ab_x + g_y_0_y_0_yyy_xxxxx[k];

                g_y_0_y_0_xyyy_xxxy[k] = -g_y_0_y_0_yyy_xxxy[k] * ab_x + g_y_0_y_0_yyy_xxxxy[k];

                g_y_0_y_0_xyyy_xxxz[k] = -g_y_0_y_0_yyy_xxxz[k] * ab_x + g_y_0_y_0_yyy_xxxxz[k];

                g_y_0_y_0_xyyy_xxyy[k] = -g_y_0_y_0_yyy_xxyy[k] * ab_x + g_y_0_y_0_yyy_xxxyy[k];

                g_y_0_y_0_xyyy_xxyz[k] = -g_y_0_y_0_yyy_xxyz[k] * ab_x + g_y_0_y_0_yyy_xxxyz[k];

                g_y_0_y_0_xyyy_xxzz[k] = -g_y_0_y_0_yyy_xxzz[k] * ab_x + g_y_0_y_0_yyy_xxxzz[k];

                g_y_0_y_0_xyyy_xyyy[k] = -g_y_0_y_0_yyy_xyyy[k] * ab_x + g_y_0_y_0_yyy_xxyyy[k];

                g_y_0_y_0_xyyy_xyyz[k] = -g_y_0_y_0_yyy_xyyz[k] * ab_x + g_y_0_y_0_yyy_xxyyz[k];

                g_y_0_y_0_xyyy_xyzz[k] = -g_y_0_y_0_yyy_xyzz[k] * ab_x + g_y_0_y_0_yyy_xxyzz[k];

                g_y_0_y_0_xyyy_xzzz[k] = -g_y_0_y_0_yyy_xzzz[k] * ab_x + g_y_0_y_0_yyy_xxzzz[k];

                g_y_0_y_0_xyyy_yyyy[k] = -g_y_0_y_0_yyy_yyyy[k] * ab_x + g_y_0_y_0_yyy_xyyyy[k];

                g_y_0_y_0_xyyy_yyyz[k] = -g_y_0_y_0_yyy_yyyz[k] * ab_x + g_y_0_y_0_yyy_xyyyz[k];

                g_y_0_y_0_xyyy_yyzz[k] = -g_y_0_y_0_yyy_yyzz[k] * ab_x + g_y_0_y_0_yyy_xyyzz[k];

                g_y_0_y_0_xyyy_yzzz[k] = -g_y_0_y_0_yyy_yzzz[k] * ab_x + g_y_0_y_0_yyy_xyzzz[k];

                g_y_0_y_0_xyyy_zzzz[k] = -g_y_0_y_0_yyy_zzzz[k] * ab_x + g_y_0_y_0_yyy_xzzzz[k];
            }

            /// Set up 1005-1020 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1005 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1006 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1007 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1008 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1009 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1010 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1011 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1012 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1013 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1014 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1015 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1016 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1017 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1018 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1019 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyz_xxxx, g_y_0_y_0_xyyz_xxxy, g_y_0_y_0_xyyz_xxxz, g_y_0_y_0_xyyz_xxyy, g_y_0_y_0_xyyz_xxyz, g_y_0_y_0_xyyz_xxzz, g_y_0_y_0_xyyz_xyyy, g_y_0_y_0_xyyz_xyyz, g_y_0_y_0_xyyz_xyzz, g_y_0_y_0_xyyz_xzzz, g_y_0_y_0_xyyz_yyyy, g_y_0_y_0_xyyz_yyyz, g_y_0_y_0_xyyz_yyzz, g_y_0_y_0_xyyz_yzzz, g_y_0_y_0_xyyz_zzzz, g_y_0_y_0_yyz_xxxx, g_y_0_y_0_yyz_xxxxx, g_y_0_y_0_yyz_xxxxy, g_y_0_y_0_yyz_xxxxz, g_y_0_y_0_yyz_xxxy, g_y_0_y_0_yyz_xxxyy, g_y_0_y_0_yyz_xxxyz, g_y_0_y_0_yyz_xxxz, g_y_0_y_0_yyz_xxxzz, g_y_0_y_0_yyz_xxyy, g_y_0_y_0_yyz_xxyyy, g_y_0_y_0_yyz_xxyyz, g_y_0_y_0_yyz_xxyz, g_y_0_y_0_yyz_xxyzz, g_y_0_y_0_yyz_xxzz, g_y_0_y_0_yyz_xxzzz, g_y_0_y_0_yyz_xyyy, g_y_0_y_0_yyz_xyyyy, g_y_0_y_0_yyz_xyyyz, g_y_0_y_0_yyz_xyyz, g_y_0_y_0_yyz_xyyzz, g_y_0_y_0_yyz_xyzz, g_y_0_y_0_yyz_xyzzz, g_y_0_y_0_yyz_xzzz, g_y_0_y_0_yyz_xzzzz, g_y_0_y_0_yyz_yyyy, g_y_0_y_0_yyz_yyyz, g_y_0_y_0_yyz_yyzz, g_y_0_y_0_yyz_yzzz, g_y_0_y_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyz_xxxx[k] = -g_y_0_y_0_yyz_xxxx[k] * ab_x + g_y_0_y_0_yyz_xxxxx[k];

                g_y_0_y_0_xyyz_xxxy[k] = -g_y_0_y_0_yyz_xxxy[k] * ab_x + g_y_0_y_0_yyz_xxxxy[k];

                g_y_0_y_0_xyyz_xxxz[k] = -g_y_0_y_0_yyz_xxxz[k] * ab_x + g_y_0_y_0_yyz_xxxxz[k];

                g_y_0_y_0_xyyz_xxyy[k] = -g_y_0_y_0_yyz_xxyy[k] * ab_x + g_y_0_y_0_yyz_xxxyy[k];

                g_y_0_y_0_xyyz_xxyz[k] = -g_y_0_y_0_yyz_xxyz[k] * ab_x + g_y_0_y_0_yyz_xxxyz[k];

                g_y_0_y_0_xyyz_xxzz[k] = -g_y_0_y_0_yyz_xxzz[k] * ab_x + g_y_0_y_0_yyz_xxxzz[k];

                g_y_0_y_0_xyyz_xyyy[k] = -g_y_0_y_0_yyz_xyyy[k] * ab_x + g_y_0_y_0_yyz_xxyyy[k];

                g_y_0_y_0_xyyz_xyyz[k] = -g_y_0_y_0_yyz_xyyz[k] * ab_x + g_y_0_y_0_yyz_xxyyz[k];

                g_y_0_y_0_xyyz_xyzz[k] = -g_y_0_y_0_yyz_xyzz[k] * ab_x + g_y_0_y_0_yyz_xxyzz[k];

                g_y_0_y_0_xyyz_xzzz[k] = -g_y_0_y_0_yyz_xzzz[k] * ab_x + g_y_0_y_0_yyz_xxzzz[k];

                g_y_0_y_0_xyyz_yyyy[k] = -g_y_0_y_0_yyz_yyyy[k] * ab_x + g_y_0_y_0_yyz_xyyyy[k];

                g_y_0_y_0_xyyz_yyyz[k] = -g_y_0_y_0_yyz_yyyz[k] * ab_x + g_y_0_y_0_yyz_xyyyz[k];

                g_y_0_y_0_xyyz_yyzz[k] = -g_y_0_y_0_yyz_yyzz[k] * ab_x + g_y_0_y_0_yyz_xyyzz[k];

                g_y_0_y_0_xyyz_yzzz[k] = -g_y_0_y_0_yyz_yzzz[k] * ab_x + g_y_0_y_0_yyz_xyzzz[k];

                g_y_0_y_0_xyyz_zzzz[k] = -g_y_0_y_0_yyz_zzzz[k] * ab_x + g_y_0_y_0_yyz_xzzzz[k];
            }

            /// Set up 1020-1035 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1020 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1021 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1022 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1023 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1024 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1025 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1026 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1027 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1028 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1029 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1030 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1031 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1032 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1033 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1034 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyzz_xxxx, g_y_0_y_0_xyzz_xxxy, g_y_0_y_0_xyzz_xxxz, g_y_0_y_0_xyzz_xxyy, g_y_0_y_0_xyzz_xxyz, g_y_0_y_0_xyzz_xxzz, g_y_0_y_0_xyzz_xyyy, g_y_0_y_0_xyzz_xyyz, g_y_0_y_0_xyzz_xyzz, g_y_0_y_0_xyzz_xzzz, g_y_0_y_0_xyzz_yyyy, g_y_0_y_0_xyzz_yyyz, g_y_0_y_0_xyzz_yyzz, g_y_0_y_0_xyzz_yzzz, g_y_0_y_0_xyzz_zzzz, g_y_0_y_0_yzz_xxxx, g_y_0_y_0_yzz_xxxxx, g_y_0_y_0_yzz_xxxxy, g_y_0_y_0_yzz_xxxxz, g_y_0_y_0_yzz_xxxy, g_y_0_y_0_yzz_xxxyy, g_y_0_y_0_yzz_xxxyz, g_y_0_y_0_yzz_xxxz, g_y_0_y_0_yzz_xxxzz, g_y_0_y_0_yzz_xxyy, g_y_0_y_0_yzz_xxyyy, g_y_0_y_0_yzz_xxyyz, g_y_0_y_0_yzz_xxyz, g_y_0_y_0_yzz_xxyzz, g_y_0_y_0_yzz_xxzz, g_y_0_y_0_yzz_xxzzz, g_y_0_y_0_yzz_xyyy, g_y_0_y_0_yzz_xyyyy, g_y_0_y_0_yzz_xyyyz, g_y_0_y_0_yzz_xyyz, g_y_0_y_0_yzz_xyyzz, g_y_0_y_0_yzz_xyzz, g_y_0_y_0_yzz_xyzzz, g_y_0_y_0_yzz_xzzz, g_y_0_y_0_yzz_xzzzz, g_y_0_y_0_yzz_yyyy, g_y_0_y_0_yzz_yyyz, g_y_0_y_0_yzz_yyzz, g_y_0_y_0_yzz_yzzz, g_y_0_y_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyzz_xxxx[k] = -g_y_0_y_0_yzz_xxxx[k] * ab_x + g_y_0_y_0_yzz_xxxxx[k];

                g_y_0_y_0_xyzz_xxxy[k] = -g_y_0_y_0_yzz_xxxy[k] * ab_x + g_y_0_y_0_yzz_xxxxy[k];

                g_y_0_y_0_xyzz_xxxz[k] = -g_y_0_y_0_yzz_xxxz[k] * ab_x + g_y_0_y_0_yzz_xxxxz[k];

                g_y_0_y_0_xyzz_xxyy[k] = -g_y_0_y_0_yzz_xxyy[k] * ab_x + g_y_0_y_0_yzz_xxxyy[k];

                g_y_0_y_0_xyzz_xxyz[k] = -g_y_0_y_0_yzz_xxyz[k] * ab_x + g_y_0_y_0_yzz_xxxyz[k];

                g_y_0_y_0_xyzz_xxzz[k] = -g_y_0_y_0_yzz_xxzz[k] * ab_x + g_y_0_y_0_yzz_xxxzz[k];

                g_y_0_y_0_xyzz_xyyy[k] = -g_y_0_y_0_yzz_xyyy[k] * ab_x + g_y_0_y_0_yzz_xxyyy[k];

                g_y_0_y_0_xyzz_xyyz[k] = -g_y_0_y_0_yzz_xyyz[k] * ab_x + g_y_0_y_0_yzz_xxyyz[k];

                g_y_0_y_0_xyzz_xyzz[k] = -g_y_0_y_0_yzz_xyzz[k] * ab_x + g_y_0_y_0_yzz_xxyzz[k];

                g_y_0_y_0_xyzz_xzzz[k] = -g_y_0_y_0_yzz_xzzz[k] * ab_x + g_y_0_y_0_yzz_xxzzz[k];

                g_y_0_y_0_xyzz_yyyy[k] = -g_y_0_y_0_yzz_yyyy[k] * ab_x + g_y_0_y_0_yzz_xyyyy[k];

                g_y_0_y_0_xyzz_yyyz[k] = -g_y_0_y_0_yzz_yyyz[k] * ab_x + g_y_0_y_0_yzz_xyyyz[k];

                g_y_0_y_0_xyzz_yyzz[k] = -g_y_0_y_0_yzz_yyzz[k] * ab_x + g_y_0_y_0_yzz_xyyzz[k];

                g_y_0_y_0_xyzz_yzzz[k] = -g_y_0_y_0_yzz_yzzz[k] * ab_x + g_y_0_y_0_yzz_xyzzz[k];

                g_y_0_y_0_xyzz_zzzz[k] = -g_y_0_y_0_yzz_zzzz[k] * ab_x + g_y_0_y_0_yzz_xzzzz[k];
            }

            /// Set up 1035-1050 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1035 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1036 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1037 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1038 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1039 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1040 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1041 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1042 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1043 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1044 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1045 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1046 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1047 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1048 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xzzz_xxxx, g_y_0_y_0_xzzz_xxxy, g_y_0_y_0_xzzz_xxxz, g_y_0_y_0_xzzz_xxyy, g_y_0_y_0_xzzz_xxyz, g_y_0_y_0_xzzz_xxzz, g_y_0_y_0_xzzz_xyyy, g_y_0_y_0_xzzz_xyyz, g_y_0_y_0_xzzz_xyzz, g_y_0_y_0_xzzz_xzzz, g_y_0_y_0_xzzz_yyyy, g_y_0_y_0_xzzz_yyyz, g_y_0_y_0_xzzz_yyzz, g_y_0_y_0_xzzz_yzzz, g_y_0_y_0_xzzz_zzzz, g_y_0_y_0_zzz_xxxx, g_y_0_y_0_zzz_xxxxx, g_y_0_y_0_zzz_xxxxy, g_y_0_y_0_zzz_xxxxz, g_y_0_y_0_zzz_xxxy, g_y_0_y_0_zzz_xxxyy, g_y_0_y_0_zzz_xxxyz, g_y_0_y_0_zzz_xxxz, g_y_0_y_0_zzz_xxxzz, g_y_0_y_0_zzz_xxyy, g_y_0_y_0_zzz_xxyyy, g_y_0_y_0_zzz_xxyyz, g_y_0_y_0_zzz_xxyz, g_y_0_y_0_zzz_xxyzz, g_y_0_y_0_zzz_xxzz, g_y_0_y_0_zzz_xxzzz, g_y_0_y_0_zzz_xyyy, g_y_0_y_0_zzz_xyyyy, g_y_0_y_0_zzz_xyyyz, g_y_0_y_0_zzz_xyyz, g_y_0_y_0_zzz_xyyzz, g_y_0_y_0_zzz_xyzz, g_y_0_y_0_zzz_xyzzz, g_y_0_y_0_zzz_xzzz, g_y_0_y_0_zzz_xzzzz, g_y_0_y_0_zzz_yyyy, g_y_0_y_0_zzz_yyyz, g_y_0_y_0_zzz_yyzz, g_y_0_y_0_zzz_yzzz, g_y_0_y_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xzzz_xxxx[k] = -g_y_0_y_0_zzz_xxxx[k] * ab_x + g_y_0_y_0_zzz_xxxxx[k];

                g_y_0_y_0_xzzz_xxxy[k] = -g_y_0_y_0_zzz_xxxy[k] * ab_x + g_y_0_y_0_zzz_xxxxy[k];

                g_y_0_y_0_xzzz_xxxz[k] = -g_y_0_y_0_zzz_xxxz[k] * ab_x + g_y_0_y_0_zzz_xxxxz[k];

                g_y_0_y_0_xzzz_xxyy[k] = -g_y_0_y_0_zzz_xxyy[k] * ab_x + g_y_0_y_0_zzz_xxxyy[k];

                g_y_0_y_0_xzzz_xxyz[k] = -g_y_0_y_0_zzz_xxyz[k] * ab_x + g_y_0_y_0_zzz_xxxyz[k];

                g_y_0_y_0_xzzz_xxzz[k] = -g_y_0_y_0_zzz_xxzz[k] * ab_x + g_y_0_y_0_zzz_xxxzz[k];

                g_y_0_y_0_xzzz_xyyy[k] = -g_y_0_y_0_zzz_xyyy[k] * ab_x + g_y_0_y_0_zzz_xxyyy[k];

                g_y_0_y_0_xzzz_xyyz[k] = -g_y_0_y_0_zzz_xyyz[k] * ab_x + g_y_0_y_0_zzz_xxyyz[k];

                g_y_0_y_0_xzzz_xyzz[k] = -g_y_0_y_0_zzz_xyzz[k] * ab_x + g_y_0_y_0_zzz_xxyzz[k];

                g_y_0_y_0_xzzz_xzzz[k] = -g_y_0_y_0_zzz_xzzz[k] * ab_x + g_y_0_y_0_zzz_xxzzz[k];

                g_y_0_y_0_xzzz_yyyy[k] = -g_y_0_y_0_zzz_yyyy[k] * ab_x + g_y_0_y_0_zzz_xyyyy[k];

                g_y_0_y_0_xzzz_yyyz[k] = -g_y_0_y_0_zzz_yyyz[k] * ab_x + g_y_0_y_0_zzz_xyyyz[k];

                g_y_0_y_0_xzzz_yyzz[k] = -g_y_0_y_0_zzz_yyzz[k] * ab_x + g_y_0_y_0_zzz_xyyzz[k];

                g_y_0_y_0_xzzz_yzzz[k] = -g_y_0_y_0_zzz_yzzz[k] * ab_x + g_y_0_y_0_zzz_xyzzz[k];

                g_y_0_y_0_xzzz_zzzz[k] = -g_y_0_y_0_zzz_zzzz[k] * ab_x + g_y_0_y_0_zzz_xzzzz[k];
            }

            /// Set up 1050-1065 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyy_xxxx = cbuffer.data(gg_geom_1010_off + 1050 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_xxxy = cbuffer.data(gg_geom_1010_off + 1051 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_xxxz = cbuffer.data(gg_geom_1010_off + 1052 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_xxyy = cbuffer.data(gg_geom_1010_off + 1053 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_xxyz = cbuffer.data(gg_geom_1010_off + 1054 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_xxzz = cbuffer.data(gg_geom_1010_off + 1055 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_xyyy = cbuffer.data(gg_geom_1010_off + 1056 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_xyyz = cbuffer.data(gg_geom_1010_off + 1057 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_xyzz = cbuffer.data(gg_geom_1010_off + 1058 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_xzzz = cbuffer.data(gg_geom_1010_off + 1059 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1060 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1061 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1062 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1063 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1064 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_yyy_xxxx, g_0_0_y_0_yyy_xxxy, g_0_0_y_0_yyy_xxxz, g_0_0_y_0_yyy_xxyy, g_0_0_y_0_yyy_xxyz, g_0_0_y_0_yyy_xxzz, g_0_0_y_0_yyy_xyyy, g_0_0_y_0_yyy_xyyz, g_0_0_y_0_yyy_xyzz, g_0_0_y_0_yyy_xzzz, g_0_0_y_0_yyy_yyyy, g_0_0_y_0_yyy_yyyz, g_0_0_y_0_yyy_yyzz, g_0_0_y_0_yyy_yzzz, g_0_0_y_0_yyy_zzzz, g_y_0_y_0_yyy_xxxx, g_y_0_y_0_yyy_xxxxy, g_y_0_y_0_yyy_xxxy, g_y_0_y_0_yyy_xxxyy, g_y_0_y_0_yyy_xxxyz, g_y_0_y_0_yyy_xxxz, g_y_0_y_0_yyy_xxyy, g_y_0_y_0_yyy_xxyyy, g_y_0_y_0_yyy_xxyyz, g_y_0_y_0_yyy_xxyz, g_y_0_y_0_yyy_xxyzz, g_y_0_y_0_yyy_xxzz, g_y_0_y_0_yyy_xyyy, g_y_0_y_0_yyy_xyyyy, g_y_0_y_0_yyy_xyyyz, g_y_0_y_0_yyy_xyyz, g_y_0_y_0_yyy_xyyzz, g_y_0_y_0_yyy_xyzz, g_y_0_y_0_yyy_xyzzz, g_y_0_y_0_yyy_xzzz, g_y_0_y_0_yyy_yyyy, g_y_0_y_0_yyy_yyyyy, g_y_0_y_0_yyy_yyyyz, g_y_0_y_0_yyy_yyyz, g_y_0_y_0_yyy_yyyzz, g_y_0_y_0_yyy_yyzz, g_y_0_y_0_yyy_yyzzz, g_y_0_y_0_yyy_yzzz, g_y_0_y_0_yyy_yzzzz, g_y_0_y_0_yyy_zzzz, g_y_0_y_0_yyyy_xxxx, g_y_0_y_0_yyyy_xxxy, g_y_0_y_0_yyyy_xxxz, g_y_0_y_0_yyyy_xxyy, g_y_0_y_0_yyyy_xxyz, g_y_0_y_0_yyyy_xxzz, g_y_0_y_0_yyyy_xyyy, g_y_0_y_0_yyyy_xyyz, g_y_0_y_0_yyyy_xyzz, g_y_0_y_0_yyyy_xzzz, g_y_0_y_0_yyyy_yyyy, g_y_0_y_0_yyyy_yyyz, g_y_0_y_0_yyyy_yyzz, g_y_0_y_0_yyyy_yzzz, g_y_0_y_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyy_xxxx[k] = -g_0_0_y_0_yyy_xxxx[k] - g_y_0_y_0_yyy_xxxx[k] * ab_y + g_y_0_y_0_yyy_xxxxy[k];

                g_y_0_y_0_yyyy_xxxy[k] = -g_0_0_y_0_yyy_xxxy[k] - g_y_0_y_0_yyy_xxxy[k] * ab_y + g_y_0_y_0_yyy_xxxyy[k];

                g_y_0_y_0_yyyy_xxxz[k] = -g_0_0_y_0_yyy_xxxz[k] - g_y_0_y_0_yyy_xxxz[k] * ab_y + g_y_0_y_0_yyy_xxxyz[k];

                g_y_0_y_0_yyyy_xxyy[k] = -g_0_0_y_0_yyy_xxyy[k] - g_y_0_y_0_yyy_xxyy[k] * ab_y + g_y_0_y_0_yyy_xxyyy[k];

                g_y_0_y_0_yyyy_xxyz[k] = -g_0_0_y_0_yyy_xxyz[k] - g_y_0_y_0_yyy_xxyz[k] * ab_y + g_y_0_y_0_yyy_xxyyz[k];

                g_y_0_y_0_yyyy_xxzz[k] = -g_0_0_y_0_yyy_xxzz[k] - g_y_0_y_0_yyy_xxzz[k] * ab_y + g_y_0_y_0_yyy_xxyzz[k];

                g_y_0_y_0_yyyy_xyyy[k] = -g_0_0_y_0_yyy_xyyy[k] - g_y_0_y_0_yyy_xyyy[k] * ab_y + g_y_0_y_0_yyy_xyyyy[k];

                g_y_0_y_0_yyyy_xyyz[k] = -g_0_0_y_0_yyy_xyyz[k] - g_y_0_y_0_yyy_xyyz[k] * ab_y + g_y_0_y_0_yyy_xyyyz[k];

                g_y_0_y_0_yyyy_xyzz[k] = -g_0_0_y_0_yyy_xyzz[k] - g_y_0_y_0_yyy_xyzz[k] * ab_y + g_y_0_y_0_yyy_xyyzz[k];

                g_y_0_y_0_yyyy_xzzz[k] = -g_0_0_y_0_yyy_xzzz[k] - g_y_0_y_0_yyy_xzzz[k] * ab_y + g_y_0_y_0_yyy_xyzzz[k];

                g_y_0_y_0_yyyy_yyyy[k] = -g_0_0_y_0_yyy_yyyy[k] - g_y_0_y_0_yyy_yyyy[k] * ab_y + g_y_0_y_0_yyy_yyyyy[k];

                g_y_0_y_0_yyyy_yyyz[k] = -g_0_0_y_0_yyy_yyyz[k] - g_y_0_y_0_yyy_yyyz[k] * ab_y + g_y_0_y_0_yyy_yyyyz[k];

                g_y_0_y_0_yyyy_yyzz[k] = -g_0_0_y_0_yyy_yyzz[k] - g_y_0_y_0_yyy_yyzz[k] * ab_y + g_y_0_y_0_yyy_yyyzz[k];

                g_y_0_y_0_yyyy_yzzz[k] = -g_0_0_y_0_yyy_yzzz[k] - g_y_0_y_0_yyy_yzzz[k] * ab_y + g_y_0_y_0_yyy_yyzzz[k];

                g_y_0_y_0_yyyy_zzzz[k] = -g_0_0_y_0_yyy_zzzz[k] - g_y_0_y_0_yyy_zzzz[k] * ab_y + g_y_0_y_0_yyy_yzzzz[k];
            }

            /// Set up 1065-1080 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1065 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1066 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1067 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1068 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1069 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1070 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1071 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1072 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1073 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1074 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1075 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1076 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1077 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1078 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyy_xxxx, g_y_0_y_0_yyy_xxxxz, g_y_0_y_0_yyy_xxxy, g_y_0_y_0_yyy_xxxyz, g_y_0_y_0_yyy_xxxz, g_y_0_y_0_yyy_xxxzz, g_y_0_y_0_yyy_xxyy, g_y_0_y_0_yyy_xxyyz, g_y_0_y_0_yyy_xxyz, g_y_0_y_0_yyy_xxyzz, g_y_0_y_0_yyy_xxzz, g_y_0_y_0_yyy_xxzzz, g_y_0_y_0_yyy_xyyy, g_y_0_y_0_yyy_xyyyz, g_y_0_y_0_yyy_xyyz, g_y_0_y_0_yyy_xyyzz, g_y_0_y_0_yyy_xyzz, g_y_0_y_0_yyy_xyzzz, g_y_0_y_0_yyy_xzzz, g_y_0_y_0_yyy_xzzzz, g_y_0_y_0_yyy_yyyy, g_y_0_y_0_yyy_yyyyz, g_y_0_y_0_yyy_yyyz, g_y_0_y_0_yyy_yyyzz, g_y_0_y_0_yyy_yyzz, g_y_0_y_0_yyy_yyzzz, g_y_0_y_0_yyy_yzzz, g_y_0_y_0_yyy_yzzzz, g_y_0_y_0_yyy_zzzz, g_y_0_y_0_yyy_zzzzz, g_y_0_y_0_yyyz_xxxx, g_y_0_y_0_yyyz_xxxy, g_y_0_y_0_yyyz_xxxz, g_y_0_y_0_yyyz_xxyy, g_y_0_y_0_yyyz_xxyz, g_y_0_y_0_yyyz_xxzz, g_y_0_y_0_yyyz_xyyy, g_y_0_y_0_yyyz_xyyz, g_y_0_y_0_yyyz_xyzz, g_y_0_y_0_yyyz_xzzz, g_y_0_y_0_yyyz_yyyy, g_y_0_y_0_yyyz_yyyz, g_y_0_y_0_yyyz_yyzz, g_y_0_y_0_yyyz_yzzz, g_y_0_y_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyz_xxxx[k] = -g_y_0_y_0_yyy_xxxx[k] * ab_z + g_y_0_y_0_yyy_xxxxz[k];

                g_y_0_y_0_yyyz_xxxy[k] = -g_y_0_y_0_yyy_xxxy[k] * ab_z + g_y_0_y_0_yyy_xxxyz[k];

                g_y_0_y_0_yyyz_xxxz[k] = -g_y_0_y_0_yyy_xxxz[k] * ab_z + g_y_0_y_0_yyy_xxxzz[k];

                g_y_0_y_0_yyyz_xxyy[k] = -g_y_0_y_0_yyy_xxyy[k] * ab_z + g_y_0_y_0_yyy_xxyyz[k];

                g_y_0_y_0_yyyz_xxyz[k] = -g_y_0_y_0_yyy_xxyz[k] * ab_z + g_y_0_y_0_yyy_xxyzz[k];

                g_y_0_y_0_yyyz_xxzz[k] = -g_y_0_y_0_yyy_xxzz[k] * ab_z + g_y_0_y_0_yyy_xxzzz[k];

                g_y_0_y_0_yyyz_xyyy[k] = -g_y_0_y_0_yyy_xyyy[k] * ab_z + g_y_0_y_0_yyy_xyyyz[k];

                g_y_0_y_0_yyyz_xyyz[k] = -g_y_0_y_0_yyy_xyyz[k] * ab_z + g_y_0_y_0_yyy_xyyzz[k];

                g_y_0_y_0_yyyz_xyzz[k] = -g_y_0_y_0_yyy_xyzz[k] * ab_z + g_y_0_y_0_yyy_xyzzz[k];

                g_y_0_y_0_yyyz_xzzz[k] = -g_y_0_y_0_yyy_xzzz[k] * ab_z + g_y_0_y_0_yyy_xzzzz[k];

                g_y_0_y_0_yyyz_yyyy[k] = -g_y_0_y_0_yyy_yyyy[k] * ab_z + g_y_0_y_0_yyy_yyyyz[k];

                g_y_0_y_0_yyyz_yyyz[k] = -g_y_0_y_0_yyy_yyyz[k] * ab_z + g_y_0_y_0_yyy_yyyzz[k];

                g_y_0_y_0_yyyz_yyzz[k] = -g_y_0_y_0_yyy_yyzz[k] * ab_z + g_y_0_y_0_yyy_yyzzz[k];

                g_y_0_y_0_yyyz_yzzz[k] = -g_y_0_y_0_yyy_yzzz[k] * ab_z + g_y_0_y_0_yyy_yzzzz[k];

                g_y_0_y_0_yyyz_zzzz[k] = -g_y_0_y_0_yyy_zzzz[k] * ab_z + g_y_0_y_0_yyy_zzzzz[k];
            }

            /// Set up 1080-1095 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1080 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1081 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1082 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1083 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1084 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1085 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1086 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1087 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1088 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1089 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1090 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1091 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1092 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1093 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1094 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyz_xxxx, g_y_0_y_0_yyz_xxxxz, g_y_0_y_0_yyz_xxxy, g_y_0_y_0_yyz_xxxyz, g_y_0_y_0_yyz_xxxz, g_y_0_y_0_yyz_xxxzz, g_y_0_y_0_yyz_xxyy, g_y_0_y_0_yyz_xxyyz, g_y_0_y_0_yyz_xxyz, g_y_0_y_0_yyz_xxyzz, g_y_0_y_0_yyz_xxzz, g_y_0_y_0_yyz_xxzzz, g_y_0_y_0_yyz_xyyy, g_y_0_y_0_yyz_xyyyz, g_y_0_y_0_yyz_xyyz, g_y_0_y_0_yyz_xyyzz, g_y_0_y_0_yyz_xyzz, g_y_0_y_0_yyz_xyzzz, g_y_0_y_0_yyz_xzzz, g_y_0_y_0_yyz_xzzzz, g_y_0_y_0_yyz_yyyy, g_y_0_y_0_yyz_yyyyz, g_y_0_y_0_yyz_yyyz, g_y_0_y_0_yyz_yyyzz, g_y_0_y_0_yyz_yyzz, g_y_0_y_0_yyz_yyzzz, g_y_0_y_0_yyz_yzzz, g_y_0_y_0_yyz_yzzzz, g_y_0_y_0_yyz_zzzz, g_y_0_y_0_yyz_zzzzz, g_y_0_y_0_yyzz_xxxx, g_y_0_y_0_yyzz_xxxy, g_y_0_y_0_yyzz_xxxz, g_y_0_y_0_yyzz_xxyy, g_y_0_y_0_yyzz_xxyz, g_y_0_y_0_yyzz_xxzz, g_y_0_y_0_yyzz_xyyy, g_y_0_y_0_yyzz_xyyz, g_y_0_y_0_yyzz_xyzz, g_y_0_y_0_yyzz_xzzz, g_y_0_y_0_yyzz_yyyy, g_y_0_y_0_yyzz_yyyz, g_y_0_y_0_yyzz_yyzz, g_y_0_y_0_yyzz_yzzz, g_y_0_y_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyzz_xxxx[k] = -g_y_0_y_0_yyz_xxxx[k] * ab_z + g_y_0_y_0_yyz_xxxxz[k];

                g_y_0_y_0_yyzz_xxxy[k] = -g_y_0_y_0_yyz_xxxy[k] * ab_z + g_y_0_y_0_yyz_xxxyz[k];

                g_y_0_y_0_yyzz_xxxz[k] = -g_y_0_y_0_yyz_xxxz[k] * ab_z + g_y_0_y_0_yyz_xxxzz[k];

                g_y_0_y_0_yyzz_xxyy[k] = -g_y_0_y_0_yyz_xxyy[k] * ab_z + g_y_0_y_0_yyz_xxyyz[k];

                g_y_0_y_0_yyzz_xxyz[k] = -g_y_0_y_0_yyz_xxyz[k] * ab_z + g_y_0_y_0_yyz_xxyzz[k];

                g_y_0_y_0_yyzz_xxzz[k] = -g_y_0_y_0_yyz_xxzz[k] * ab_z + g_y_0_y_0_yyz_xxzzz[k];

                g_y_0_y_0_yyzz_xyyy[k] = -g_y_0_y_0_yyz_xyyy[k] * ab_z + g_y_0_y_0_yyz_xyyyz[k];

                g_y_0_y_0_yyzz_xyyz[k] = -g_y_0_y_0_yyz_xyyz[k] * ab_z + g_y_0_y_0_yyz_xyyzz[k];

                g_y_0_y_0_yyzz_xyzz[k] = -g_y_0_y_0_yyz_xyzz[k] * ab_z + g_y_0_y_0_yyz_xyzzz[k];

                g_y_0_y_0_yyzz_xzzz[k] = -g_y_0_y_0_yyz_xzzz[k] * ab_z + g_y_0_y_0_yyz_xzzzz[k];

                g_y_0_y_0_yyzz_yyyy[k] = -g_y_0_y_0_yyz_yyyy[k] * ab_z + g_y_0_y_0_yyz_yyyyz[k];

                g_y_0_y_0_yyzz_yyyz[k] = -g_y_0_y_0_yyz_yyyz[k] * ab_z + g_y_0_y_0_yyz_yyyzz[k];

                g_y_0_y_0_yyzz_yyzz[k] = -g_y_0_y_0_yyz_yyzz[k] * ab_z + g_y_0_y_0_yyz_yyzzz[k];

                g_y_0_y_0_yyzz_yzzz[k] = -g_y_0_y_0_yyz_yzzz[k] * ab_z + g_y_0_y_0_yyz_yzzzz[k];

                g_y_0_y_0_yyzz_zzzz[k] = -g_y_0_y_0_yyz_zzzz[k] * ab_z + g_y_0_y_0_yyz_zzzzz[k];
            }

            /// Set up 1095-1110 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1095 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1096 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1097 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1098 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1099 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1100 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1101 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1102 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1103 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1104 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1105 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1106 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1107 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1108 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1109 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yzz_xxxx, g_y_0_y_0_yzz_xxxxz, g_y_0_y_0_yzz_xxxy, g_y_0_y_0_yzz_xxxyz, g_y_0_y_0_yzz_xxxz, g_y_0_y_0_yzz_xxxzz, g_y_0_y_0_yzz_xxyy, g_y_0_y_0_yzz_xxyyz, g_y_0_y_0_yzz_xxyz, g_y_0_y_0_yzz_xxyzz, g_y_0_y_0_yzz_xxzz, g_y_0_y_0_yzz_xxzzz, g_y_0_y_0_yzz_xyyy, g_y_0_y_0_yzz_xyyyz, g_y_0_y_0_yzz_xyyz, g_y_0_y_0_yzz_xyyzz, g_y_0_y_0_yzz_xyzz, g_y_0_y_0_yzz_xyzzz, g_y_0_y_0_yzz_xzzz, g_y_0_y_0_yzz_xzzzz, g_y_0_y_0_yzz_yyyy, g_y_0_y_0_yzz_yyyyz, g_y_0_y_0_yzz_yyyz, g_y_0_y_0_yzz_yyyzz, g_y_0_y_0_yzz_yyzz, g_y_0_y_0_yzz_yyzzz, g_y_0_y_0_yzz_yzzz, g_y_0_y_0_yzz_yzzzz, g_y_0_y_0_yzz_zzzz, g_y_0_y_0_yzz_zzzzz, g_y_0_y_0_yzzz_xxxx, g_y_0_y_0_yzzz_xxxy, g_y_0_y_0_yzzz_xxxz, g_y_0_y_0_yzzz_xxyy, g_y_0_y_0_yzzz_xxyz, g_y_0_y_0_yzzz_xxzz, g_y_0_y_0_yzzz_xyyy, g_y_0_y_0_yzzz_xyyz, g_y_0_y_0_yzzz_xyzz, g_y_0_y_0_yzzz_xzzz, g_y_0_y_0_yzzz_yyyy, g_y_0_y_0_yzzz_yyyz, g_y_0_y_0_yzzz_yyzz, g_y_0_y_0_yzzz_yzzz, g_y_0_y_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yzzz_xxxx[k] = -g_y_0_y_0_yzz_xxxx[k] * ab_z + g_y_0_y_0_yzz_xxxxz[k];

                g_y_0_y_0_yzzz_xxxy[k] = -g_y_0_y_0_yzz_xxxy[k] * ab_z + g_y_0_y_0_yzz_xxxyz[k];

                g_y_0_y_0_yzzz_xxxz[k] = -g_y_0_y_0_yzz_xxxz[k] * ab_z + g_y_0_y_0_yzz_xxxzz[k];

                g_y_0_y_0_yzzz_xxyy[k] = -g_y_0_y_0_yzz_xxyy[k] * ab_z + g_y_0_y_0_yzz_xxyyz[k];

                g_y_0_y_0_yzzz_xxyz[k] = -g_y_0_y_0_yzz_xxyz[k] * ab_z + g_y_0_y_0_yzz_xxyzz[k];

                g_y_0_y_0_yzzz_xxzz[k] = -g_y_0_y_0_yzz_xxzz[k] * ab_z + g_y_0_y_0_yzz_xxzzz[k];

                g_y_0_y_0_yzzz_xyyy[k] = -g_y_0_y_0_yzz_xyyy[k] * ab_z + g_y_0_y_0_yzz_xyyyz[k];

                g_y_0_y_0_yzzz_xyyz[k] = -g_y_0_y_0_yzz_xyyz[k] * ab_z + g_y_0_y_0_yzz_xyyzz[k];

                g_y_0_y_0_yzzz_xyzz[k] = -g_y_0_y_0_yzz_xyzz[k] * ab_z + g_y_0_y_0_yzz_xyzzz[k];

                g_y_0_y_0_yzzz_xzzz[k] = -g_y_0_y_0_yzz_xzzz[k] * ab_z + g_y_0_y_0_yzz_xzzzz[k];

                g_y_0_y_0_yzzz_yyyy[k] = -g_y_0_y_0_yzz_yyyy[k] * ab_z + g_y_0_y_0_yzz_yyyyz[k];

                g_y_0_y_0_yzzz_yyyz[k] = -g_y_0_y_0_yzz_yyyz[k] * ab_z + g_y_0_y_0_yzz_yyyzz[k];

                g_y_0_y_0_yzzz_yyzz[k] = -g_y_0_y_0_yzz_yyzz[k] * ab_z + g_y_0_y_0_yzz_yyzzz[k];

                g_y_0_y_0_yzzz_yzzz[k] = -g_y_0_y_0_yzz_yzzz[k] * ab_z + g_y_0_y_0_yzz_yzzzz[k];

                g_y_0_y_0_yzzz_zzzz[k] = -g_y_0_y_0_yzz_zzzz[k] * ab_z + g_y_0_y_0_yzz_zzzzz[k];
            }

            /// Set up 1110-1125 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_zzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1110 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1111 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1112 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1113 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1114 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1115 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1116 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1117 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1118 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1119 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1120 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1121 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1122 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1123 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1124 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_zzz_xxxx, g_y_0_y_0_zzz_xxxxz, g_y_0_y_0_zzz_xxxy, g_y_0_y_0_zzz_xxxyz, g_y_0_y_0_zzz_xxxz, g_y_0_y_0_zzz_xxxzz, g_y_0_y_0_zzz_xxyy, g_y_0_y_0_zzz_xxyyz, g_y_0_y_0_zzz_xxyz, g_y_0_y_0_zzz_xxyzz, g_y_0_y_0_zzz_xxzz, g_y_0_y_0_zzz_xxzzz, g_y_0_y_0_zzz_xyyy, g_y_0_y_0_zzz_xyyyz, g_y_0_y_0_zzz_xyyz, g_y_0_y_0_zzz_xyyzz, g_y_0_y_0_zzz_xyzz, g_y_0_y_0_zzz_xyzzz, g_y_0_y_0_zzz_xzzz, g_y_0_y_0_zzz_xzzzz, g_y_0_y_0_zzz_yyyy, g_y_0_y_0_zzz_yyyyz, g_y_0_y_0_zzz_yyyz, g_y_0_y_0_zzz_yyyzz, g_y_0_y_0_zzz_yyzz, g_y_0_y_0_zzz_yyzzz, g_y_0_y_0_zzz_yzzz, g_y_0_y_0_zzz_yzzzz, g_y_0_y_0_zzz_zzzz, g_y_0_y_0_zzz_zzzzz, g_y_0_y_0_zzzz_xxxx, g_y_0_y_0_zzzz_xxxy, g_y_0_y_0_zzzz_xxxz, g_y_0_y_0_zzzz_xxyy, g_y_0_y_0_zzzz_xxyz, g_y_0_y_0_zzzz_xxzz, g_y_0_y_0_zzzz_xyyy, g_y_0_y_0_zzzz_xyyz, g_y_0_y_0_zzzz_xyzz, g_y_0_y_0_zzzz_xzzz, g_y_0_y_0_zzzz_yyyy, g_y_0_y_0_zzzz_yyyz, g_y_0_y_0_zzzz_yyzz, g_y_0_y_0_zzzz_yzzz, g_y_0_y_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_zzzz_xxxx[k] = -g_y_0_y_0_zzz_xxxx[k] * ab_z + g_y_0_y_0_zzz_xxxxz[k];

                g_y_0_y_0_zzzz_xxxy[k] = -g_y_0_y_0_zzz_xxxy[k] * ab_z + g_y_0_y_0_zzz_xxxyz[k];

                g_y_0_y_0_zzzz_xxxz[k] = -g_y_0_y_0_zzz_xxxz[k] * ab_z + g_y_0_y_0_zzz_xxxzz[k];

                g_y_0_y_0_zzzz_xxyy[k] = -g_y_0_y_0_zzz_xxyy[k] * ab_z + g_y_0_y_0_zzz_xxyyz[k];

                g_y_0_y_0_zzzz_xxyz[k] = -g_y_0_y_0_zzz_xxyz[k] * ab_z + g_y_0_y_0_zzz_xxyzz[k];

                g_y_0_y_0_zzzz_xxzz[k] = -g_y_0_y_0_zzz_xxzz[k] * ab_z + g_y_0_y_0_zzz_xxzzz[k];

                g_y_0_y_0_zzzz_xyyy[k] = -g_y_0_y_0_zzz_xyyy[k] * ab_z + g_y_0_y_0_zzz_xyyyz[k];

                g_y_0_y_0_zzzz_xyyz[k] = -g_y_0_y_0_zzz_xyyz[k] * ab_z + g_y_0_y_0_zzz_xyyzz[k];

                g_y_0_y_0_zzzz_xyzz[k] = -g_y_0_y_0_zzz_xyzz[k] * ab_z + g_y_0_y_0_zzz_xyzzz[k];

                g_y_0_y_0_zzzz_xzzz[k] = -g_y_0_y_0_zzz_xzzz[k] * ab_z + g_y_0_y_0_zzz_xzzzz[k];

                g_y_0_y_0_zzzz_yyyy[k] = -g_y_0_y_0_zzz_yyyy[k] * ab_z + g_y_0_y_0_zzz_yyyyz[k];

                g_y_0_y_0_zzzz_yyyz[k] = -g_y_0_y_0_zzz_yyyz[k] * ab_z + g_y_0_y_0_zzz_yyyzz[k];

                g_y_0_y_0_zzzz_yyzz[k] = -g_y_0_y_0_zzz_yyzz[k] * ab_z + g_y_0_y_0_zzz_yyzzz[k];

                g_y_0_y_0_zzzz_yzzz[k] = -g_y_0_y_0_zzz_yzzz[k] * ab_z + g_y_0_y_0_zzz_yzzzz[k];

                g_y_0_y_0_zzzz_zzzz[k] = -g_y_0_y_0_zzz_zzzz[k] * ab_z + g_y_0_y_0_zzz_zzzzz[k];
            }

            /// Set up 1125-1140 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxx_xxxx = cbuffer.data(gg_geom_1010_off + 1125 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_xxxy = cbuffer.data(gg_geom_1010_off + 1126 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_xxxz = cbuffer.data(gg_geom_1010_off + 1127 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_xxyy = cbuffer.data(gg_geom_1010_off + 1128 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_xxyz = cbuffer.data(gg_geom_1010_off + 1129 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_xxzz = cbuffer.data(gg_geom_1010_off + 1130 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_xyyy = cbuffer.data(gg_geom_1010_off + 1131 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_xyyz = cbuffer.data(gg_geom_1010_off + 1132 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_xyzz = cbuffer.data(gg_geom_1010_off + 1133 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_xzzz = cbuffer.data(gg_geom_1010_off + 1134 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_yyyy = cbuffer.data(gg_geom_1010_off + 1135 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_yyyz = cbuffer.data(gg_geom_1010_off + 1136 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_yyzz = cbuffer.data(gg_geom_1010_off + 1137 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_yzzz = cbuffer.data(gg_geom_1010_off + 1138 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_zzzz = cbuffer.data(gg_geom_1010_off + 1139 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxx_xxxx, g_y_0_z_0_xxx_xxxxx, g_y_0_z_0_xxx_xxxxy, g_y_0_z_0_xxx_xxxxz, g_y_0_z_0_xxx_xxxy, g_y_0_z_0_xxx_xxxyy, g_y_0_z_0_xxx_xxxyz, g_y_0_z_0_xxx_xxxz, g_y_0_z_0_xxx_xxxzz, g_y_0_z_0_xxx_xxyy, g_y_0_z_0_xxx_xxyyy, g_y_0_z_0_xxx_xxyyz, g_y_0_z_0_xxx_xxyz, g_y_0_z_0_xxx_xxyzz, g_y_0_z_0_xxx_xxzz, g_y_0_z_0_xxx_xxzzz, g_y_0_z_0_xxx_xyyy, g_y_0_z_0_xxx_xyyyy, g_y_0_z_0_xxx_xyyyz, g_y_0_z_0_xxx_xyyz, g_y_0_z_0_xxx_xyyzz, g_y_0_z_0_xxx_xyzz, g_y_0_z_0_xxx_xyzzz, g_y_0_z_0_xxx_xzzz, g_y_0_z_0_xxx_xzzzz, g_y_0_z_0_xxx_yyyy, g_y_0_z_0_xxx_yyyz, g_y_0_z_0_xxx_yyzz, g_y_0_z_0_xxx_yzzz, g_y_0_z_0_xxx_zzzz, g_y_0_z_0_xxxx_xxxx, g_y_0_z_0_xxxx_xxxy, g_y_0_z_0_xxxx_xxxz, g_y_0_z_0_xxxx_xxyy, g_y_0_z_0_xxxx_xxyz, g_y_0_z_0_xxxx_xxzz, g_y_0_z_0_xxxx_xyyy, g_y_0_z_0_xxxx_xyyz, g_y_0_z_0_xxxx_xyzz, g_y_0_z_0_xxxx_xzzz, g_y_0_z_0_xxxx_yyyy, g_y_0_z_0_xxxx_yyyz, g_y_0_z_0_xxxx_yyzz, g_y_0_z_0_xxxx_yzzz, g_y_0_z_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxx_xxxx[k] = -g_y_0_z_0_xxx_xxxx[k] * ab_x + g_y_0_z_0_xxx_xxxxx[k];

                g_y_0_z_0_xxxx_xxxy[k] = -g_y_0_z_0_xxx_xxxy[k] * ab_x + g_y_0_z_0_xxx_xxxxy[k];

                g_y_0_z_0_xxxx_xxxz[k] = -g_y_0_z_0_xxx_xxxz[k] * ab_x + g_y_0_z_0_xxx_xxxxz[k];

                g_y_0_z_0_xxxx_xxyy[k] = -g_y_0_z_0_xxx_xxyy[k] * ab_x + g_y_0_z_0_xxx_xxxyy[k];

                g_y_0_z_0_xxxx_xxyz[k] = -g_y_0_z_0_xxx_xxyz[k] * ab_x + g_y_0_z_0_xxx_xxxyz[k];

                g_y_0_z_0_xxxx_xxzz[k] = -g_y_0_z_0_xxx_xxzz[k] * ab_x + g_y_0_z_0_xxx_xxxzz[k];

                g_y_0_z_0_xxxx_xyyy[k] = -g_y_0_z_0_xxx_xyyy[k] * ab_x + g_y_0_z_0_xxx_xxyyy[k];

                g_y_0_z_0_xxxx_xyyz[k] = -g_y_0_z_0_xxx_xyyz[k] * ab_x + g_y_0_z_0_xxx_xxyyz[k];

                g_y_0_z_0_xxxx_xyzz[k] = -g_y_0_z_0_xxx_xyzz[k] * ab_x + g_y_0_z_0_xxx_xxyzz[k];

                g_y_0_z_0_xxxx_xzzz[k] = -g_y_0_z_0_xxx_xzzz[k] * ab_x + g_y_0_z_0_xxx_xxzzz[k];

                g_y_0_z_0_xxxx_yyyy[k] = -g_y_0_z_0_xxx_yyyy[k] * ab_x + g_y_0_z_0_xxx_xyyyy[k];

                g_y_0_z_0_xxxx_yyyz[k] = -g_y_0_z_0_xxx_yyyz[k] * ab_x + g_y_0_z_0_xxx_xyyyz[k];

                g_y_0_z_0_xxxx_yyzz[k] = -g_y_0_z_0_xxx_yyzz[k] * ab_x + g_y_0_z_0_xxx_xyyzz[k];

                g_y_0_z_0_xxxx_yzzz[k] = -g_y_0_z_0_xxx_yzzz[k] * ab_x + g_y_0_z_0_xxx_xyzzz[k];

                g_y_0_z_0_xxxx_zzzz[k] = -g_y_0_z_0_xxx_zzzz[k] * ab_x + g_y_0_z_0_xxx_xzzzz[k];
            }

            /// Set up 1140-1155 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxy_xxxx = cbuffer.data(gg_geom_1010_off + 1140 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_xxxy = cbuffer.data(gg_geom_1010_off + 1141 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_xxxz = cbuffer.data(gg_geom_1010_off + 1142 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_xxyy = cbuffer.data(gg_geom_1010_off + 1143 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_xxyz = cbuffer.data(gg_geom_1010_off + 1144 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_xxzz = cbuffer.data(gg_geom_1010_off + 1145 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_xyyy = cbuffer.data(gg_geom_1010_off + 1146 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_xyyz = cbuffer.data(gg_geom_1010_off + 1147 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_xyzz = cbuffer.data(gg_geom_1010_off + 1148 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_xzzz = cbuffer.data(gg_geom_1010_off + 1149 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_yyyy = cbuffer.data(gg_geom_1010_off + 1150 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_yyyz = cbuffer.data(gg_geom_1010_off + 1151 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_yyzz = cbuffer.data(gg_geom_1010_off + 1152 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_yzzz = cbuffer.data(gg_geom_1010_off + 1153 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_zzzz = cbuffer.data(gg_geom_1010_off + 1154 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxy_xxxx, g_y_0_z_0_xxxy_xxxy, g_y_0_z_0_xxxy_xxxz, g_y_0_z_0_xxxy_xxyy, g_y_0_z_0_xxxy_xxyz, g_y_0_z_0_xxxy_xxzz, g_y_0_z_0_xxxy_xyyy, g_y_0_z_0_xxxy_xyyz, g_y_0_z_0_xxxy_xyzz, g_y_0_z_0_xxxy_xzzz, g_y_0_z_0_xxxy_yyyy, g_y_0_z_0_xxxy_yyyz, g_y_0_z_0_xxxy_yyzz, g_y_0_z_0_xxxy_yzzz, g_y_0_z_0_xxxy_zzzz, g_y_0_z_0_xxy_xxxx, g_y_0_z_0_xxy_xxxxx, g_y_0_z_0_xxy_xxxxy, g_y_0_z_0_xxy_xxxxz, g_y_0_z_0_xxy_xxxy, g_y_0_z_0_xxy_xxxyy, g_y_0_z_0_xxy_xxxyz, g_y_0_z_0_xxy_xxxz, g_y_0_z_0_xxy_xxxzz, g_y_0_z_0_xxy_xxyy, g_y_0_z_0_xxy_xxyyy, g_y_0_z_0_xxy_xxyyz, g_y_0_z_0_xxy_xxyz, g_y_0_z_0_xxy_xxyzz, g_y_0_z_0_xxy_xxzz, g_y_0_z_0_xxy_xxzzz, g_y_0_z_0_xxy_xyyy, g_y_0_z_0_xxy_xyyyy, g_y_0_z_0_xxy_xyyyz, g_y_0_z_0_xxy_xyyz, g_y_0_z_0_xxy_xyyzz, g_y_0_z_0_xxy_xyzz, g_y_0_z_0_xxy_xyzzz, g_y_0_z_0_xxy_xzzz, g_y_0_z_0_xxy_xzzzz, g_y_0_z_0_xxy_yyyy, g_y_0_z_0_xxy_yyyz, g_y_0_z_0_xxy_yyzz, g_y_0_z_0_xxy_yzzz, g_y_0_z_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxy_xxxx[k] = -g_y_0_z_0_xxy_xxxx[k] * ab_x + g_y_0_z_0_xxy_xxxxx[k];

                g_y_0_z_0_xxxy_xxxy[k] = -g_y_0_z_0_xxy_xxxy[k] * ab_x + g_y_0_z_0_xxy_xxxxy[k];

                g_y_0_z_0_xxxy_xxxz[k] = -g_y_0_z_0_xxy_xxxz[k] * ab_x + g_y_0_z_0_xxy_xxxxz[k];

                g_y_0_z_0_xxxy_xxyy[k] = -g_y_0_z_0_xxy_xxyy[k] * ab_x + g_y_0_z_0_xxy_xxxyy[k];

                g_y_0_z_0_xxxy_xxyz[k] = -g_y_0_z_0_xxy_xxyz[k] * ab_x + g_y_0_z_0_xxy_xxxyz[k];

                g_y_0_z_0_xxxy_xxzz[k] = -g_y_0_z_0_xxy_xxzz[k] * ab_x + g_y_0_z_0_xxy_xxxzz[k];

                g_y_0_z_0_xxxy_xyyy[k] = -g_y_0_z_0_xxy_xyyy[k] * ab_x + g_y_0_z_0_xxy_xxyyy[k];

                g_y_0_z_0_xxxy_xyyz[k] = -g_y_0_z_0_xxy_xyyz[k] * ab_x + g_y_0_z_0_xxy_xxyyz[k];

                g_y_0_z_0_xxxy_xyzz[k] = -g_y_0_z_0_xxy_xyzz[k] * ab_x + g_y_0_z_0_xxy_xxyzz[k];

                g_y_0_z_0_xxxy_xzzz[k] = -g_y_0_z_0_xxy_xzzz[k] * ab_x + g_y_0_z_0_xxy_xxzzz[k];

                g_y_0_z_0_xxxy_yyyy[k] = -g_y_0_z_0_xxy_yyyy[k] * ab_x + g_y_0_z_0_xxy_xyyyy[k];

                g_y_0_z_0_xxxy_yyyz[k] = -g_y_0_z_0_xxy_yyyz[k] * ab_x + g_y_0_z_0_xxy_xyyyz[k];

                g_y_0_z_0_xxxy_yyzz[k] = -g_y_0_z_0_xxy_yyzz[k] * ab_x + g_y_0_z_0_xxy_xyyzz[k];

                g_y_0_z_0_xxxy_yzzz[k] = -g_y_0_z_0_xxy_yzzz[k] * ab_x + g_y_0_z_0_xxy_xyzzz[k];

                g_y_0_z_0_xxxy_zzzz[k] = -g_y_0_z_0_xxy_zzzz[k] * ab_x + g_y_0_z_0_xxy_xzzzz[k];
            }

            /// Set up 1155-1170 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxz_xxxx = cbuffer.data(gg_geom_1010_off + 1155 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_xxxy = cbuffer.data(gg_geom_1010_off + 1156 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_xxxz = cbuffer.data(gg_geom_1010_off + 1157 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_xxyy = cbuffer.data(gg_geom_1010_off + 1158 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_xxyz = cbuffer.data(gg_geom_1010_off + 1159 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_xxzz = cbuffer.data(gg_geom_1010_off + 1160 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_xyyy = cbuffer.data(gg_geom_1010_off + 1161 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_xyyz = cbuffer.data(gg_geom_1010_off + 1162 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_xyzz = cbuffer.data(gg_geom_1010_off + 1163 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_xzzz = cbuffer.data(gg_geom_1010_off + 1164 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_yyyy = cbuffer.data(gg_geom_1010_off + 1165 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_yyyz = cbuffer.data(gg_geom_1010_off + 1166 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_yyzz = cbuffer.data(gg_geom_1010_off + 1167 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_yzzz = cbuffer.data(gg_geom_1010_off + 1168 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_zzzz = cbuffer.data(gg_geom_1010_off + 1169 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxz_xxxx, g_y_0_z_0_xxxz_xxxy, g_y_0_z_0_xxxz_xxxz, g_y_0_z_0_xxxz_xxyy, g_y_0_z_0_xxxz_xxyz, g_y_0_z_0_xxxz_xxzz, g_y_0_z_0_xxxz_xyyy, g_y_0_z_0_xxxz_xyyz, g_y_0_z_0_xxxz_xyzz, g_y_0_z_0_xxxz_xzzz, g_y_0_z_0_xxxz_yyyy, g_y_0_z_0_xxxz_yyyz, g_y_0_z_0_xxxz_yyzz, g_y_0_z_0_xxxz_yzzz, g_y_0_z_0_xxxz_zzzz, g_y_0_z_0_xxz_xxxx, g_y_0_z_0_xxz_xxxxx, g_y_0_z_0_xxz_xxxxy, g_y_0_z_0_xxz_xxxxz, g_y_0_z_0_xxz_xxxy, g_y_0_z_0_xxz_xxxyy, g_y_0_z_0_xxz_xxxyz, g_y_0_z_0_xxz_xxxz, g_y_0_z_0_xxz_xxxzz, g_y_0_z_0_xxz_xxyy, g_y_0_z_0_xxz_xxyyy, g_y_0_z_0_xxz_xxyyz, g_y_0_z_0_xxz_xxyz, g_y_0_z_0_xxz_xxyzz, g_y_0_z_0_xxz_xxzz, g_y_0_z_0_xxz_xxzzz, g_y_0_z_0_xxz_xyyy, g_y_0_z_0_xxz_xyyyy, g_y_0_z_0_xxz_xyyyz, g_y_0_z_0_xxz_xyyz, g_y_0_z_0_xxz_xyyzz, g_y_0_z_0_xxz_xyzz, g_y_0_z_0_xxz_xyzzz, g_y_0_z_0_xxz_xzzz, g_y_0_z_0_xxz_xzzzz, g_y_0_z_0_xxz_yyyy, g_y_0_z_0_xxz_yyyz, g_y_0_z_0_xxz_yyzz, g_y_0_z_0_xxz_yzzz, g_y_0_z_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxz_xxxx[k] = -g_y_0_z_0_xxz_xxxx[k] * ab_x + g_y_0_z_0_xxz_xxxxx[k];

                g_y_0_z_0_xxxz_xxxy[k] = -g_y_0_z_0_xxz_xxxy[k] * ab_x + g_y_0_z_0_xxz_xxxxy[k];

                g_y_0_z_0_xxxz_xxxz[k] = -g_y_0_z_0_xxz_xxxz[k] * ab_x + g_y_0_z_0_xxz_xxxxz[k];

                g_y_0_z_0_xxxz_xxyy[k] = -g_y_0_z_0_xxz_xxyy[k] * ab_x + g_y_0_z_0_xxz_xxxyy[k];

                g_y_0_z_0_xxxz_xxyz[k] = -g_y_0_z_0_xxz_xxyz[k] * ab_x + g_y_0_z_0_xxz_xxxyz[k];

                g_y_0_z_0_xxxz_xxzz[k] = -g_y_0_z_0_xxz_xxzz[k] * ab_x + g_y_0_z_0_xxz_xxxzz[k];

                g_y_0_z_0_xxxz_xyyy[k] = -g_y_0_z_0_xxz_xyyy[k] * ab_x + g_y_0_z_0_xxz_xxyyy[k];

                g_y_0_z_0_xxxz_xyyz[k] = -g_y_0_z_0_xxz_xyyz[k] * ab_x + g_y_0_z_0_xxz_xxyyz[k];

                g_y_0_z_0_xxxz_xyzz[k] = -g_y_0_z_0_xxz_xyzz[k] * ab_x + g_y_0_z_0_xxz_xxyzz[k];

                g_y_0_z_0_xxxz_xzzz[k] = -g_y_0_z_0_xxz_xzzz[k] * ab_x + g_y_0_z_0_xxz_xxzzz[k];

                g_y_0_z_0_xxxz_yyyy[k] = -g_y_0_z_0_xxz_yyyy[k] * ab_x + g_y_0_z_0_xxz_xyyyy[k];

                g_y_0_z_0_xxxz_yyyz[k] = -g_y_0_z_0_xxz_yyyz[k] * ab_x + g_y_0_z_0_xxz_xyyyz[k];

                g_y_0_z_0_xxxz_yyzz[k] = -g_y_0_z_0_xxz_yyzz[k] * ab_x + g_y_0_z_0_xxz_xyyzz[k];

                g_y_0_z_0_xxxz_yzzz[k] = -g_y_0_z_0_xxz_yzzz[k] * ab_x + g_y_0_z_0_xxz_xyzzz[k];

                g_y_0_z_0_xxxz_zzzz[k] = -g_y_0_z_0_xxz_zzzz[k] * ab_x + g_y_0_z_0_xxz_xzzzz[k];
            }

            /// Set up 1170-1185 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyy_xxxx = cbuffer.data(gg_geom_1010_off + 1170 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_xxxy = cbuffer.data(gg_geom_1010_off + 1171 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_xxxz = cbuffer.data(gg_geom_1010_off + 1172 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_xxyy = cbuffer.data(gg_geom_1010_off + 1173 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_xxyz = cbuffer.data(gg_geom_1010_off + 1174 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_xxzz = cbuffer.data(gg_geom_1010_off + 1175 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_xyyy = cbuffer.data(gg_geom_1010_off + 1176 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_xyyz = cbuffer.data(gg_geom_1010_off + 1177 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_xyzz = cbuffer.data(gg_geom_1010_off + 1178 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_xzzz = cbuffer.data(gg_geom_1010_off + 1179 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_yyyy = cbuffer.data(gg_geom_1010_off + 1180 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_yyyz = cbuffer.data(gg_geom_1010_off + 1181 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_yyzz = cbuffer.data(gg_geom_1010_off + 1182 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_yzzz = cbuffer.data(gg_geom_1010_off + 1183 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_zzzz = cbuffer.data(gg_geom_1010_off + 1184 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyy_xxxx, g_y_0_z_0_xxyy_xxxy, g_y_0_z_0_xxyy_xxxz, g_y_0_z_0_xxyy_xxyy, g_y_0_z_0_xxyy_xxyz, g_y_0_z_0_xxyy_xxzz, g_y_0_z_0_xxyy_xyyy, g_y_0_z_0_xxyy_xyyz, g_y_0_z_0_xxyy_xyzz, g_y_0_z_0_xxyy_xzzz, g_y_0_z_0_xxyy_yyyy, g_y_0_z_0_xxyy_yyyz, g_y_0_z_0_xxyy_yyzz, g_y_0_z_0_xxyy_yzzz, g_y_0_z_0_xxyy_zzzz, g_y_0_z_0_xyy_xxxx, g_y_0_z_0_xyy_xxxxx, g_y_0_z_0_xyy_xxxxy, g_y_0_z_0_xyy_xxxxz, g_y_0_z_0_xyy_xxxy, g_y_0_z_0_xyy_xxxyy, g_y_0_z_0_xyy_xxxyz, g_y_0_z_0_xyy_xxxz, g_y_0_z_0_xyy_xxxzz, g_y_0_z_0_xyy_xxyy, g_y_0_z_0_xyy_xxyyy, g_y_0_z_0_xyy_xxyyz, g_y_0_z_0_xyy_xxyz, g_y_0_z_0_xyy_xxyzz, g_y_0_z_0_xyy_xxzz, g_y_0_z_0_xyy_xxzzz, g_y_0_z_0_xyy_xyyy, g_y_0_z_0_xyy_xyyyy, g_y_0_z_0_xyy_xyyyz, g_y_0_z_0_xyy_xyyz, g_y_0_z_0_xyy_xyyzz, g_y_0_z_0_xyy_xyzz, g_y_0_z_0_xyy_xyzzz, g_y_0_z_0_xyy_xzzz, g_y_0_z_0_xyy_xzzzz, g_y_0_z_0_xyy_yyyy, g_y_0_z_0_xyy_yyyz, g_y_0_z_0_xyy_yyzz, g_y_0_z_0_xyy_yzzz, g_y_0_z_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyy_xxxx[k] = -g_y_0_z_0_xyy_xxxx[k] * ab_x + g_y_0_z_0_xyy_xxxxx[k];

                g_y_0_z_0_xxyy_xxxy[k] = -g_y_0_z_0_xyy_xxxy[k] * ab_x + g_y_0_z_0_xyy_xxxxy[k];

                g_y_0_z_0_xxyy_xxxz[k] = -g_y_0_z_0_xyy_xxxz[k] * ab_x + g_y_0_z_0_xyy_xxxxz[k];

                g_y_0_z_0_xxyy_xxyy[k] = -g_y_0_z_0_xyy_xxyy[k] * ab_x + g_y_0_z_0_xyy_xxxyy[k];

                g_y_0_z_0_xxyy_xxyz[k] = -g_y_0_z_0_xyy_xxyz[k] * ab_x + g_y_0_z_0_xyy_xxxyz[k];

                g_y_0_z_0_xxyy_xxzz[k] = -g_y_0_z_0_xyy_xxzz[k] * ab_x + g_y_0_z_0_xyy_xxxzz[k];

                g_y_0_z_0_xxyy_xyyy[k] = -g_y_0_z_0_xyy_xyyy[k] * ab_x + g_y_0_z_0_xyy_xxyyy[k];

                g_y_0_z_0_xxyy_xyyz[k] = -g_y_0_z_0_xyy_xyyz[k] * ab_x + g_y_0_z_0_xyy_xxyyz[k];

                g_y_0_z_0_xxyy_xyzz[k] = -g_y_0_z_0_xyy_xyzz[k] * ab_x + g_y_0_z_0_xyy_xxyzz[k];

                g_y_0_z_0_xxyy_xzzz[k] = -g_y_0_z_0_xyy_xzzz[k] * ab_x + g_y_0_z_0_xyy_xxzzz[k];

                g_y_0_z_0_xxyy_yyyy[k] = -g_y_0_z_0_xyy_yyyy[k] * ab_x + g_y_0_z_0_xyy_xyyyy[k];

                g_y_0_z_0_xxyy_yyyz[k] = -g_y_0_z_0_xyy_yyyz[k] * ab_x + g_y_0_z_0_xyy_xyyyz[k];

                g_y_0_z_0_xxyy_yyzz[k] = -g_y_0_z_0_xyy_yyzz[k] * ab_x + g_y_0_z_0_xyy_xyyzz[k];

                g_y_0_z_0_xxyy_yzzz[k] = -g_y_0_z_0_xyy_yzzz[k] * ab_x + g_y_0_z_0_xyy_xyzzz[k];

                g_y_0_z_0_xxyy_zzzz[k] = -g_y_0_z_0_xyy_zzzz[k] * ab_x + g_y_0_z_0_xyy_xzzzz[k];
            }

            /// Set up 1185-1200 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyz_xxxx = cbuffer.data(gg_geom_1010_off + 1185 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_xxxy = cbuffer.data(gg_geom_1010_off + 1186 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_xxxz = cbuffer.data(gg_geom_1010_off + 1187 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_xxyy = cbuffer.data(gg_geom_1010_off + 1188 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_xxyz = cbuffer.data(gg_geom_1010_off + 1189 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_xxzz = cbuffer.data(gg_geom_1010_off + 1190 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_xyyy = cbuffer.data(gg_geom_1010_off + 1191 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_xyyz = cbuffer.data(gg_geom_1010_off + 1192 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_xyzz = cbuffer.data(gg_geom_1010_off + 1193 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_xzzz = cbuffer.data(gg_geom_1010_off + 1194 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_yyyy = cbuffer.data(gg_geom_1010_off + 1195 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_yyyz = cbuffer.data(gg_geom_1010_off + 1196 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_yyzz = cbuffer.data(gg_geom_1010_off + 1197 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_yzzz = cbuffer.data(gg_geom_1010_off + 1198 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_zzzz = cbuffer.data(gg_geom_1010_off + 1199 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyz_xxxx, g_y_0_z_0_xxyz_xxxy, g_y_0_z_0_xxyz_xxxz, g_y_0_z_0_xxyz_xxyy, g_y_0_z_0_xxyz_xxyz, g_y_0_z_0_xxyz_xxzz, g_y_0_z_0_xxyz_xyyy, g_y_0_z_0_xxyz_xyyz, g_y_0_z_0_xxyz_xyzz, g_y_0_z_0_xxyz_xzzz, g_y_0_z_0_xxyz_yyyy, g_y_0_z_0_xxyz_yyyz, g_y_0_z_0_xxyz_yyzz, g_y_0_z_0_xxyz_yzzz, g_y_0_z_0_xxyz_zzzz, g_y_0_z_0_xyz_xxxx, g_y_0_z_0_xyz_xxxxx, g_y_0_z_0_xyz_xxxxy, g_y_0_z_0_xyz_xxxxz, g_y_0_z_0_xyz_xxxy, g_y_0_z_0_xyz_xxxyy, g_y_0_z_0_xyz_xxxyz, g_y_0_z_0_xyz_xxxz, g_y_0_z_0_xyz_xxxzz, g_y_0_z_0_xyz_xxyy, g_y_0_z_0_xyz_xxyyy, g_y_0_z_0_xyz_xxyyz, g_y_0_z_0_xyz_xxyz, g_y_0_z_0_xyz_xxyzz, g_y_0_z_0_xyz_xxzz, g_y_0_z_0_xyz_xxzzz, g_y_0_z_0_xyz_xyyy, g_y_0_z_0_xyz_xyyyy, g_y_0_z_0_xyz_xyyyz, g_y_0_z_0_xyz_xyyz, g_y_0_z_0_xyz_xyyzz, g_y_0_z_0_xyz_xyzz, g_y_0_z_0_xyz_xyzzz, g_y_0_z_0_xyz_xzzz, g_y_0_z_0_xyz_xzzzz, g_y_0_z_0_xyz_yyyy, g_y_0_z_0_xyz_yyyz, g_y_0_z_0_xyz_yyzz, g_y_0_z_0_xyz_yzzz, g_y_0_z_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyz_xxxx[k] = -g_y_0_z_0_xyz_xxxx[k] * ab_x + g_y_0_z_0_xyz_xxxxx[k];

                g_y_0_z_0_xxyz_xxxy[k] = -g_y_0_z_0_xyz_xxxy[k] * ab_x + g_y_0_z_0_xyz_xxxxy[k];

                g_y_0_z_0_xxyz_xxxz[k] = -g_y_0_z_0_xyz_xxxz[k] * ab_x + g_y_0_z_0_xyz_xxxxz[k];

                g_y_0_z_0_xxyz_xxyy[k] = -g_y_0_z_0_xyz_xxyy[k] * ab_x + g_y_0_z_0_xyz_xxxyy[k];

                g_y_0_z_0_xxyz_xxyz[k] = -g_y_0_z_0_xyz_xxyz[k] * ab_x + g_y_0_z_0_xyz_xxxyz[k];

                g_y_0_z_0_xxyz_xxzz[k] = -g_y_0_z_0_xyz_xxzz[k] * ab_x + g_y_0_z_0_xyz_xxxzz[k];

                g_y_0_z_0_xxyz_xyyy[k] = -g_y_0_z_0_xyz_xyyy[k] * ab_x + g_y_0_z_0_xyz_xxyyy[k];

                g_y_0_z_0_xxyz_xyyz[k] = -g_y_0_z_0_xyz_xyyz[k] * ab_x + g_y_0_z_0_xyz_xxyyz[k];

                g_y_0_z_0_xxyz_xyzz[k] = -g_y_0_z_0_xyz_xyzz[k] * ab_x + g_y_0_z_0_xyz_xxyzz[k];

                g_y_0_z_0_xxyz_xzzz[k] = -g_y_0_z_0_xyz_xzzz[k] * ab_x + g_y_0_z_0_xyz_xxzzz[k];

                g_y_0_z_0_xxyz_yyyy[k] = -g_y_0_z_0_xyz_yyyy[k] * ab_x + g_y_0_z_0_xyz_xyyyy[k];

                g_y_0_z_0_xxyz_yyyz[k] = -g_y_0_z_0_xyz_yyyz[k] * ab_x + g_y_0_z_0_xyz_xyyyz[k];

                g_y_0_z_0_xxyz_yyzz[k] = -g_y_0_z_0_xyz_yyzz[k] * ab_x + g_y_0_z_0_xyz_xyyzz[k];

                g_y_0_z_0_xxyz_yzzz[k] = -g_y_0_z_0_xyz_yzzz[k] * ab_x + g_y_0_z_0_xyz_xyzzz[k];

                g_y_0_z_0_xxyz_zzzz[k] = -g_y_0_z_0_xyz_zzzz[k] * ab_x + g_y_0_z_0_xyz_xzzzz[k];
            }

            /// Set up 1200-1215 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxzz_xxxx = cbuffer.data(gg_geom_1010_off + 1200 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_xxxy = cbuffer.data(gg_geom_1010_off + 1201 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_xxxz = cbuffer.data(gg_geom_1010_off + 1202 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_xxyy = cbuffer.data(gg_geom_1010_off + 1203 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_xxyz = cbuffer.data(gg_geom_1010_off + 1204 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_xxzz = cbuffer.data(gg_geom_1010_off + 1205 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_xyyy = cbuffer.data(gg_geom_1010_off + 1206 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_xyyz = cbuffer.data(gg_geom_1010_off + 1207 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_xyzz = cbuffer.data(gg_geom_1010_off + 1208 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_xzzz = cbuffer.data(gg_geom_1010_off + 1209 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_yyyy = cbuffer.data(gg_geom_1010_off + 1210 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_yyyz = cbuffer.data(gg_geom_1010_off + 1211 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_yyzz = cbuffer.data(gg_geom_1010_off + 1212 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_yzzz = cbuffer.data(gg_geom_1010_off + 1213 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_zzzz = cbuffer.data(gg_geom_1010_off + 1214 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxzz_xxxx, g_y_0_z_0_xxzz_xxxy, g_y_0_z_0_xxzz_xxxz, g_y_0_z_0_xxzz_xxyy, g_y_0_z_0_xxzz_xxyz, g_y_0_z_0_xxzz_xxzz, g_y_0_z_0_xxzz_xyyy, g_y_0_z_0_xxzz_xyyz, g_y_0_z_0_xxzz_xyzz, g_y_0_z_0_xxzz_xzzz, g_y_0_z_0_xxzz_yyyy, g_y_0_z_0_xxzz_yyyz, g_y_0_z_0_xxzz_yyzz, g_y_0_z_0_xxzz_yzzz, g_y_0_z_0_xxzz_zzzz, g_y_0_z_0_xzz_xxxx, g_y_0_z_0_xzz_xxxxx, g_y_0_z_0_xzz_xxxxy, g_y_0_z_0_xzz_xxxxz, g_y_0_z_0_xzz_xxxy, g_y_0_z_0_xzz_xxxyy, g_y_0_z_0_xzz_xxxyz, g_y_0_z_0_xzz_xxxz, g_y_0_z_0_xzz_xxxzz, g_y_0_z_0_xzz_xxyy, g_y_0_z_0_xzz_xxyyy, g_y_0_z_0_xzz_xxyyz, g_y_0_z_0_xzz_xxyz, g_y_0_z_0_xzz_xxyzz, g_y_0_z_0_xzz_xxzz, g_y_0_z_0_xzz_xxzzz, g_y_0_z_0_xzz_xyyy, g_y_0_z_0_xzz_xyyyy, g_y_0_z_0_xzz_xyyyz, g_y_0_z_0_xzz_xyyz, g_y_0_z_0_xzz_xyyzz, g_y_0_z_0_xzz_xyzz, g_y_0_z_0_xzz_xyzzz, g_y_0_z_0_xzz_xzzz, g_y_0_z_0_xzz_xzzzz, g_y_0_z_0_xzz_yyyy, g_y_0_z_0_xzz_yyyz, g_y_0_z_0_xzz_yyzz, g_y_0_z_0_xzz_yzzz, g_y_0_z_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxzz_xxxx[k] = -g_y_0_z_0_xzz_xxxx[k] * ab_x + g_y_0_z_0_xzz_xxxxx[k];

                g_y_0_z_0_xxzz_xxxy[k] = -g_y_0_z_0_xzz_xxxy[k] * ab_x + g_y_0_z_0_xzz_xxxxy[k];

                g_y_0_z_0_xxzz_xxxz[k] = -g_y_0_z_0_xzz_xxxz[k] * ab_x + g_y_0_z_0_xzz_xxxxz[k];

                g_y_0_z_0_xxzz_xxyy[k] = -g_y_0_z_0_xzz_xxyy[k] * ab_x + g_y_0_z_0_xzz_xxxyy[k];

                g_y_0_z_0_xxzz_xxyz[k] = -g_y_0_z_0_xzz_xxyz[k] * ab_x + g_y_0_z_0_xzz_xxxyz[k];

                g_y_0_z_0_xxzz_xxzz[k] = -g_y_0_z_0_xzz_xxzz[k] * ab_x + g_y_0_z_0_xzz_xxxzz[k];

                g_y_0_z_0_xxzz_xyyy[k] = -g_y_0_z_0_xzz_xyyy[k] * ab_x + g_y_0_z_0_xzz_xxyyy[k];

                g_y_0_z_0_xxzz_xyyz[k] = -g_y_0_z_0_xzz_xyyz[k] * ab_x + g_y_0_z_0_xzz_xxyyz[k];

                g_y_0_z_0_xxzz_xyzz[k] = -g_y_0_z_0_xzz_xyzz[k] * ab_x + g_y_0_z_0_xzz_xxyzz[k];

                g_y_0_z_0_xxzz_xzzz[k] = -g_y_0_z_0_xzz_xzzz[k] * ab_x + g_y_0_z_0_xzz_xxzzz[k];

                g_y_0_z_0_xxzz_yyyy[k] = -g_y_0_z_0_xzz_yyyy[k] * ab_x + g_y_0_z_0_xzz_xyyyy[k];

                g_y_0_z_0_xxzz_yyyz[k] = -g_y_0_z_0_xzz_yyyz[k] * ab_x + g_y_0_z_0_xzz_xyyyz[k];

                g_y_0_z_0_xxzz_yyzz[k] = -g_y_0_z_0_xzz_yyzz[k] * ab_x + g_y_0_z_0_xzz_xyyzz[k];

                g_y_0_z_0_xxzz_yzzz[k] = -g_y_0_z_0_xzz_yzzz[k] * ab_x + g_y_0_z_0_xzz_xyzzz[k];

                g_y_0_z_0_xxzz_zzzz[k] = -g_y_0_z_0_xzz_zzzz[k] * ab_x + g_y_0_z_0_xzz_xzzzz[k];
            }

            /// Set up 1215-1230 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyy_xxxx = cbuffer.data(gg_geom_1010_off + 1215 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_xxxy = cbuffer.data(gg_geom_1010_off + 1216 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_xxxz = cbuffer.data(gg_geom_1010_off + 1217 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_xxyy = cbuffer.data(gg_geom_1010_off + 1218 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_xxyz = cbuffer.data(gg_geom_1010_off + 1219 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_xxzz = cbuffer.data(gg_geom_1010_off + 1220 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_xyyy = cbuffer.data(gg_geom_1010_off + 1221 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_xyyz = cbuffer.data(gg_geom_1010_off + 1222 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_xyzz = cbuffer.data(gg_geom_1010_off + 1223 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_xzzz = cbuffer.data(gg_geom_1010_off + 1224 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1225 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1226 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1227 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1228 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1229 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyy_xxxx, g_y_0_z_0_xyyy_xxxy, g_y_0_z_0_xyyy_xxxz, g_y_0_z_0_xyyy_xxyy, g_y_0_z_0_xyyy_xxyz, g_y_0_z_0_xyyy_xxzz, g_y_0_z_0_xyyy_xyyy, g_y_0_z_0_xyyy_xyyz, g_y_0_z_0_xyyy_xyzz, g_y_0_z_0_xyyy_xzzz, g_y_0_z_0_xyyy_yyyy, g_y_0_z_0_xyyy_yyyz, g_y_0_z_0_xyyy_yyzz, g_y_0_z_0_xyyy_yzzz, g_y_0_z_0_xyyy_zzzz, g_y_0_z_0_yyy_xxxx, g_y_0_z_0_yyy_xxxxx, g_y_0_z_0_yyy_xxxxy, g_y_0_z_0_yyy_xxxxz, g_y_0_z_0_yyy_xxxy, g_y_0_z_0_yyy_xxxyy, g_y_0_z_0_yyy_xxxyz, g_y_0_z_0_yyy_xxxz, g_y_0_z_0_yyy_xxxzz, g_y_0_z_0_yyy_xxyy, g_y_0_z_0_yyy_xxyyy, g_y_0_z_0_yyy_xxyyz, g_y_0_z_0_yyy_xxyz, g_y_0_z_0_yyy_xxyzz, g_y_0_z_0_yyy_xxzz, g_y_0_z_0_yyy_xxzzz, g_y_0_z_0_yyy_xyyy, g_y_0_z_0_yyy_xyyyy, g_y_0_z_0_yyy_xyyyz, g_y_0_z_0_yyy_xyyz, g_y_0_z_0_yyy_xyyzz, g_y_0_z_0_yyy_xyzz, g_y_0_z_0_yyy_xyzzz, g_y_0_z_0_yyy_xzzz, g_y_0_z_0_yyy_xzzzz, g_y_0_z_0_yyy_yyyy, g_y_0_z_0_yyy_yyyz, g_y_0_z_0_yyy_yyzz, g_y_0_z_0_yyy_yzzz, g_y_0_z_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyy_xxxx[k] = -g_y_0_z_0_yyy_xxxx[k] * ab_x + g_y_0_z_0_yyy_xxxxx[k];

                g_y_0_z_0_xyyy_xxxy[k] = -g_y_0_z_0_yyy_xxxy[k] * ab_x + g_y_0_z_0_yyy_xxxxy[k];

                g_y_0_z_0_xyyy_xxxz[k] = -g_y_0_z_0_yyy_xxxz[k] * ab_x + g_y_0_z_0_yyy_xxxxz[k];

                g_y_0_z_0_xyyy_xxyy[k] = -g_y_0_z_0_yyy_xxyy[k] * ab_x + g_y_0_z_0_yyy_xxxyy[k];

                g_y_0_z_0_xyyy_xxyz[k] = -g_y_0_z_0_yyy_xxyz[k] * ab_x + g_y_0_z_0_yyy_xxxyz[k];

                g_y_0_z_0_xyyy_xxzz[k] = -g_y_0_z_0_yyy_xxzz[k] * ab_x + g_y_0_z_0_yyy_xxxzz[k];

                g_y_0_z_0_xyyy_xyyy[k] = -g_y_0_z_0_yyy_xyyy[k] * ab_x + g_y_0_z_0_yyy_xxyyy[k];

                g_y_0_z_0_xyyy_xyyz[k] = -g_y_0_z_0_yyy_xyyz[k] * ab_x + g_y_0_z_0_yyy_xxyyz[k];

                g_y_0_z_0_xyyy_xyzz[k] = -g_y_0_z_0_yyy_xyzz[k] * ab_x + g_y_0_z_0_yyy_xxyzz[k];

                g_y_0_z_0_xyyy_xzzz[k] = -g_y_0_z_0_yyy_xzzz[k] * ab_x + g_y_0_z_0_yyy_xxzzz[k];

                g_y_0_z_0_xyyy_yyyy[k] = -g_y_0_z_0_yyy_yyyy[k] * ab_x + g_y_0_z_0_yyy_xyyyy[k];

                g_y_0_z_0_xyyy_yyyz[k] = -g_y_0_z_0_yyy_yyyz[k] * ab_x + g_y_0_z_0_yyy_xyyyz[k];

                g_y_0_z_0_xyyy_yyzz[k] = -g_y_0_z_0_yyy_yyzz[k] * ab_x + g_y_0_z_0_yyy_xyyzz[k];

                g_y_0_z_0_xyyy_yzzz[k] = -g_y_0_z_0_yyy_yzzz[k] * ab_x + g_y_0_z_0_yyy_xyzzz[k];

                g_y_0_z_0_xyyy_zzzz[k] = -g_y_0_z_0_yyy_zzzz[k] * ab_x + g_y_0_z_0_yyy_xzzzz[k];
            }

            /// Set up 1230-1245 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1230 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1231 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1232 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1233 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1234 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1235 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1236 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1237 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1238 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1239 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1240 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1241 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1242 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1243 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1244 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyz_xxxx, g_y_0_z_0_xyyz_xxxy, g_y_0_z_0_xyyz_xxxz, g_y_0_z_0_xyyz_xxyy, g_y_0_z_0_xyyz_xxyz, g_y_0_z_0_xyyz_xxzz, g_y_0_z_0_xyyz_xyyy, g_y_0_z_0_xyyz_xyyz, g_y_0_z_0_xyyz_xyzz, g_y_0_z_0_xyyz_xzzz, g_y_0_z_0_xyyz_yyyy, g_y_0_z_0_xyyz_yyyz, g_y_0_z_0_xyyz_yyzz, g_y_0_z_0_xyyz_yzzz, g_y_0_z_0_xyyz_zzzz, g_y_0_z_0_yyz_xxxx, g_y_0_z_0_yyz_xxxxx, g_y_0_z_0_yyz_xxxxy, g_y_0_z_0_yyz_xxxxz, g_y_0_z_0_yyz_xxxy, g_y_0_z_0_yyz_xxxyy, g_y_0_z_0_yyz_xxxyz, g_y_0_z_0_yyz_xxxz, g_y_0_z_0_yyz_xxxzz, g_y_0_z_0_yyz_xxyy, g_y_0_z_0_yyz_xxyyy, g_y_0_z_0_yyz_xxyyz, g_y_0_z_0_yyz_xxyz, g_y_0_z_0_yyz_xxyzz, g_y_0_z_0_yyz_xxzz, g_y_0_z_0_yyz_xxzzz, g_y_0_z_0_yyz_xyyy, g_y_0_z_0_yyz_xyyyy, g_y_0_z_0_yyz_xyyyz, g_y_0_z_0_yyz_xyyz, g_y_0_z_0_yyz_xyyzz, g_y_0_z_0_yyz_xyzz, g_y_0_z_0_yyz_xyzzz, g_y_0_z_0_yyz_xzzz, g_y_0_z_0_yyz_xzzzz, g_y_0_z_0_yyz_yyyy, g_y_0_z_0_yyz_yyyz, g_y_0_z_0_yyz_yyzz, g_y_0_z_0_yyz_yzzz, g_y_0_z_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyz_xxxx[k] = -g_y_0_z_0_yyz_xxxx[k] * ab_x + g_y_0_z_0_yyz_xxxxx[k];

                g_y_0_z_0_xyyz_xxxy[k] = -g_y_0_z_0_yyz_xxxy[k] * ab_x + g_y_0_z_0_yyz_xxxxy[k];

                g_y_0_z_0_xyyz_xxxz[k] = -g_y_0_z_0_yyz_xxxz[k] * ab_x + g_y_0_z_0_yyz_xxxxz[k];

                g_y_0_z_0_xyyz_xxyy[k] = -g_y_0_z_0_yyz_xxyy[k] * ab_x + g_y_0_z_0_yyz_xxxyy[k];

                g_y_0_z_0_xyyz_xxyz[k] = -g_y_0_z_0_yyz_xxyz[k] * ab_x + g_y_0_z_0_yyz_xxxyz[k];

                g_y_0_z_0_xyyz_xxzz[k] = -g_y_0_z_0_yyz_xxzz[k] * ab_x + g_y_0_z_0_yyz_xxxzz[k];

                g_y_0_z_0_xyyz_xyyy[k] = -g_y_0_z_0_yyz_xyyy[k] * ab_x + g_y_0_z_0_yyz_xxyyy[k];

                g_y_0_z_0_xyyz_xyyz[k] = -g_y_0_z_0_yyz_xyyz[k] * ab_x + g_y_0_z_0_yyz_xxyyz[k];

                g_y_0_z_0_xyyz_xyzz[k] = -g_y_0_z_0_yyz_xyzz[k] * ab_x + g_y_0_z_0_yyz_xxyzz[k];

                g_y_0_z_0_xyyz_xzzz[k] = -g_y_0_z_0_yyz_xzzz[k] * ab_x + g_y_0_z_0_yyz_xxzzz[k];

                g_y_0_z_0_xyyz_yyyy[k] = -g_y_0_z_0_yyz_yyyy[k] * ab_x + g_y_0_z_0_yyz_xyyyy[k];

                g_y_0_z_0_xyyz_yyyz[k] = -g_y_0_z_0_yyz_yyyz[k] * ab_x + g_y_0_z_0_yyz_xyyyz[k];

                g_y_0_z_0_xyyz_yyzz[k] = -g_y_0_z_0_yyz_yyzz[k] * ab_x + g_y_0_z_0_yyz_xyyzz[k];

                g_y_0_z_0_xyyz_yzzz[k] = -g_y_0_z_0_yyz_yzzz[k] * ab_x + g_y_0_z_0_yyz_xyzzz[k];

                g_y_0_z_0_xyyz_zzzz[k] = -g_y_0_z_0_yyz_zzzz[k] * ab_x + g_y_0_z_0_yyz_xzzzz[k];
            }

            /// Set up 1245-1260 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1245 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1246 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1247 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1248 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1249 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1250 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1251 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1252 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1253 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1254 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1255 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1256 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1257 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1258 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyzz_xxxx, g_y_0_z_0_xyzz_xxxy, g_y_0_z_0_xyzz_xxxz, g_y_0_z_0_xyzz_xxyy, g_y_0_z_0_xyzz_xxyz, g_y_0_z_0_xyzz_xxzz, g_y_0_z_0_xyzz_xyyy, g_y_0_z_0_xyzz_xyyz, g_y_0_z_0_xyzz_xyzz, g_y_0_z_0_xyzz_xzzz, g_y_0_z_0_xyzz_yyyy, g_y_0_z_0_xyzz_yyyz, g_y_0_z_0_xyzz_yyzz, g_y_0_z_0_xyzz_yzzz, g_y_0_z_0_xyzz_zzzz, g_y_0_z_0_yzz_xxxx, g_y_0_z_0_yzz_xxxxx, g_y_0_z_0_yzz_xxxxy, g_y_0_z_0_yzz_xxxxz, g_y_0_z_0_yzz_xxxy, g_y_0_z_0_yzz_xxxyy, g_y_0_z_0_yzz_xxxyz, g_y_0_z_0_yzz_xxxz, g_y_0_z_0_yzz_xxxzz, g_y_0_z_0_yzz_xxyy, g_y_0_z_0_yzz_xxyyy, g_y_0_z_0_yzz_xxyyz, g_y_0_z_0_yzz_xxyz, g_y_0_z_0_yzz_xxyzz, g_y_0_z_0_yzz_xxzz, g_y_0_z_0_yzz_xxzzz, g_y_0_z_0_yzz_xyyy, g_y_0_z_0_yzz_xyyyy, g_y_0_z_0_yzz_xyyyz, g_y_0_z_0_yzz_xyyz, g_y_0_z_0_yzz_xyyzz, g_y_0_z_0_yzz_xyzz, g_y_0_z_0_yzz_xyzzz, g_y_0_z_0_yzz_xzzz, g_y_0_z_0_yzz_xzzzz, g_y_0_z_0_yzz_yyyy, g_y_0_z_0_yzz_yyyz, g_y_0_z_0_yzz_yyzz, g_y_0_z_0_yzz_yzzz, g_y_0_z_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyzz_xxxx[k] = -g_y_0_z_0_yzz_xxxx[k] * ab_x + g_y_0_z_0_yzz_xxxxx[k];

                g_y_0_z_0_xyzz_xxxy[k] = -g_y_0_z_0_yzz_xxxy[k] * ab_x + g_y_0_z_0_yzz_xxxxy[k];

                g_y_0_z_0_xyzz_xxxz[k] = -g_y_0_z_0_yzz_xxxz[k] * ab_x + g_y_0_z_0_yzz_xxxxz[k];

                g_y_0_z_0_xyzz_xxyy[k] = -g_y_0_z_0_yzz_xxyy[k] * ab_x + g_y_0_z_0_yzz_xxxyy[k];

                g_y_0_z_0_xyzz_xxyz[k] = -g_y_0_z_0_yzz_xxyz[k] * ab_x + g_y_0_z_0_yzz_xxxyz[k];

                g_y_0_z_0_xyzz_xxzz[k] = -g_y_0_z_0_yzz_xxzz[k] * ab_x + g_y_0_z_0_yzz_xxxzz[k];

                g_y_0_z_0_xyzz_xyyy[k] = -g_y_0_z_0_yzz_xyyy[k] * ab_x + g_y_0_z_0_yzz_xxyyy[k];

                g_y_0_z_0_xyzz_xyyz[k] = -g_y_0_z_0_yzz_xyyz[k] * ab_x + g_y_0_z_0_yzz_xxyyz[k];

                g_y_0_z_0_xyzz_xyzz[k] = -g_y_0_z_0_yzz_xyzz[k] * ab_x + g_y_0_z_0_yzz_xxyzz[k];

                g_y_0_z_0_xyzz_xzzz[k] = -g_y_0_z_0_yzz_xzzz[k] * ab_x + g_y_0_z_0_yzz_xxzzz[k];

                g_y_0_z_0_xyzz_yyyy[k] = -g_y_0_z_0_yzz_yyyy[k] * ab_x + g_y_0_z_0_yzz_xyyyy[k];

                g_y_0_z_0_xyzz_yyyz[k] = -g_y_0_z_0_yzz_yyyz[k] * ab_x + g_y_0_z_0_yzz_xyyyz[k];

                g_y_0_z_0_xyzz_yyzz[k] = -g_y_0_z_0_yzz_yyzz[k] * ab_x + g_y_0_z_0_yzz_xyyzz[k];

                g_y_0_z_0_xyzz_yzzz[k] = -g_y_0_z_0_yzz_yzzz[k] * ab_x + g_y_0_z_0_yzz_xyzzz[k];

                g_y_0_z_0_xyzz_zzzz[k] = -g_y_0_z_0_yzz_zzzz[k] * ab_x + g_y_0_z_0_yzz_xzzzz[k];
            }

            /// Set up 1260-1275 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1260 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1261 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1262 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1263 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1264 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1265 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1266 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1267 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1268 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1269 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1270 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1271 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1272 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1273 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1274 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xzzz_xxxx, g_y_0_z_0_xzzz_xxxy, g_y_0_z_0_xzzz_xxxz, g_y_0_z_0_xzzz_xxyy, g_y_0_z_0_xzzz_xxyz, g_y_0_z_0_xzzz_xxzz, g_y_0_z_0_xzzz_xyyy, g_y_0_z_0_xzzz_xyyz, g_y_0_z_0_xzzz_xyzz, g_y_0_z_0_xzzz_xzzz, g_y_0_z_0_xzzz_yyyy, g_y_0_z_0_xzzz_yyyz, g_y_0_z_0_xzzz_yyzz, g_y_0_z_0_xzzz_yzzz, g_y_0_z_0_xzzz_zzzz, g_y_0_z_0_zzz_xxxx, g_y_0_z_0_zzz_xxxxx, g_y_0_z_0_zzz_xxxxy, g_y_0_z_0_zzz_xxxxz, g_y_0_z_0_zzz_xxxy, g_y_0_z_0_zzz_xxxyy, g_y_0_z_0_zzz_xxxyz, g_y_0_z_0_zzz_xxxz, g_y_0_z_0_zzz_xxxzz, g_y_0_z_0_zzz_xxyy, g_y_0_z_0_zzz_xxyyy, g_y_0_z_0_zzz_xxyyz, g_y_0_z_0_zzz_xxyz, g_y_0_z_0_zzz_xxyzz, g_y_0_z_0_zzz_xxzz, g_y_0_z_0_zzz_xxzzz, g_y_0_z_0_zzz_xyyy, g_y_0_z_0_zzz_xyyyy, g_y_0_z_0_zzz_xyyyz, g_y_0_z_0_zzz_xyyz, g_y_0_z_0_zzz_xyyzz, g_y_0_z_0_zzz_xyzz, g_y_0_z_0_zzz_xyzzz, g_y_0_z_0_zzz_xzzz, g_y_0_z_0_zzz_xzzzz, g_y_0_z_0_zzz_yyyy, g_y_0_z_0_zzz_yyyz, g_y_0_z_0_zzz_yyzz, g_y_0_z_0_zzz_yzzz, g_y_0_z_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xzzz_xxxx[k] = -g_y_0_z_0_zzz_xxxx[k] * ab_x + g_y_0_z_0_zzz_xxxxx[k];

                g_y_0_z_0_xzzz_xxxy[k] = -g_y_0_z_0_zzz_xxxy[k] * ab_x + g_y_0_z_0_zzz_xxxxy[k];

                g_y_0_z_0_xzzz_xxxz[k] = -g_y_0_z_0_zzz_xxxz[k] * ab_x + g_y_0_z_0_zzz_xxxxz[k];

                g_y_0_z_0_xzzz_xxyy[k] = -g_y_0_z_0_zzz_xxyy[k] * ab_x + g_y_0_z_0_zzz_xxxyy[k];

                g_y_0_z_0_xzzz_xxyz[k] = -g_y_0_z_0_zzz_xxyz[k] * ab_x + g_y_0_z_0_zzz_xxxyz[k];

                g_y_0_z_0_xzzz_xxzz[k] = -g_y_0_z_0_zzz_xxzz[k] * ab_x + g_y_0_z_0_zzz_xxxzz[k];

                g_y_0_z_0_xzzz_xyyy[k] = -g_y_0_z_0_zzz_xyyy[k] * ab_x + g_y_0_z_0_zzz_xxyyy[k];

                g_y_0_z_0_xzzz_xyyz[k] = -g_y_0_z_0_zzz_xyyz[k] * ab_x + g_y_0_z_0_zzz_xxyyz[k];

                g_y_0_z_0_xzzz_xyzz[k] = -g_y_0_z_0_zzz_xyzz[k] * ab_x + g_y_0_z_0_zzz_xxyzz[k];

                g_y_0_z_0_xzzz_xzzz[k] = -g_y_0_z_0_zzz_xzzz[k] * ab_x + g_y_0_z_0_zzz_xxzzz[k];

                g_y_0_z_0_xzzz_yyyy[k] = -g_y_0_z_0_zzz_yyyy[k] * ab_x + g_y_0_z_0_zzz_xyyyy[k];

                g_y_0_z_0_xzzz_yyyz[k] = -g_y_0_z_0_zzz_yyyz[k] * ab_x + g_y_0_z_0_zzz_xyyyz[k];

                g_y_0_z_0_xzzz_yyzz[k] = -g_y_0_z_0_zzz_yyzz[k] * ab_x + g_y_0_z_0_zzz_xyyzz[k];

                g_y_0_z_0_xzzz_yzzz[k] = -g_y_0_z_0_zzz_yzzz[k] * ab_x + g_y_0_z_0_zzz_xyzzz[k];

                g_y_0_z_0_xzzz_zzzz[k] = -g_y_0_z_0_zzz_zzzz[k] * ab_x + g_y_0_z_0_zzz_xzzzz[k];
            }

            /// Set up 1275-1290 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyy_xxxx = cbuffer.data(gg_geom_1010_off + 1275 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_xxxy = cbuffer.data(gg_geom_1010_off + 1276 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_xxxz = cbuffer.data(gg_geom_1010_off + 1277 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_xxyy = cbuffer.data(gg_geom_1010_off + 1278 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_xxyz = cbuffer.data(gg_geom_1010_off + 1279 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_xxzz = cbuffer.data(gg_geom_1010_off + 1280 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_xyyy = cbuffer.data(gg_geom_1010_off + 1281 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_xyyz = cbuffer.data(gg_geom_1010_off + 1282 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_xyzz = cbuffer.data(gg_geom_1010_off + 1283 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_xzzz = cbuffer.data(gg_geom_1010_off + 1284 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1285 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1286 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1287 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1288 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1289 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_yyy_xxxx, g_0_0_z_0_yyy_xxxy, g_0_0_z_0_yyy_xxxz, g_0_0_z_0_yyy_xxyy, g_0_0_z_0_yyy_xxyz, g_0_0_z_0_yyy_xxzz, g_0_0_z_0_yyy_xyyy, g_0_0_z_0_yyy_xyyz, g_0_0_z_0_yyy_xyzz, g_0_0_z_0_yyy_xzzz, g_0_0_z_0_yyy_yyyy, g_0_0_z_0_yyy_yyyz, g_0_0_z_0_yyy_yyzz, g_0_0_z_0_yyy_yzzz, g_0_0_z_0_yyy_zzzz, g_y_0_z_0_yyy_xxxx, g_y_0_z_0_yyy_xxxxy, g_y_0_z_0_yyy_xxxy, g_y_0_z_0_yyy_xxxyy, g_y_0_z_0_yyy_xxxyz, g_y_0_z_0_yyy_xxxz, g_y_0_z_0_yyy_xxyy, g_y_0_z_0_yyy_xxyyy, g_y_0_z_0_yyy_xxyyz, g_y_0_z_0_yyy_xxyz, g_y_0_z_0_yyy_xxyzz, g_y_0_z_0_yyy_xxzz, g_y_0_z_0_yyy_xyyy, g_y_0_z_0_yyy_xyyyy, g_y_0_z_0_yyy_xyyyz, g_y_0_z_0_yyy_xyyz, g_y_0_z_0_yyy_xyyzz, g_y_0_z_0_yyy_xyzz, g_y_0_z_0_yyy_xyzzz, g_y_0_z_0_yyy_xzzz, g_y_0_z_0_yyy_yyyy, g_y_0_z_0_yyy_yyyyy, g_y_0_z_0_yyy_yyyyz, g_y_0_z_0_yyy_yyyz, g_y_0_z_0_yyy_yyyzz, g_y_0_z_0_yyy_yyzz, g_y_0_z_0_yyy_yyzzz, g_y_0_z_0_yyy_yzzz, g_y_0_z_0_yyy_yzzzz, g_y_0_z_0_yyy_zzzz, g_y_0_z_0_yyyy_xxxx, g_y_0_z_0_yyyy_xxxy, g_y_0_z_0_yyyy_xxxz, g_y_0_z_0_yyyy_xxyy, g_y_0_z_0_yyyy_xxyz, g_y_0_z_0_yyyy_xxzz, g_y_0_z_0_yyyy_xyyy, g_y_0_z_0_yyyy_xyyz, g_y_0_z_0_yyyy_xyzz, g_y_0_z_0_yyyy_xzzz, g_y_0_z_0_yyyy_yyyy, g_y_0_z_0_yyyy_yyyz, g_y_0_z_0_yyyy_yyzz, g_y_0_z_0_yyyy_yzzz, g_y_0_z_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyy_xxxx[k] = -g_0_0_z_0_yyy_xxxx[k] - g_y_0_z_0_yyy_xxxx[k] * ab_y + g_y_0_z_0_yyy_xxxxy[k];

                g_y_0_z_0_yyyy_xxxy[k] = -g_0_0_z_0_yyy_xxxy[k] - g_y_0_z_0_yyy_xxxy[k] * ab_y + g_y_0_z_0_yyy_xxxyy[k];

                g_y_0_z_0_yyyy_xxxz[k] = -g_0_0_z_0_yyy_xxxz[k] - g_y_0_z_0_yyy_xxxz[k] * ab_y + g_y_0_z_0_yyy_xxxyz[k];

                g_y_0_z_0_yyyy_xxyy[k] = -g_0_0_z_0_yyy_xxyy[k] - g_y_0_z_0_yyy_xxyy[k] * ab_y + g_y_0_z_0_yyy_xxyyy[k];

                g_y_0_z_0_yyyy_xxyz[k] = -g_0_0_z_0_yyy_xxyz[k] - g_y_0_z_0_yyy_xxyz[k] * ab_y + g_y_0_z_0_yyy_xxyyz[k];

                g_y_0_z_0_yyyy_xxzz[k] = -g_0_0_z_0_yyy_xxzz[k] - g_y_0_z_0_yyy_xxzz[k] * ab_y + g_y_0_z_0_yyy_xxyzz[k];

                g_y_0_z_0_yyyy_xyyy[k] = -g_0_0_z_0_yyy_xyyy[k] - g_y_0_z_0_yyy_xyyy[k] * ab_y + g_y_0_z_0_yyy_xyyyy[k];

                g_y_0_z_0_yyyy_xyyz[k] = -g_0_0_z_0_yyy_xyyz[k] - g_y_0_z_0_yyy_xyyz[k] * ab_y + g_y_0_z_0_yyy_xyyyz[k];

                g_y_0_z_0_yyyy_xyzz[k] = -g_0_0_z_0_yyy_xyzz[k] - g_y_0_z_0_yyy_xyzz[k] * ab_y + g_y_0_z_0_yyy_xyyzz[k];

                g_y_0_z_0_yyyy_xzzz[k] = -g_0_0_z_0_yyy_xzzz[k] - g_y_0_z_0_yyy_xzzz[k] * ab_y + g_y_0_z_0_yyy_xyzzz[k];

                g_y_0_z_0_yyyy_yyyy[k] = -g_0_0_z_0_yyy_yyyy[k] - g_y_0_z_0_yyy_yyyy[k] * ab_y + g_y_0_z_0_yyy_yyyyy[k];

                g_y_0_z_0_yyyy_yyyz[k] = -g_0_0_z_0_yyy_yyyz[k] - g_y_0_z_0_yyy_yyyz[k] * ab_y + g_y_0_z_0_yyy_yyyyz[k];

                g_y_0_z_0_yyyy_yyzz[k] = -g_0_0_z_0_yyy_yyzz[k] - g_y_0_z_0_yyy_yyzz[k] * ab_y + g_y_0_z_0_yyy_yyyzz[k];

                g_y_0_z_0_yyyy_yzzz[k] = -g_0_0_z_0_yyy_yzzz[k] - g_y_0_z_0_yyy_yzzz[k] * ab_y + g_y_0_z_0_yyy_yyzzz[k];

                g_y_0_z_0_yyyy_zzzz[k] = -g_0_0_z_0_yyy_zzzz[k] - g_y_0_z_0_yyy_zzzz[k] * ab_y + g_y_0_z_0_yyy_yzzzz[k];
            }

            /// Set up 1290-1305 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1290 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1291 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1292 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1293 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1294 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1295 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1296 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1297 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1298 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1299 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1300 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1301 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1302 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1303 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1304 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyy_xxxx, g_y_0_z_0_yyy_xxxxz, g_y_0_z_0_yyy_xxxy, g_y_0_z_0_yyy_xxxyz, g_y_0_z_0_yyy_xxxz, g_y_0_z_0_yyy_xxxzz, g_y_0_z_0_yyy_xxyy, g_y_0_z_0_yyy_xxyyz, g_y_0_z_0_yyy_xxyz, g_y_0_z_0_yyy_xxyzz, g_y_0_z_0_yyy_xxzz, g_y_0_z_0_yyy_xxzzz, g_y_0_z_0_yyy_xyyy, g_y_0_z_0_yyy_xyyyz, g_y_0_z_0_yyy_xyyz, g_y_0_z_0_yyy_xyyzz, g_y_0_z_0_yyy_xyzz, g_y_0_z_0_yyy_xyzzz, g_y_0_z_0_yyy_xzzz, g_y_0_z_0_yyy_xzzzz, g_y_0_z_0_yyy_yyyy, g_y_0_z_0_yyy_yyyyz, g_y_0_z_0_yyy_yyyz, g_y_0_z_0_yyy_yyyzz, g_y_0_z_0_yyy_yyzz, g_y_0_z_0_yyy_yyzzz, g_y_0_z_0_yyy_yzzz, g_y_0_z_0_yyy_yzzzz, g_y_0_z_0_yyy_zzzz, g_y_0_z_0_yyy_zzzzz, g_y_0_z_0_yyyz_xxxx, g_y_0_z_0_yyyz_xxxy, g_y_0_z_0_yyyz_xxxz, g_y_0_z_0_yyyz_xxyy, g_y_0_z_0_yyyz_xxyz, g_y_0_z_0_yyyz_xxzz, g_y_0_z_0_yyyz_xyyy, g_y_0_z_0_yyyz_xyyz, g_y_0_z_0_yyyz_xyzz, g_y_0_z_0_yyyz_xzzz, g_y_0_z_0_yyyz_yyyy, g_y_0_z_0_yyyz_yyyz, g_y_0_z_0_yyyz_yyzz, g_y_0_z_0_yyyz_yzzz, g_y_0_z_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyz_xxxx[k] = -g_y_0_z_0_yyy_xxxx[k] * ab_z + g_y_0_z_0_yyy_xxxxz[k];

                g_y_0_z_0_yyyz_xxxy[k] = -g_y_0_z_0_yyy_xxxy[k] * ab_z + g_y_0_z_0_yyy_xxxyz[k];

                g_y_0_z_0_yyyz_xxxz[k] = -g_y_0_z_0_yyy_xxxz[k] * ab_z + g_y_0_z_0_yyy_xxxzz[k];

                g_y_0_z_0_yyyz_xxyy[k] = -g_y_0_z_0_yyy_xxyy[k] * ab_z + g_y_0_z_0_yyy_xxyyz[k];

                g_y_0_z_0_yyyz_xxyz[k] = -g_y_0_z_0_yyy_xxyz[k] * ab_z + g_y_0_z_0_yyy_xxyzz[k];

                g_y_0_z_0_yyyz_xxzz[k] = -g_y_0_z_0_yyy_xxzz[k] * ab_z + g_y_0_z_0_yyy_xxzzz[k];

                g_y_0_z_0_yyyz_xyyy[k] = -g_y_0_z_0_yyy_xyyy[k] * ab_z + g_y_0_z_0_yyy_xyyyz[k];

                g_y_0_z_0_yyyz_xyyz[k] = -g_y_0_z_0_yyy_xyyz[k] * ab_z + g_y_0_z_0_yyy_xyyzz[k];

                g_y_0_z_0_yyyz_xyzz[k] = -g_y_0_z_0_yyy_xyzz[k] * ab_z + g_y_0_z_0_yyy_xyzzz[k];

                g_y_0_z_0_yyyz_xzzz[k] = -g_y_0_z_0_yyy_xzzz[k] * ab_z + g_y_0_z_0_yyy_xzzzz[k];

                g_y_0_z_0_yyyz_yyyy[k] = -g_y_0_z_0_yyy_yyyy[k] * ab_z + g_y_0_z_0_yyy_yyyyz[k];

                g_y_0_z_0_yyyz_yyyz[k] = -g_y_0_z_0_yyy_yyyz[k] * ab_z + g_y_0_z_0_yyy_yyyzz[k];

                g_y_0_z_0_yyyz_yyzz[k] = -g_y_0_z_0_yyy_yyzz[k] * ab_z + g_y_0_z_0_yyy_yyzzz[k];

                g_y_0_z_0_yyyz_yzzz[k] = -g_y_0_z_0_yyy_yzzz[k] * ab_z + g_y_0_z_0_yyy_yzzzz[k];

                g_y_0_z_0_yyyz_zzzz[k] = -g_y_0_z_0_yyy_zzzz[k] * ab_z + g_y_0_z_0_yyy_zzzzz[k];
            }

            /// Set up 1305-1320 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1305 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1306 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1307 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1308 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1309 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1310 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1311 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1312 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1313 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1314 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1315 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1316 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1317 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1318 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1319 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyz_xxxx, g_y_0_z_0_yyz_xxxxz, g_y_0_z_0_yyz_xxxy, g_y_0_z_0_yyz_xxxyz, g_y_0_z_0_yyz_xxxz, g_y_0_z_0_yyz_xxxzz, g_y_0_z_0_yyz_xxyy, g_y_0_z_0_yyz_xxyyz, g_y_0_z_0_yyz_xxyz, g_y_0_z_0_yyz_xxyzz, g_y_0_z_0_yyz_xxzz, g_y_0_z_0_yyz_xxzzz, g_y_0_z_0_yyz_xyyy, g_y_0_z_0_yyz_xyyyz, g_y_0_z_0_yyz_xyyz, g_y_0_z_0_yyz_xyyzz, g_y_0_z_0_yyz_xyzz, g_y_0_z_0_yyz_xyzzz, g_y_0_z_0_yyz_xzzz, g_y_0_z_0_yyz_xzzzz, g_y_0_z_0_yyz_yyyy, g_y_0_z_0_yyz_yyyyz, g_y_0_z_0_yyz_yyyz, g_y_0_z_0_yyz_yyyzz, g_y_0_z_0_yyz_yyzz, g_y_0_z_0_yyz_yyzzz, g_y_0_z_0_yyz_yzzz, g_y_0_z_0_yyz_yzzzz, g_y_0_z_0_yyz_zzzz, g_y_0_z_0_yyz_zzzzz, g_y_0_z_0_yyzz_xxxx, g_y_0_z_0_yyzz_xxxy, g_y_0_z_0_yyzz_xxxz, g_y_0_z_0_yyzz_xxyy, g_y_0_z_0_yyzz_xxyz, g_y_0_z_0_yyzz_xxzz, g_y_0_z_0_yyzz_xyyy, g_y_0_z_0_yyzz_xyyz, g_y_0_z_0_yyzz_xyzz, g_y_0_z_0_yyzz_xzzz, g_y_0_z_0_yyzz_yyyy, g_y_0_z_0_yyzz_yyyz, g_y_0_z_0_yyzz_yyzz, g_y_0_z_0_yyzz_yzzz, g_y_0_z_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyzz_xxxx[k] = -g_y_0_z_0_yyz_xxxx[k] * ab_z + g_y_0_z_0_yyz_xxxxz[k];

                g_y_0_z_0_yyzz_xxxy[k] = -g_y_0_z_0_yyz_xxxy[k] * ab_z + g_y_0_z_0_yyz_xxxyz[k];

                g_y_0_z_0_yyzz_xxxz[k] = -g_y_0_z_0_yyz_xxxz[k] * ab_z + g_y_0_z_0_yyz_xxxzz[k];

                g_y_0_z_0_yyzz_xxyy[k] = -g_y_0_z_0_yyz_xxyy[k] * ab_z + g_y_0_z_0_yyz_xxyyz[k];

                g_y_0_z_0_yyzz_xxyz[k] = -g_y_0_z_0_yyz_xxyz[k] * ab_z + g_y_0_z_0_yyz_xxyzz[k];

                g_y_0_z_0_yyzz_xxzz[k] = -g_y_0_z_0_yyz_xxzz[k] * ab_z + g_y_0_z_0_yyz_xxzzz[k];

                g_y_0_z_0_yyzz_xyyy[k] = -g_y_0_z_0_yyz_xyyy[k] * ab_z + g_y_0_z_0_yyz_xyyyz[k];

                g_y_0_z_0_yyzz_xyyz[k] = -g_y_0_z_0_yyz_xyyz[k] * ab_z + g_y_0_z_0_yyz_xyyzz[k];

                g_y_0_z_0_yyzz_xyzz[k] = -g_y_0_z_0_yyz_xyzz[k] * ab_z + g_y_0_z_0_yyz_xyzzz[k];

                g_y_0_z_0_yyzz_xzzz[k] = -g_y_0_z_0_yyz_xzzz[k] * ab_z + g_y_0_z_0_yyz_xzzzz[k];

                g_y_0_z_0_yyzz_yyyy[k] = -g_y_0_z_0_yyz_yyyy[k] * ab_z + g_y_0_z_0_yyz_yyyyz[k];

                g_y_0_z_0_yyzz_yyyz[k] = -g_y_0_z_0_yyz_yyyz[k] * ab_z + g_y_0_z_0_yyz_yyyzz[k];

                g_y_0_z_0_yyzz_yyzz[k] = -g_y_0_z_0_yyz_yyzz[k] * ab_z + g_y_0_z_0_yyz_yyzzz[k];

                g_y_0_z_0_yyzz_yzzz[k] = -g_y_0_z_0_yyz_yzzz[k] * ab_z + g_y_0_z_0_yyz_yzzzz[k];

                g_y_0_z_0_yyzz_zzzz[k] = -g_y_0_z_0_yyz_zzzz[k] * ab_z + g_y_0_z_0_yyz_zzzzz[k];
            }

            /// Set up 1320-1335 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1320 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1321 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1322 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1323 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1324 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1325 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1326 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1327 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1328 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1329 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1330 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1331 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1332 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1333 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1334 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yzz_xxxx, g_y_0_z_0_yzz_xxxxz, g_y_0_z_0_yzz_xxxy, g_y_0_z_0_yzz_xxxyz, g_y_0_z_0_yzz_xxxz, g_y_0_z_0_yzz_xxxzz, g_y_0_z_0_yzz_xxyy, g_y_0_z_0_yzz_xxyyz, g_y_0_z_0_yzz_xxyz, g_y_0_z_0_yzz_xxyzz, g_y_0_z_0_yzz_xxzz, g_y_0_z_0_yzz_xxzzz, g_y_0_z_0_yzz_xyyy, g_y_0_z_0_yzz_xyyyz, g_y_0_z_0_yzz_xyyz, g_y_0_z_0_yzz_xyyzz, g_y_0_z_0_yzz_xyzz, g_y_0_z_0_yzz_xyzzz, g_y_0_z_0_yzz_xzzz, g_y_0_z_0_yzz_xzzzz, g_y_0_z_0_yzz_yyyy, g_y_0_z_0_yzz_yyyyz, g_y_0_z_0_yzz_yyyz, g_y_0_z_0_yzz_yyyzz, g_y_0_z_0_yzz_yyzz, g_y_0_z_0_yzz_yyzzz, g_y_0_z_0_yzz_yzzz, g_y_0_z_0_yzz_yzzzz, g_y_0_z_0_yzz_zzzz, g_y_0_z_0_yzz_zzzzz, g_y_0_z_0_yzzz_xxxx, g_y_0_z_0_yzzz_xxxy, g_y_0_z_0_yzzz_xxxz, g_y_0_z_0_yzzz_xxyy, g_y_0_z_0_yzzz_xxyz, g_y_0_z_0_yzzz_xxzz, g_y_0_z_0_yzzz_xyyy, g_y_0_z_0_yzzz_xyyz, g_y_0_z_0_yzzz_xyzz, g_y_0_z_0_yzzz_xzzz, g_y_0_z_0_yzzz_yyyy, g_y_0_z_0_yzzz_yyyz, g_y_0_z_0_yzzz_yyzz, g_y_0_z_0_yzzz_yzzz, g_y_0_z_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yzzz_xxxx[k] = -g_y_0_z_0_yzz_xxxx[k] * ab_z + g_y_0_z_0_yzz_xxxxz[k];

                g_y_0_z_0_yzzz_xxxy[k] = -g_y_0_z_0_yzz_xxxy[k] * ab_z + g_y_0_z_0_yzz_xxxyz[k];

                g_y_0_z_0_yzzz_xxxz[k] = -g_y_0_z_0_yzz_xxxz[k] * ab_z + g_y_0_z_0_yzz_xxxzz[k];

                g_y_0_z_0_yzzz_xxyy[k] = -g_y_0_z_0_yzz_xxyy[k] * ab_z + g_y_0_z_0_yzz_xxyyz[k];

                g_y_0_z_0_yzzz_xxyz[k] = -g_y_0_z_0_yzz_xxyz[k] * ab_z + g_y_0_z_0_yzz_xxyzz[k];

                g_y_0_z_0_yzzz_xxzz[k] = -g_y_0_z_0_yzz_xxzz[k] * ab_z + g_y_0_z_0_yzz_xxzzz[k];

                g_y_0_z_0_yzzz_xyyy[k] = -g_y_0_z_0_yzz_xyyy[k] * ab_z + g_y_0_z_0_yzz_xyyyz[k];

                g_y_0_z_0_yzzz_xyyz[k] = -g_y_0_z_0_yzz_xyyz[k] * ab_z + g_y_0_z_0_yzz_xyyzz[k];

                g_y_0_z_0_yzzz_xyzz[k] = -g_y_0_z_0_yzz_xyzz[k] * ab_z + g_y_0_z_0_yzz_xyzzz[k];

                g_y_0_z_0_yzzz_xzzz[k] = -g_y_0_z_0_yzz_xzzz[k] * ab_z + g_y_0_z_0_yzz_xzzzz[k];

                g_y_0_z_0_yzzz_yyyy[k] = -g_y_0_z_0_yzz_yyyy[k] * ab_z + g_y_0_z_0_yzz_yyyyz[k];

                g_y_0_z_0_yzzz_yyyz[k] = -g_y_0_z_0_yzz_yyyz[k] * ab_z + g_y_0_z_0_yzz_yyyzz[k];

                g_y_0_z_0_yzzz_yyzz[k] = -g_y_0_z_0_yzz_yyzz[k] * ab_z + g_y_0_z_0_yzz_yyzzz[k];

                g_y_0_z_0_yzzz_yzzz[k] = -g_y_0_z_0_yzz_yzzz[k] * ab_z + g_y_0_z_0_yzz_yzzzz[k];

                g_y_0_z_0_yzzz_zzzz[k] = -g_y_0_z_0_yzz_zzzz[k] * ab_z + g_y_0_z_0_yzz_zzzzz[k];
            }

            /// Set up 1335-1350 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_zzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1335 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1336 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1337 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1338 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1339 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1340 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1341 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1342 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1343 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1344 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1345 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1346 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1347 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1348 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1349 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_zzz_xxxx, g_y_0_z_0_zzz_xxxxz, g_y_0_z_0_zzz_xxxy, g_y_0_z_0_zzz_xxxyz, g_y_0_z_0_zzz_xxxz, g_y_0_z_0_zzz_xxxzz, g_y_0_z_0_zzz_xxyy, g_y_0_z_0_zzz_xxyyz, g_y_0_z_0_zzz_xxyz, g_y_0_z_0_zzz_xxyzz, g_y_0_z_0_zzz_xxzz, g_y_0_z_0_zzz_xxzzz, g_y_0_z_0_zzz_xyyy, g_y_0_z_0_zzz_xyyyz, g_y_0_z_0_zzz_xyyz, g_y_0_z_0_zzz_xyyzz, g_y_0_z_0_zzz_xyzz, g_y_0_z_0_zzz_xyzzz, g_y_0_z_0_zzz_xzzz, g_y_0_z_0_zzz_xzzzz, g_y_0_z_0_zzz_yyyy, g_y_0_z_0_zzz_yyyyz, g_y_0_z_0_zzz_yyyz, g_y_0_z_0_zzz_yyyzz, g_y_0_z_0_zzz_yyzz, g_y_0_z_0_zzz_yyzzz, g_y_0_z_0_zzz_yzzz, g_y_0_z_0_zzz_yzzzz, g_y_0_z_0_zzz_zzzz, g_y_0_z_0_zzz_zzzzz, g_y_0_z_0_zzzz_xxxx, g_y_0_z_0_zzzz_xxxy, g_y_0_z_0_zzzz_xxxz, g_y_0_z_0_zzzz_xxyy, g_y_0_z_0_zzzz_xxyz, g_y_0_z_0_zzzz_xxzz, g_y_0_z_0_zzzz_xyyy, g_y_0_z_0_zzzz_xyyz, g_y_0_z_0_zzzz_xyzz, g_y_0_z_0_zzzz_xzzz, g_y_0_z_0_zzzz_yyyy, g_y_0_z_0_zzzz_yyyz, g_y_0_z_0_zzzz_yyzz, g_y_0_z_0_zzzz_yzzz, g_y_0_z_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_zzzz_xxxx[k] = -g_y_0_z_0_zzz_xxxx[k] * ab_z + g_y_0_z_0_zzz_xxxxz[k];

                g_y_0_z_0_zzzz_xxxy[k] = -g_y_0_z_0_zzz_xxxy[k] * ab_z + g_y_0_z_0_zzz_xxxyz[k];

                g_y_0_z_0_zzzz_xxxz[k] = -g_y_0_z_0_zzz_xxxz[k] * ab_z + g_y_0_z_0_zzz_xxxzz[k];

                g_y_0_z_0_zzzz_xxyy[k] = -g_y_0_z_0_zzz_xxyy[k] * ab_z + g_y_0_z_0_zzz_xxyyz[k];

                g_y_0_z_0_zzzz_xxyz[k] = -g_y_0_z_0_zzz_xxyz[k] * ab_z + g_y_0_z_0_zzz_xxyzz[k];

                g_y_0_z_0_zzzz_xxzz[k] = -g_y_0_z_0_zzz_xxzz[k] * ab_z + g_y_0_z_0_zzz_xxzzz[k];

                g_y_0_z_0_zzzz_xyyy[k] = -g_y_0_z_0_zzz_xyyy[k] * ab_z + g_y_0_z_0_zzz_xyyyz[k];

                g_y_0_z_0_zzzz_xyyz[k] = -g_y_0_z_0_zzz_xyyz[k] * ab_z + g_y_0_z_0_zzz_xyyzz[k];

                g_y_0_z_0_zzzz_xyzz[k] = -g_y_0_z_0_zzz_xyzz[k] * ab_z + g_y_0_z_0_zzz_xyzzz[k];

                g_y_0_z_0_zzzz_xzzz[k] = -g_y_0_z_0_zzz_xzzz[k] * ab_z + g_y_0_z_0_zzz_xzzzz[k];

                g_y_0_z_0_zzzz_yyyy[k] = -g_y_0_z_0_zzz_yyyy[k] * ab_z + g_y_0_z_0_zzz_yyyyz[k];

                g_y_0_z_0_zzzz_yyyz[k] = -g_y_0_z_0_zzz_yyyz[k] * ab_z + g_y_0_z_0_zzz_yyyzz[k];

                g_y_0_z_0_zzzz_yyzz[k] = -g_y_0_z_0_zzz_yyzz[k] * ab_z + g_y_0_z_0_zzz_yyzzz[k];

                g_y_0_z_0_zzzz_yzzz[k] = -g_y_0_z_0_zzz_yzzz[k] * ab_z + g_y_0_z_0_zzz_yzzzz[k];

                g_y_0_z_0_zzzz_zzzz[k] = -g_y_0_z_0_zzz_zzzz[k] * ab_z + g_y_0_z_0_zzz_zzzzz[k];
            }

            /// Set up 1350-1365 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxx_xxxx = cbuffer.data(gg_geom_1010_off + 1350 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_xxxy = cbuffer.data(gg_geom_1010_off + 1351 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_xxxz = cbuffer.data(gg_geom_1010_off + 1352 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_xxyy = cbuffer.data(gg_geom_1010_off + 1353 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_xxyz = cbuffer.data(gg_geom_1010_off + 1354 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_xxzz = cbuffer.data(gg_geom_1010_off + 1355 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_xyyy = cbuffer.data(gg_geom_1010_off + 1356 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_xyyz = cbuffer.data(gg_geom_1010_off + 1357 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_xyzz = cbuffer.data(gg_geom_1010_off + 1358 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_xzzz = cbuffer.data(gg_geom_1010_off + 1359 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_yyyy = cbuffer.data(gg_geom_1010_off + 1360 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_yyyz = cbuffer.data(gg_geom_1010_off + 1361 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_yyzz = cbuffer.data(gg_geom_1010_off + 1362 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_yzzz = cbuffer.data(gg_geom_1010_off + 1363 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_zzzz = cbuffer.data(gg_geom_1010_off + 1364 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxx_xxxx, g_z_0_x_0_xxx_xxxxx, g_z_0_x_0_xxx_xxxxy, g_z_0_x_0_xxx_xxxxz, g_z_0_x_0_xxx_xxxy, g_z_0_x_0_xxx_xxxyy, g_z_0_x_0_xxx_xxxyz, g_z_0_x_0_xxx_xxxz, g_z_0_x_0_xxx_xxxzz, g_z_0_x_0_xxx_xxyy, g_z_0_x_0_xxx_xxyyy, g_z_0_x_0_xxx_xxyyz, g_z_0_x_0_xxx_xxyz, g_z_0_x_0_xxx_xxyzz, g_z_0_x_0_xxx_xxzz, g_z_0_x_0_xxx_xxzzz, g_z_0_x_0_xxx_xyyy, g_z_0_x_0_xxx_xyyyy, g_z_0_x_0_xxx_xyyyz, g_z_0_x_0_xxx_xyyz, g_z_0_x_0_xxx_xyyzz, g_z_0_x_0_xxx_xyzz, g_z_0_x_0_xxx_xyzzz, g_z_0_x_0_xxx_xzzz, g_z_0_x_0_xxx_xzzzz, g_z_0_x_0_xxx_yyyy, g_z_0_x_0_xxx_yyyz, g_z_0_x_0_xxx_yyzz, g_z_0_x_0_xxx_yzzz, g_z_0_x_0_xxx_zzzz, g_z_0_x_0_xxxx_xxxx, g_z_0_x_0_xxxx_xxxy, g_z_0_x_0_xxxx_xxxz, g_z_0_x_0_xxxx_xxyy, g_z_0_x_0_xxxx_xxyz, g_z_0_x_0_xxxx_xxzz, g_z_0_x_0_xxxx_xyyy, g_z_0_x_0_xxxx_xyyz, g_z_0_x_0_xxxx_xyzz, g_z_0_x_0_xxxx_xzzz, g_z_0_x_0_xxxx_yyyy, g_z_0_x_0_xxxx_yyyz, g_z_0_x_0_xxxx_yyzz, g_z_0_x_0_xxxx_yzzz, g_z_0_x_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxx_xxxx[k] = -g_z_0_x_0_xxx_xxxx[k] * ab_x + g_z_0_x_0_xxx_xxxxx[k];

                g_z_0_x_0_xxxx_xxxy[k] = -g_z_0_x_0_xxx_xxxy[k] * ab_x + g_z_0_x_0_xxx_xxxxy[k];

                g_z_0_x_0_xxxx_xxxz[k] = -g_z_0_x_0_xxx_xxxz[k] * ab_x + g_z_0_x_0_xxx_xxxxz[k];

                g_z_0_x_0_xxxx_xxyy[k] = -g_z_0_x_0_xxx_xxyy[k] * ab_x + g_z_0_x_0_xxx_xxxyy[k];

                g_z_0_x_0_xxxx_xxyz[k] = -g_z_0_x_0_xxx_xxyz[k] * ab_x + g_z_0_x_0_xxx_xxxyz[k];

                g_z_0_x_0_xxxx_xxzz[k] = -g_z_0_x_0_xxx_xxzz[k] * ab_x + g_z_0_x_0_xxx_xxxzz[k];

                g_z_0_x_0_xxxx_xyyy[k] = -g_z_0_x_0_xxx_xyyy[k] * ab_x + g_z_0_x_0_xxx_xxyyy[k];

                g_z_0_x_0_xxxx_xyyz[k] = -g_z_0_x_0_xxx_xyyz[k] * ab_x + g_z_0_x_0_xxx_xxyyz[k];

                g_z_0_x_0_xxxx_xyzz[k] = -g_z_0_x_0_xxx_xyzz[k] * ab_x + g_z_0_x_0_xxx_xxyzz[k];

                g_z_0_x_0_xxxx_xzzz[k] = -g_z_0_x_0_xxx_xzzz[k] * ab_x + g_z_0_x_0_xxx_xxzzz[k];

                g_z_0_x_0_xxxx_yyyy[k] = -g_z_0_x_0_xxx_yyyy[k] * ab_x + g_z_0_x_0_xxx_xyyyy[k];

                g_z_0_x_0_xxxx_yyyz[k] = -g_z_0_x_0_xxx_yyyz[k] * ab_x + g_z_0_x_0_xxx_xyyyz[k];

                g_z_0_x_0_xxxx_yyzz[k] = -g_z_0_x_0_xxx_yyzz[k] * ab_x + g_z_0_x_0_xxx_xyyzz[k];

                g_z_0_x_0_xxxx_yzzz[k] = -g_z_0_x_0_xxx_yzzz[k] * ab_x + g_z_0_x_0_xxx_xyzzz[k];

                g_z_0_x_0_xxxx_zzzz[k] = -g_z_0_x_0_xxx_zzzz[k] * ab_x + g_z_0_x_0_xxx_xzzzz[k];
            }

            /// Set up 1365-1380 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxy_xxxx = cbuffer.data(gg_geom_1010_off + 1365 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_xxxy = cbuffer.data(gg_geom_1010_off + 1366 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_xxxz = cbuffer.data(gg_geom_1010_off + 1367 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_xxyy = cbuffer.data(gg_geom_1010_off + 1368 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_xxyz = cbuffer.data(gg_geom_1010_off + 1369 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_xxzz = cbuffer.data(gg_geom_1010_off + 1370 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_xyyy = cbuffer.data(gg_geom_1010_off + 1371 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_xyyz = cbuffer.data(gg_geom_1010_off + 1372 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_xyzz = cbuffer.data(gg_geom_1010_off + 1373 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_xzzz = cbuffer.data(gg_geom_1010_off + 1374 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_yyyy = cbuffer.data(gg_geom_1010_off + 1375 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_yyyz = cbuffer.data(gg_geom_1010_off + 1376 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_yyzz = cbuffer.data(gg_geom_1010_off + 1377 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_yzzz = cbuffer.data(gg_geom_1010_off + 1378 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_zzzz = cbuffer.data(gg_geom_1010_off + 1379 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxy_xxxx, g_z_0_x_0_xxxy_xxxy, g_z_0_x_0_xxxy_xxxz, g_z_0_x_0_xxxy_xxyy, g_z_0_x_0_xxxy_xxyz, g_z_0_x_0_xxxy_xxzz, g_z_0_x_0_xxxy_xyyy, g_z_0_x_0_xxxy_xyyz, g_z_0_x_0_xxxy_xyzz, g_z_0_x_0_xxxy_xzzz, g_z_0_x_0_xxxy_yyyy, g_z_0_x_0_xxxy_yyyz, g_z_0_x_0_xxxy_yyzz, g_z_0_x_0_xxxy_yzzz, g_z_0_x_0_xxxy_zzzz, g_z_0_x_0_xxy_xxxx, g_z_0_x_0_xxy_xxxxx, g_z_0_x_0_xxy_xxxxy, g_z_0_x_0_xxy_xxxxz, g_z_0_x_0_xxy_xxxy, g_z_0_x_0_xxy_xxxyy, g_z_0_x_0_xxy_xxxyz, g_z_0_x_0_xxy_xxxz, g_z_0_x_0_xxy_xxxzz, g_z_0_x_0_xxy_xxyy, g_z_0_x_0_xxy_xxyyy, g_z_0_x_0_xxy_xxyyz, g_z_0_x_0_xxy_xxyz, g_z_0_x_0_xxy_xxyzz, g_z_0_x_0_xxy_xxzz, g_z_0_x_0_xxy_xxzzz, g_z_0_x_0_xxy_xyyy, g_z_0_x_0_xxy_xyyyy, g_z_0_x_0_xxy_xyyyz, g_z_0_x_0_xxy_xyyz, g_z_0_x_0_xxy_xyyzz, g_z_0_x_0_xxy_xyzz, g_z_0_x_0_xxy_xyzzz, g_z_0_x_0_xxy_xzzz, g_z_0_x_0_xxy_xzzzz, g_z_0_x_0_xxy_yyyy, g_z_0_x_0_xxy_yyyz, g_z_0_x_0_xxy_yyzz, g_z_0_x_0_xxy_yzzz, g_z_0_x_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxy_xxxx[k] = -g_z_0_x_0_xxy_xxxx[k] * ab_x + g_z_0_x_0_xxy_xxxxx[k];

                g_z_0_x_0_xxxy_xxxy[k] = -g_z_0_x_0_xxy_xxxy[k] * ab_x + g_z_0_x_0_xxy_xxxxy[k];

                g_z_0_x_0_xxxy_xxxz[k] = -g_z_0_x_0_xxy_xxxz[k] * ab_x + g_z_0_x_0_xxy_xxxxz[k];

                g_z_0_x_0_xxxy_xxyy[k] = -g_z_0_x_0_xxy_xxyy[k] * ab_x + g_z_0_x_0_xxy_xxxyy[k];

                g_z_0_x_0_xxxy_xxyz[k] = -g_z_0_x_0_xxy_xxyz[k] * ab_x + g_z_0_x_0_xxy_xxxyz[k];

                g_z_0_x_0_xxxy_xxzz[k] = -g_z_0_x_0_xxy_xxzz[k] * ab_x + g_z_0_x_0_xxy_xxxzz[k];

                g_z_0_x_0_xxxy_xyyy[k] = -g_z_0_x_0_xxy_xyyy[k] * ab_x + g_z_0_x_0_xxy_xxyyy[k];

                g_z_0_x_0_xxxy_xyyz[k] = -g_z_0_x_0_xxy_xyyz[k] * ab_x + g_z_0_x_0_xxy_xxyyz[k];

                g_z_0_x_0_xxxy_xyzz[k] = -g_z_0_x_0_xxy_xyzz[k] * ab_x + g_z_0_x_0_xxy_xxyzz[k];

                g_z_0_x_0_xxxy_xzzz[k] = -g_z_0_x_0_xxy_xzzz[k] * ab_x + g_z_0_x_0_xxy_xxzzz[k];

                g_z_0_x_0_xxxy_yyyy[k] = -g_z_0_x_0_xxy_yyyy[k] * ab_x + g_z_0_x_0_xxy_xyyyy[k];

                g_z_0_x_0_xxxy_yyyz[k] = -g_z_0_x_0_xxy_yyyz[k] * ab_x + g_z_0_x_0_xxy_xyyyz[k];

                g_z_0_x_0_xxxy_yyzz[k] = -g_z_0_x_0_xxy_yyzz[k] * ab_x + g_z_0_x_0_xxy_xyyzz[k];

                g_z_0_x_0_xxxy_yzzz[k] = -g_z_0_x_0_xxy_yzzz[k] * ab_x + g_z_0_x_0_xxy_xyzzz[k];

                g_z_0_x_0_xxxy_zzzz[k] = -g_z_0_x_0_xxy_zzzz[k] * ab_x + g_z_0_x_0_xxy_xzzzz[k];
            }

            /// Set up 1380-1395 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxz_xxxx = cbuffer.data(gg_geom_1010_off + 1380 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_xxxy = cbuffer.data(gg_geom_1010_off + 1381 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_xxxz = cbuffer.data(gg_geom_1010_off + 1382 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_xxyy = cbuffer.data(gg_geom_1010_off + 1383 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_xxyz = cbuffer.data(gg_geom_1010_off + 1384 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_xxzz = cbuffer.data(gg_geom_1010_off + 1385 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_xyyy = cbuffer.data(gg_geom_1010_off + 1386 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_xyyz = cbuffer.data(gg_geom_1010_off + 1387 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_xyzz = cbuffer.data(gg_geom_1010_off + 1388 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_xzzz = cbuffer.data(gg_geom_1010_off + 1389 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_yyyy = cbuffer.data(gg_geom_1010_off + 1390 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_yyyz = cbuffer.data(gg_geom_1010_off + 1391 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_yyzz = cbuffer.data(gg_geom_1010_off + 1392 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_yzzz = cbuffer.data(gg_geom_1010_off + 1393 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_zzzz = cbuffer.data(gg_geom_1010_off + 1394 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxz_xxxx, g_z_0_x_0_xxxz_xxxy, g_z_0_x_0_xxxz_xxxz, g_z_0_x_0_xxxz_xxyy, g_z_0_x_0_xxxz_xxyz, g_z_0_x_0_xxxz_xxzz, g_z_0_x_0_xxxz_xyyy, g_z_0_x_0_xxxz_xyyz, g_z_0_x_0_xxxz_xyzz, g_z_0_x_0_xxxz_xzzz, g_z_0_x_0_xxxz_yyyy, g_z_0_x_0_xxxz_yyyz, g_z_0_x_0_xxxz_yyzz, g_z_0_x_0_xxxz_yzzz, g_z_0_x_0_xxxz_zzzz, g_z_0_x_0_xxz_xxxx, g_z_0_x_0_xxz_xxxxx, g_z_0_x_0_xxz_xxxxy, g_z_0_x_0_xxz_xxxxz, g_z_0_x_0_xxz_xxxy, g_z_0_x_0_xxz_xxxyy, g_z_0_x_0_xxz_xxxyz, g_z_0_x_0_xxz_xxxz, g_z_0_x_0_xxz_xxxzz, g_z_0_x_0_xxz_xxyy, g_z_0_x_0_xxz_xxyyy, g_z_0_x_0_xxz_xxyyz, g_z_0_x_0_xxz_xxyz, g_z_0_x_0_xxz_xxyzz, g_z_0_x_0_xxz_xxzz, g_z_0_x_0_xxz_xxzzz, g_z_0_x_0_xxz_xyyy, g_z_0_x_0_xxz_xyyyy, g_z_0_x_0_xxz_xyyyz, g_z_0_x_0_xxz_xyyz, g_z_0_x_0_xxz_xyyzz, g_z_0_x_0_xxz_xyzz, g_z_0_x_0_xxz_xyzzz, g_z_0_x_0_xxz_xzzz, g_z_0_x_0_xxz_xzzzz, g_z_0_x_0_xxz_yyyy, g_z_0_x_0_xxz_yyyz, g_z_0_x_0_xxz_yyzz, g_z_0_x_0_xxz_yzzz, g_z_0_x_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxz_xxxx[k] = -g_z_0_x_0_xxz_xxxx[k] * ab_x + g_z_0_x_0_xxz_xxxxx[k];

                g_z_0_x_0_xxxz_xxxy[k] = -g_z_0_x_0_xxz_xxxy[k] * ab_x + g_z_0_x_0_xxz_xxxxy[k];

                g_z_0_x_0_xxxz_xxxz[k] = -g_z_0_x_0_xxz_xxxz[k] * ab_x + g_z_0_x_0_xxz_xxxxz[k];

                g_z_0_x_0_xxxz_xxyy[k] = -g_z_0_x_0_xxz_xxyy[k] * ab_x + g_z_0_x_0_xxz_xxxyy[k];

                g_z_0_x_0_xxxz_xxyz[k] = -g_z_0_x_0_xxz_xxyz[k] * ab_x + g_z_0_x_0_xxz_xxxyz[k];

                g_z_0_x_0_xxxz_xxzz[k] = -g_z_0_x_0_xxz_xxzz[k] * ab_x + g_z_0_x_0_xxz_xxxzz[k];

                g_z_0_x_0_xxxz_xyyy[k] = -g_z_0_x_0_xxz_xyyy[k] * ab_x + g_z_0_x_0_xxz_xxyyy[k];

                g_z_0_x_0_xxxz_xyyz[k] = -g_z_0_x_0_xxz_xyyz[k] * ab_x + g_z_0_x_0_xxz_xxyyz[k];

                g_z_0_x_0_xxxz_xyzz[k] = -g_z_0_x_0_xxz_xyzz[k] * ab_x + g_z_0_x_0_xxz_xxyzz[k];

                g_z_0_x_0_xxxz_xzzz[k] = -g_z_0_x_0_xxz_xzzz[k] * ab_x + g_z_0_x_0_xxz_xxzzz[k];

                g_z_0_x_0_xxxz_yyyy[k] = -g_z_0_x_0_xxz_yyyy[k] * ab_x + g_z_0_x_0_xxz_xyyyy[k];

                g_z_0_x_0_xxxz_yyyz[k] = -g_z_0_x_0_xxz_yyyz[k] * ab_x + g_z_0_x_0_xxz_xyyyz[k];

                g_z_0_x_0_xxxz_yyzz[k] = -g_z_0_x_0_xxz_yyzz[k] * ab_x + g_z_0_x_0_xxz_xyyzz[k];

                g_z_0_x_0_xxxz_yzzz[k] = -g_z_0_x_0_xxz_yzzz[k] * ab_x + g_z_0_x_0_xxz_xyzzz[k];

                g_z_0_x_0_xxxz_zzzz[k] = -g_z_0_x_0_xxz_zzzz[k] * ab_x + g_z_0_x_0_xxz_xzzzz[k];
            }

            /// Set up 1395-1410 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyy_xxxx = cbuffer.data(gg_geom_1010_off + 1395 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_xxxy = cbuffer.data(gg_geom_1010_off + 1396 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_xxxz = cbuffer.data(gg_geom_1010_off + 1397 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_xxyy = cbuffer.data(gg_geom_1010_off + 1398 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_xxyz = cbuffer.data(gg_geom_1010_off + 1399 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_xxzz = cbuffer.data(gg_geom_1010_off + 1400 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_xyyy = cbuffer.data(gg_geom_1010_off + 1401 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_xyyz = cbuffer.data(gg_geom_1010_off + 1402 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_xyzz = cbuffer.data(gg_geom_1010_off + 1403 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_xzzz = cbuffer.data(gg_geom_1010_off + 1404 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_yyyy = cbuffer.data(gg_geom_1010_off + 1405 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_yyyz = cbuffer.data(gg_geom_1010_off + 1406 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_yyzz = cbuffer.data(gg_geom_1010_off + 1407 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_yzzz = cbuffer.data(gg_geom_1010_off + 1408 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_zzzz = cbuffer.data(gg_geom_1010_off + 1409 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyy_xxxx, g_z_0_x_0_xxyy_xxxy, g_z_0_x_0_xxyy_xxxz, g_z_0_x_0_xxyy_xxyy, g_z_0_x_0_xxyy_xxyz, g_z_0_x_0_xxyy_xxzz, g_z_0_x_0_xxyy_xyyy, g_z_0_x_0_xxyy_xyyz, g_z_0_x_0_xxyy_xyzz, g_z_0_x_0_xxyy_xzzz, g_z_0_x_0_xxyy_yyyy, g_z_0_x_0_xxyy_yyyz, g_z_0_x_0_xxyy_yyzz, g_z_0_x_0_xxyy_yzzz, g_z_0_x_0_xxyy_zzzz, g_z_0_x_0_xyy_xxxx, g_z_0_x_0_xyy_xxxxx, g_z_0_x_0_xyy_xxxxy, g_z_0_x_0_xyy_xxxxz, g_z_0_x_0_xyy_xxxy, g_z_0_x_0_xyy_xxxyy, g_z_0_x_0_xyy_xxxyz, g_z_0_x_0_xyy_xxxz, g_z_0_x_0_xyy_xxxzz, g_z_0_x_0_xyy_xxyy, g_z_0_x_0_xyy_xxyyy, g_z_0_x_0_xyy_xxyyz, g_z_0_x_0_xyy_xxyz, g_z_0_x_0_xyy_xxyzz, g_z_0_x_0_xyy_xxzz, g_z_0_x_0_xyy_xxzzz, g_z_0_x_0_xyy_xyyy, g_z_0_x_0_xyy_xyyyy, g_z_0_x_0_xyy_xyyyz, g_z_0_x_0_xyy_xyyz, g_z_0_x_0_xyy_xyyzz, g_z_0_x_0_xyy_xyzz, g_z_0_x_0_xyy_xyzzz, g_z_0_x_0_xyy_xzzz, g_z_0_x_0_xyy_xzzzz, g_z_0_x_0_xyy_yyyy, g_z_0_x_0_xyy_yyyz, g_z_0_x_0_xyy_yyzz, g_z_0_x_0_xyy_yzzz, g_z_0_x_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyy_xxxx[k] = -g_z_0_x_0_xyy_xxxx[k] * ab_x + g_z_0_x_0_xyy_xxxxx[k];

                g_z_0_x_0_xxyy_xxxy[k] = -g_z_0_x_0_xyy_xxxy[k] * ab_x + g_z_0_x_0_xyy_xxxxy[k];

                g_z_0_x_0_xxyy_xxxz[k] = -g_z_0_x_0_xyy_xxxz[k] * ab_x + g_z_0_x_0_xyy_xxxxz[k];

                g_z_0_x_0_xxyy_xxyy[k] = -g_z_0_x_0_xyy_xxyy[k] * ab_x + g_z_0_x_0_xyy_xxxyy[k];

                g_z_0_x_0_xxyy_xxyz[k] = -g_z_0_x_0_xyy_xxyz[k] * ab_x + g_z_0_x_0_xyy_xxxyz[k];

                g_z_0_x_0_xxyy_xxzz[k] = -g_z_0_x_0_xyy_xxzz[k] * ab_x + g_z_0_x_0_xyy_xxxzz[k];

                g_z_0_x_0_xxyy_xyyy[k] = -g_z_0_x_0_xyy_xyyy[k] * ab_x + g_z_0_x_0_xyy_xxyyy[k];

                g_z_0_x_0_xxyy_xyyz[k] = -g_z_0_x_0_xyy_xyyz[k] * ab_x + g_z_0_x_0_xyy_xxyyz[k];

                g_z_0_x_0_xxyy_xyzz[k] = -g_z_0_x_0_xyy_xyzz[k] * ab_x + g_z_0_x_0_xyy_xxyzz[k];

                g_z_0_x_0_xxyy_xzzz[k] = -g_z_0_x_0_xyy_xzzz[k] * ab_x + g_z_0_x_0_xyy_xxzzz[k];

                g_z_0_x_0_xxyy_yyyy[k] = -g_z_0_x_0_xyy_yyyy[k] * ab_x + g_z_0_x_0_xyy_xyyyy[k];

                g_z_0_x_0_xxyy_yyyz[k] = -g_z_0_x_0_xyy_yyyz[k] * ab_x + g_z_0_x_0_xyy_xyyyz[k];

                g_z_0_x_0_xxyy_yyzz[k] = -g_z_0_x_0_xyy_yyzz[k] * ab_x + g_z_0_x_0_xyy_xyyzz[k];

                g_z_0_x_0_xxyy_yzzz[k] = -g_z_0_x_0_xyy_yzzz[k] * ab_x + g_z_0_x_0_xyy_xyzzz[k];

                g_z_0_x_0_xxyy_zzzz[k] = -g_z_0_x_0_xyy_zzzz[k] * ab_x + g_z_0_x_0_xyy_xzzzz[k];
            }

            /// Set up 1410-1425 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyz_xxxx = cbuffer.data(gg_geom_1010_off + 1410 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_xxxy = cbuffer.data(gg_geom_1010_off + 1411 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_xxxz = cbuffer.data(gg_geom_1010_off + 1412 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_xxyy = cbuffer.data(gg_geom_1010_off + 1413 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_xxyz = cbuffer.data(gg_geom_1010_off + 1414 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_xxzz = cbuffer.data(gg_geom_1010_off + 1415 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_xyyy = cbuffer.data(gg_geom_1010_off + 1416 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_xyyz = cbuffer.data(gg_geom_1010_off + 1417 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_xyzz = cbuffer.data(gg_geom_1010_off + 1418 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_xzzz = cbuffer.data(gg_geom_1010_off + 1419 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_yyyy = cbuffer.data(gg_geom_1010_off + 1420 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_yyyz = cbuffer.data(gg_geom_1010_off + 1421 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_yyzz = cbuffer.data(gg_geom_1010_off + 1422 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_yzzz = cbuffer.data(gg_geom_1010_off + 1423 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_zzzz = cbuffer.data(gg_geom_1010_off + 1424 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyz_xxxx, g_z_0_x_0_xxyz_xxxy, g_z_0_x_0_xxyz_xxxz, g_z_0_x_0_xxyz_xxyy, g_z_0_x_0_xxyz_xxyz, g_z_0_x_0_xxyz_xxzz, g_z_0_x_0_xxyz_xyyy, g_z_0_x_0_xxyz_xyyz, g_z_0_x_0_xxyz_xyzz, g_z_0_x_0_xxyz_xzzz, g_z_0_x_0_xxyz_yyyy, g_z_0_x_0_xxyz_yyyz, g_z_0_x_0_xxyz_yyzz, g_z_0_x_0_xxyz_yzzz, g_z_0_x_0_xxyz_zzzz, g_z_0_x_0_xyz_xxxx, g_z_0_x_0_xyz_xxxxx, g_z_0_x_0_xyz_xxxxy, g_z_0_x_0_xyz_xxxxz, g_z_0_x_0_xyz_xxxy, g_z_0_x_0_xyz_xxxyy, g_z_0_x_0_xyz_xxxyz, g_z_0_x_0_xyz_xxxz, g_z_0_x_0_xyz_xxxzz, g_z_0_x_0_xyz_xxyy, g_z_0_x_0_xyz_xxyyy, g_z_0_x_0_xyz_xxyyz, g_z_0_x_0_xyz_xxyz, g_z_0_x_0_xyz_xxyzz, g_z_0_x_0_xyz_xxzz, g_z_0_x_0_xyz_xxzzz, g_z_0_x_0_xyz_xyyy, g_z_0_x_0_xyz_xyyyy, g_z_0_x_0_xyz_xyyyz, g_z_0_x_0_xyz_xyyz, g_z_0_x_0_xyz_xyyzz, g_z_0_x_0_xyz_xyzz, g_z_0_x_0_xyz_xyzzz, g_z_0_x_0_xyz_xzzz, g_z_0_x_0_xyz_xzzzz, g_z_0_x_0_xyz_yyyy, g_z_0_x_0_xyz_yyyz, g_z_0_x_0_xyz_yyzz, g_z_0_x_0_xyz_yzzz, g_z_0_x_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyz_xxxx[k] = -g_z_0_x_0_xyz_xxxx[k] * ab_x + g_z_0_x_0_xyz_xxxxx[k];

                g_z_0_x_0_xxyz_xxxy[k] = -g_z_0_x_0_xyz_xxxy[k] * ab_x + g_z_0_x_0_xyz_xxxxy[k];

                g_z_0_x_0_xxyz_xxxz[k] = -g_z_0_x_0_xyz_xxxz[k] * ab_x + g_z_0_x_0_xyz_xxxxz[k];

                g_z_0_x_0_xxyz_xxyy[k] = -g_z_0_x_0_xyz_xxyy[k] * ab_x + g_z_0_x_0_xyz_xxxyy[k];

                g_z_0_x_0_xxyz_xxyz[k] = -g_z_0_x_0_xyz_xxyz[k] * ab_x + g_z_0_x_0_xyz_xxxyz[k];

                g_z_0_x_0_xxyz_xxzz[k] = -g_z_0_x_0_xyz_xxzz[k] * ab_x + g_z_0_x_0_xyz_xxxzz[k];

                g_z_0_x_0_xxyz_xyyy[k] = -g_z_0_x_0_xyz_xyyy[k] * ab_x + g_z_0_x_0_xyz_xxyyy[k];

                g_z_0_x_0_xxyz_xyyz[k] = -g_z_0_x_0_xyz_xyyz[k] * ab_x + g_z_0_x_0_xyz_xxyyz[k];

                g_z_0_x_0_xxyz_xyzz[k] = -g_z_0_x_0_xyz_xyzz[k] * ab_x + g_z_0_x_0_xyz_xxyzz[k];

                g_z_0_x_0_xxyz_xzzz[k] = -g_z_0_x_0_xyz_xzzz[k] * ab_x + g_z_0_x_0_xyz_xxzzz[k];

                g_z_0_x_0_xxyz_yyyy[k] = -g_z_0_x_0_xyz_yyyy[k] * ab_x + g_z_0_x_0_xyz_xyyyy[k];

                g_z_0_x_0_xxyz_yyyz[k] = -g_z_0_x_0_xyz_yyyz[k] * ab_x + g_z_0_x_0_xyz_xyyyz[k];

                g_z_0_x_0_xxyz_yyzz[k] = -g_z_0_x_0_xyz_yyzz[k] * ab_x + g_z_0_x_0_xyz_xyyzz[k];

                g_z_0_x_0_xxyz_yzzz[k] = -g_z_0_x_0_xyz_yzzz[k] * ab_x + g_z_0_x_0_xyz_xyzzz[k];

                g_z_0_x_0_xxyz_zzzz[k] = -g_z_0_x_0_xyz_zzzz[k] * ab_x + g_z_0_x_0_xyz_xzzzz[k];
            }

            /// Set up 1425-1440 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxzz_xxxx = cbuffer.data(gg_geom_1010_off + 1425 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_xxxy = cbuffer.data(gg_geom_1010_off + 1426 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_xxxz = cbuffer.data(gg_geom_1010_off + 1427 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_xxyy = cbuffer.data(gg_geom_1010_off + 1428 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_xxyz = cbuffer.data(gg_geom_1010_off + 1429 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_xxzz = cbuffer.data(gg_geom_1010_off + 1430 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_xyyy = cbuffer.data(gg_geom_1010_off + 1431 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_xyyz = cbuffer.data(gg_geom_1010_off + 1432 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_xyzz = cbuffer.data(gg_geom_1010_off + 1433 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_xzzz = cbuffer.data(gg_geom_1010_off + 1434 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_yyyy = cbuffer.data(gg_geom_1010_off + 1435 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_yyyz = cbuffer.data(gg_geom_1010_off + 1436 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_yyzz = cbuffer.data(gg_geom_1010_off + 1437 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_yzzz = cbuffer.data(gg_geom_1010_off + 1438 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_zzzz = cbuffer.data(gg_geom_1010_off + 1439 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxzz_xxxx, g_z_0_x_0_xxzz_xxxy, g_z_0_x_0_xxzz_xxxz, g_z_0_x_0_xxzz_xxyy, g_z_0_x_0_xxzz_xxyz, g_z_0_x_0_xxzz_xxzz, g_z_0_x_0_xxzz_xyyy, g_z_0_x_0_xxzz_xyyz, g_z_0_x_0_xxzz_xyzz, g_z_0_x_0_xxzz_xzzz, g_z_0_x_0_xxzz_yyyy, g_z_0_x_0_xxzz_yyyz, g_z_0_x_0_xxzz_yyzz, g_z_0_x_0_xxzz_yzzz, g_z_0_x_0_xxzz_zzzz, g_z_0_x_0_xzz_xxxx, g_z_0_x_0_xzz_xxxxx, g_z_0_x_0_xzz_xxxxy, g_z_0_x_0_xzz_xxxxz, g_z_0_x_0_xzz_xxxy, g_z_0_x_0_xzz_xxxyy, g_z_0_x_0_xzz_xxxyz, g_z_0_x_0_xzz_xxxz, g_z_0_x_0_xzz_xxxzz, g_z_0_x_0_xzz_xxyy, g_z_0_x_0_xzz_xxyyy, g_z_0_x_0_xzz_xxyyz, g_z_0_x_0_xzz_xxyz, g_z_0_x_0_xzz_xxyzz, g_z_0_x_0_xzz_xxzz, g_z_0_x_0_xzz_xxzzz, g_z_0_x_0_xzz_xyyy, g_z_0_x_0_xzz_xyyyy, g_z_0_x_0_xzz_xyyyz, g_z_0_x_0_xzz_xyyz, g_z_0_x_0_xzz_xyyzz, g_z_0_x_0_xzz_xyzz, g_z_0_x_0_xzz_xyzzz, g_z_0_x_0_xzz_xzzz, g_z_0_x_0_xzz_xzzzz, g_z_0_x_0_xzz_yyyy, g_z_0_x_0_xzz_yyyz, g_z_0_x_0_xzz_yyzz, g_z_0_x_0_xzz_yzzz, g_z_0_x_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxzz_xxxx[k] = -g_z_0_x_0_xzz_xxxx[k] * ab_x + g_z_0_x_0_xzz_xxxxx[k];

                g_z_0_x_0_xxzz_xxxy[k] = -g_z_0_x_0_xzz_xxxy[k] * ab_x + g_z_0_x_0_xzz_xxxxy[k];

                g_z_0_x_0_xxzz_xxxz[k] = -g_z_0_x_0_xzz_xxxz[k] * ab_x + g_z_0_x_0_xzz_xxxxz[k];

                g_z_0_x_0_xxzz_xxyy[k] = -g_z_0_x_0_xzz_xxyy[k] * ab_x + g_z_0_x_0_xzz_xxxyy[k];

                g_z_0_x_0_xxzz_xxyz[k] = -g_z_0_x_0_xzz_xxyz[k] * ab_x + g_z_0_x_0_xzz_xxxyz[k];

                g_z_0_x_0_xxzz_xxzz[k] = -g_z_0_x_0_xzz_xxzz[k] * ab_x + g_z_0_x_0_xzz_xxxzz[k];

                g_z_0_x_0_xxzz_xyyy[k] = -g_z_0_x_0_xzz_xyyy[k] * ab_x + g_z_0_x_0_xzz_xxyyy[k];

                g_z_0_x_0_xxzz_xyyz[k] = -g_z_0_x_0_xzz_xyyz[k] * ab_x + g_z_0_x_0_xzz_xxyyz[k];

                g_z_0_x_0_xxzz_xyzz[k] = -g_z_0_x_0_xzz_xyzz[k] * ab_x + g_z_0_x_0_xzz_xxyzz[k];

                g_z_0_x_0_xxzz_xzzz[k] = -g_z_0_x_0_xzz_xzzz[k] * ab_x + g_z_0_x_0_xzz_xxzzz[k];

                g_z_0_x_0_xxzz_yyyy[k] = -g_z_0_x_0_xzz_yyyy[k] * ab_x + g_z_0_x_0_xzz_xyyyy[k];

                g_z_0_x_0_xxzz_yyyz[k] = -g_z_0_x_0_xzz_yyyz[k] * ab_x + g_z_0_x_0_xzz_xyyyz[k];

                g_z_0_x_0_xxzz_yyzz[k] = -g_z_0_x_0_xzz_yyzz[k] * ab_x + g_z_0_x_0_xzz_xyyzz[k];

                g_z_0_x_0_xxzz_yzzz[k] = -g_z_0_x_0_xzz_yzzz[k] * ab_x + g_z_0_x_0_xzz_xyzzz[k];

                g_z_0_x_0_xxzz_zzzz[k] = -g_z_0_x_0_xzz_zzzz[k] * ab_x + g_z_0_x_0_xzz_xzzzz[k];
            }

            /// Set up 1440-1455 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyy_xxxx = cbuffer.data(gg_geom_1010_off + 1440 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_xxxy = cbuffer.data(gg_geom_1010_off + 1441 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_xxxz = cbuffer.data(gg_geom_1010_off + 1442 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_xxyy = cbuffer.data(gg_geom_1010_off + 1443 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_xxyz = cbuffer.data(gg_geom_1010_off + 1444 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_xxzz = cbuffer.data(gg_geom_1010_off + 1445 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_xyyy = cbuffer.data(gg_geom_1010_off + 1446 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_xyyz = cbuffer.data(gg_geom_1010_off + 1447 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_xyzz = cbuffer.data(gg_geom_1010_off + 1448 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_xzzz = cbuffer.data(gg_geom_1010_off + 1449 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1450 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1451 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1452 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1453 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1454 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyy_xxxx, g_z_0_x_0_xyyy_xxxy, g_z_0_x_0_xyyy_xxxz, g_z_0_x_0_xyyy_xxyy, g_z_0_x_0_xyyy_xxyz, g_z_0_x_0_xyyy_xxzz, g_z_0_x_0_xyyy_xyyy, g_z_0_x_0_xyyy_xyyz, g_z_0_x_0_xyyy_xyzz, g_z_0_x_0_xyyy_xzzz, g_z_0_x_0_xyyy_yyyy, g_z_0_x_0_xyyy_yyyz, g_z_0_x_0_xyyy_yyzz, g_z_0_x_0_xyyy_yzzz, g_z_0_x_0_xyyy_zzzz, g_z_0_x_0_yyy_xxxx, g_z_0_x_0_yyy_xxxxx, g_z_0_x_0_yyy_xxxxy, g_z_0_x_0_yyy_xxxxz, g_z_0_x_0_yyy_xxxy, g_z_0_x_0_yyy_xxxyy, g_z_0_x_0_yyy_xxxyz, g_z_0_x_0_yyy_xxxz, g_z_0_x_0_yyy_xxxzz, g_z_0_x_0_yyy_xxyy, g_z_0_x_0_yyy_xxyyy, g_z_0_x_0_yyy_xxyyz, g_z_0_x_0_yyy_xxyz, g_z_0_x_0_yyy_xxyzz, g_z_0_x_0_yyy_xxzz, g_z_0_x_0_yyy_xxzzz, g_z_0_x_0_yyy_xyyy, g_z_0_x_0_yyy_xyyyy, g_z_0_x_0_yyy_xyyyz, g_z_0_x_0_yyy_xyyz, g_z_0_x_0_yyy_xyyzz, g_z_0_x_0_yyy_xyzz, g_z_0_x_0_yyy_xyzzz, g_z_0_x_0_yyy_xzzz, g_z_0_x_0_yyy_xzzzz, g_z_0_x_0_yyy_yyyy, g_z_0_x_0_yyy_yyyz, g_z_0_x_0_yyy_yyzz, g_z_0_x_0_yyy_yzzz, g_z_0_x_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyy_xxxx[k] = -g_z_0_x_0_yyy_xxxx[k] * ab_x + g_z_0_x_0_yyy_xxxxx[k];

                g_z_0_x_0_xyyy_xxxy[k] = -g_z_0_x_0_yyy_xxxy[k] * ab_x + g_z_0_x_0_yyy_xxxxy[k];

                g_z_0_x_0_xyyy_xxxz[k] = -g_z_0_x_0_yyy_xxxz[k] * ab_x + g_z_0_x_0_yyy_xxxxz[k];

                g_z_0_x_0_xyyy_xxyy[k] = -g_z_0_x_0_yyy_xxyy[k] * ab_x + g_z_0_x_0_yyy_xxxyy[k];

                g_z_0_x_0_xyyy_xxyz[k] = -g_z_0_x_0_yyy_xxyz[k] * ab_x + g_z_0_x_0_yyy_xxxyz[k];

                g_z_0_x_0_xyyy_xxzz[k] = -g_z_0_x_0_yyy_xxzz[k] * ab_x + g_z_0_x_0_yyy_xxxzz[k];

                g_z_0_x_0_xyyy_xyyy[k] = -g_z_0_x_0_yyy_xyyy[k] * ab_x + g_z_0_x_0_yyy_xxyyy[k];

                g_z_0_x_0_xyyy_xyyz[k] = -g_z_0_x_0_yyy_xyyz[k] * ab_x + g_z_0_x_0_yyy_xxyyz[k];

                g_z_0_x_0_xyyy_xyzz[k] = -g_z_0_x_0_yyy_xyzz[k] * ab_x + g_z_0_x_0_yyy_xxyzz[k];

                g_z_0_x_0_xyyy_xzzz[k] = -g_z_0_x_0_yyy_xzzz[k] * ab_x + g_z_0_x_0_yyy_xxzzz[k];

                g_z_0_x_0_xyyy_yyyy[k] = -g_z_0_x_0_yyy_yyyy[k] * ab_x + g_z_0_x_0_yyy_xyyyy[k];

                g_z_0_x_0_xyyy_yyyz[k] = -g_z_0_x_0_yyy_yyyz[k] * ab_x + g_z_0_x_0_yyy_xyyyz[k];

                g_z_0_x_0_xyyy_yyzz[k] = -g_z_0_x_0_yyy_yyzz[k] * ab_x + g_z_0_x_0_yyy_xyyzz[k];

                g_z_0_x_0_xyyy_yzzz[k] = -g_z_0_x_0_yyy_yzzz[k] * ab_x + g_z_0_x_0_yyy_xyzzz[k];

                g_z_0_x_0_xyyy_zzzz[k] = -g_z_0_x_0_yyy_zzzz[k] * ab_x + g_z_0_x_0_yyy_xzzzz[k];
            }

            /// Set up 1455-1470 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1455 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1456 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1457 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1458 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1459 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1460 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1461 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1462 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1463 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1464 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1465 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1466 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1467 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1468 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1469 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyz_xxxx, g_z_0_x_0_xyyz_xxxy, g_z_0_x_0_xyyz_xxxz, g_z_0_x_0_xyyz_xxyy, g_z_0_x_0_xyyz_xxyz, g_z_0_x_0_xyyz_xxzz, g_z_0_x_0_xyyz_xyyy, g_z_0_x_0_xyyz_xyyz, g_z_0_x_0_xyyz_xyzz, g_z_0_x_0_xyyz_xzzz, g_z_0_x_0_xyyz_yyyy, g_z_0_x_0_xyyz_yyyz, g_z_0_x_0_xyyz_yyzz, g_z_0_x_0_xyyz_yzzz, g_z_0_x_0_xyyz_zzzz, g_z_0_x_0_yyz_xxxx, g_z_0_x_0_yyz_xxxxx, g_z_0_x_0_yyz_xxxxy, g_z_0_x_0_yyz_xxxxz, g_z_0_x_0_yyz_xxxy, g_z_0_x_0_yyz_xxxyy, g_z_0_x_0_yyz_xxxyz, g_z_0_x_0_yyz_xxxz, g_z_0_x_0_yyz_xxxzz, g_z_0_x_0_yyz_xxyy, g_z_0_x_0_yyz_xxyyy, g_z_0_x_0_yyz_xxyyz, g_z_0_x_0_yyz_xxyz, g_z_0_x_0_yyz_xxyzz, g_z_0_x_0_yyz_xxzz, g_z_0_x_0_yyz_xxzzz, g_z_0_x_0_yyz_xyyy, g_z_0_x_0_yyz_xyyyy, g_z_0_x_0_yyz_xyyyz, g_z_0_x_0_yyz_xyyz, g_z_0_x_0_yyz_xyyzz, g_z_0_x_0_yyz_xyzz, g_z_0_x_0_yyz_xyzzz, g_z_0_x_0_yyz_xzzz, g_z_0_x_0_yyz_xzzzz, g_z_0_x_0_yyz_yyyy, g_z_0_x_0_yyz_yyyz, g_z_0_x_0_yyz_yyzz, g_z_0_x_0_yyz_yzzz, g_z_0_x_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyz_xxxx[k] = -g_z_0_x_0_yyz_xxxx[k] * ab_x + g_z_0_x_0_yyz_xxxxx[k];

                g_z_0_x_0_xyyz_xxxy[k] = -g_z_0_x_0_yyz_xxxy[k] * ab_x + g_z_0_x_0_yyz_xxxxy[k];

                g_z_0_x_0_xyyz_xxxz[k] = -g_z_0_x_0_yyz_xxxz[k] * ab_x + g_z_0_x_0_yyz_xxxxz[k];

                g_z_0_x_0_xyyz_xxyy[k] = -g_z_0_x_0_yyz_xxyy[k] * ab_x + g_z_0_x_0_yyz_xxxyy[k];

                g_z_0_x_0_xyyz_xxyz[k] = -g_z_0_x_0_yyz_xxyz[k] * ab_x + g_z_0_x_0_yyz_xxxyz[k];

                g_z_0_x_0_xyyz_xxzz[k] = -g_z_0_x_0_yyz_xxzz[k] * ab_x + g_z_0_x_0_yyz_xxxzz[k];

                g_z_0_x_0_xyyz_xyyy[k] = -g_z_0_x_0_yyz_xyyy[k] * ab_x + g_z_0_x_0_yyz_xxyyy[k];

                g_z_0_x_0_xyyz_xyyz[k] = -g_z_0_x_0_yyz_xyyz[k] * ab_x + g_z_0_x_0_yyz_xxyyz[k];

                g_z_0_x_0_xyyz_xyzz[k] = -g_z_0_x_0_yyz_xyzz[k] * ab_x + g_z_0_x_0_yyz_xxyzz[k];

                g_z_0_x_0_xyyz_xzzz[k] = -g_z_0_x_0_yyz_xzzz[k] * ab_x + g_z_0_x_0_yyz_xxzzz[k];

                g_z_0_x_0_xyyz_yyyy[k] = -g_z_0_x_0_yyz_yyyy[k] * ab_x + g_z_0_x_0_yyz_xyyyy[k];

                g_z_0_x_0_xyyz_yyyz[k] = -g_z_0_x_0_yyz_yyyz[k] * ab_x + g_z_0_x_0_yyz_xyyyz[k];

                g_z_0_x_0_xyyz_yyzz[k] = -g_z_0_x_0_yyz_yyzz[k] * ab_x + g_z_0_x_0_yyz_xyyzz[k];

                g_z_0_x_0_xyyz_yzzz[k] = -g_z_0_x_0_yyz_yzzz[k] * ab_x + g_z_0_x_0_yyz_xyzzz[k];

                g_z_0_x_0_xyyz_zzzz[k] = -g_z_0_x_0_yyz_zzzz[k] * ab_x + g_z_0_x_0_yyz_xzzzz[k];
            }

            /// Set up 1470-1485 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1470 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1471 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1472 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1473 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1474 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1475 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1476 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1477 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1478 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1479 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1480 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1481 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1482 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1483 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1484 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyzz_xxxx, g_z_0_x_0_xyzz_xxxy, g_z_0_x_0_xyzz_xxxz, g_z_0_x_0_xyzz_xxyy, g_z_0_x_0_xyzz_xxyz, g_z_0_x_0_xyzz_xxzz, g_z_0_x_0_xyzz_xyyy, g_z_0_x_0_xyzz_xyyz, g_z_0_x_0_xyzz_xyzz, g_z_0_x_0_xyzz_xzzz, g_z_0_x_0_xyzz_yyyy, g_z_0_x_0_xyzz_yyyz, g_z_0_x_0_xyzz_yyzz, g_z_0_x_0_xyzz_yzzz, g_z_0_x_0_xyzz_zzzz, g_z_0_x_0_yzz_xxxx, g_z_0_x_0_yzz_xxxxx, g_z_0_x_0_yzz_xxxxy, g_z_0_x_0_yzz_xxxxz, g_z_0_x_0_yzz_xxxy, g_z_0_x_0_yzz_xxxyy, g_z_0_x_0_yzz_xxxyz, g_z_0_x_0_yzz_xxxz, g_z_0_x_0_yzz_xxxzz, g_z_0_x_0_yzz_xxyy, g_z_0_x_0_yzz_xxyyy, g_z_0_x_0_yzz_xxyyz, g_z_0_x_0_yzz_xxyz, g_z_0_x_0_yzz_xxyzz, g_z_0_x_0_yzz_xxzz, g_z_0_x_0_yzz_xxzzz, g_z_0_x_0_yzz_xyyy, g_z_0_x_0_yzz_xyyyy, g_z_0_x_0_yzz_xyyyz, g_z_0_x_0_yzz_xyyz, g_z_0_x_0_yzz_xyyzz, g_z_0_x_0_yzz_xyzz, g_z_0_x_0_yzz_xyzzz, g_z_0_x_0_yzz_xzzz, g_z_0_x_0_yzz_xzzzz, g_z_0_x_0_yzz_yyyy, g_z_0_x_0_yzz_yyyz, g_z_0_x_0_yzz_yyzz, g_z_0_x_0_yzz_yzzz, g_z_0_x_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyzz_xxxx[k] = -g_z_0_x_0_yzz_xxxx[k] * ab_x + g_z_0_x_0_yzz_xxxxx[k];

                g_z_0_x_0_xyzz_xxxy[k] = -g_z_0_x_0_yzz_xxxy[k] * ab_x + g_z_0_x_0_yzz_xxxxy[k];

                g_z_0_x_0_xyzz_xxxz[k] = -g_z_0_x_0_yzz_xxxz[k] * ab_x + g_z_0_x_0_yzz_xxxxz[k];

                g_z_0_x_0_xyzz_xxyy[k] = -g_z_0_x_0_yzz_xxyy[k] * ab_x + g_z_0_x_0_yzz_xxxyy[k];

                g_z_0_x_0_xyzz_xxyz[k] = -g_z_0_x_0_yzz_xxyz[k] * ab_x + g_z_0_x_0_yzz_xxxyz[k];

                g_z_0_x_0_xyzz_xxzz[k] = -g_z_0_x_0_yzz_xxzz[k] * ab_x + g_z_0_x_0_yzz_xxxzz[k];

                g_z_0_x_0_xyzz_xyyy[k] = -g_z_0_x_0_yzz_xyyy[k] * ab_x + g_z_0_x_0_yzz_xxyyy[k];

                g_z_0_x_0_xyzz_xyyz[k] = -g_z_0_x_0_yzz_xyyz[k] * ab_x + g_z_0_x_0_yzz_xxyyz[k];

                g_z_0_x_0_xyzz_xyzz[k] = -g_z_0_x_0_yzz_xyzz[k] * ab_x + g_z_0_x_0_yzz_xxyzz[k];

                g_z_0_x_0_xyzz_xzzz[k] = -g_z_0_x_0_yzz_xzzz[k] * ab_x + g_z_0_x_0_yzz_xxzzz[k];

                g_z_0_x_0_xyzz_yyyy[k] = -g_z_0_x_0_yzz_yyyy[k] * ab_x + g_z_0_x_0_yzz_xyyyy[k];

                g_z_0_x_0_xyzz_yyyz[k] = -g_z_0_x_0_yzz_yyyz[k] * ab_x + g_z_0_x_0_yzz_xyyyz[k];

                g_z_0_x_0_xyzz_yyzz[k] = -g_z_0_x_0_yzz_yyzz[k] * ab_x + g_z_0_x_0_yzz_xyyzz[k];

                g_z_0_x_0_xyzz_yzzz[k] = -g_z_0_x_0_yzz_yzzz[k] * ab_x + g_z_0_x_0_yzz_xyzzz[k];

                g_z_0_x_0_xyzz_zzzz[k] = -g_z_0_x_0_yzz_zzzz[k] * ab_x + g_z_0_x_0_yzz_xzzzz[k];
            }

            /// Set up 1485-1500 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1485 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1486 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1487 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1488 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1489 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1490 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1491 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1492 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1493 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1494 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1495 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1496 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1497 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1498 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1499 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xzzz_xxxx, g_z_0_x_0_xzzz_xxxy, g_z_0_x_0_xzzz_xxxz, g_z_0_x_0_xzzz_xxyy, g_z_0_x_0_xzzz_xxyz, g_z_0_x_0_xzzz_xxzz, g_z_0_x_0_xzzz_xyyy, g_z_0_x_0_xzzz_xyyz, g_z_0_x_0_xzzz_xyzz, g_z_0_x_0_xzzz_xzzz, g_z_0_x_0_xzzz_yyyy, g_z_0_x_0_xzzz_yyyz, g_z_0_x_0_xzzz_yyzz, g_z_0_x_0_xzzz_yzzz, g_z_0_x_0_xzzz_zzzz, g_z_0_x_0_zzz_xxxx, g_z_0_x_0_zzz_xxxxx, g_z_0_x_0_zzz_xxxxy, g_z_0_x_0_zzz_xxxxz, g_z_0_x_0_zzz_xxxy, g_z_0_x_0_zzz_xxxyy, g_z_0_x_0_zzz_xxxyz, g_z_0_x_0_zzz_xxxz, g_z_0_x_0_zzz_xxxzz, g_z_0_x_0_zzz_xxyy, g_z_0_x_0_zzz_xxyyy, g_z_0_x_0_zzz_xxyyz, g_z_0_x_0_zzz_xxyz, g_z_0_x_0_zzz_xxyzz, g_z_0_x_0_zzz_xxzz, g_z_0_x_0_zzz_xxzzz, g_z_0_x_0_zzz_xyyy, g_z_0_x_0_zzz_xyyyy, g_z_0_x_0_zzz_xyyyz, g_z_0_x_0_zzz_xyyz, g_z_0_x_0_zzz_xyyzz, g_z_0_x_0_zzz_xyzz, g_z_0_x_0_zzz_xyzzz, g_z_0_x_0_zzz_xzzz, g_z_0_x_0_zzz_xzzzz, g_z_0_x_0_zzz_yyyy, g_z_0_x_0_zzz_yyyz, g_z_0_x_0_zzz_yyzz, g_z_0_x_0_zzz_yzzz, g_z_0_x_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xzzz_xxxx[k] = -g_z_0_x_0_zzz_xxxx[k] * ab_x + g_z_0_x_0_zzz_xxxxx[k];

                g_z_0_x_0_xzzz_xxxy[k] = -g_z_0_x_0_zzz_xxxy[k] * ab_x + g_z_0_x_0_zzz_xxxxy[k];

                g_z_0_x_0_xzzz_xxxz[k] = -g_z_0_x_0_zzz_xxxz[k] * ab_x + g_z_0_x_0_zzz_xxxxz[k];

                g_z_0_x_0_xzzz_xxyy[k] = -g_z_0_x_0_zzz_xxyy[k] * ab_x + g_z_0_x_0_zzz_xxxyy[k];

                g_z_0_x_0_xzzz_xxyz[k] = -g_z_0_x_0_zzz_xxyz[k] * ab_x + g_z_0_x_0_zzz_xxxyz[k];

                g_z_0_x_0_xzzz_xxzz[k] = -g_z_0_x_0_zzz_xxzz[k] * ab_x + g_z_0_x_0_zzz_xxxzz[k];

                g_z_0_x_0_xzzz_xyyy[k] = -g_z_0_x_0_zzz_xyyy[k] * ab_x + g_z_0_x_0_zzz_xxyyy[k];

                g_z_0_x_0_xzzz_xyyz[k] = -g_z_0_x_0_zzz_xyyz[k] * ab_x + g_z_0_x_0_zzz_xxyyz[k];

                g_z_0_x_0_xzzz_xyzz[k] = -g_z_0_x_0_zzz_xyzz[k] * ab_x + g_z_0_x_0_zzz_xxyzz[k];

                g_z_0_x_0_xzzz_xzzz[k] = -g_z_0_x_0_zzz_xzzz[k] * ab_x + g_z_0_x_0_zzz_xxzzz[k];

                g_z_0_x_0_xzzz_yyyy[k] = -g_z_0_x_0_zzz_yyyy[k] * ab_x + g_z_0_x_0_zzz_xyyyy[k];

                g_z_0_x_0_xzzz_yyyz[k] = -g_z_0_x_0_zzz_yyyz[k] * ab_x + g_z_0_x_0_zzz_xyyyz[k];

                g_z_0_x_0_xzzz_yyzz[k] = -g_z_0_x_0_zzz_yyzz[k] * ab_x + g_z_0_x_0_zzz_xyyzz[k];

                g_z_0_x_0_xzzz_yzzz[k] = -g_z_0_x_0_zzz_yzzz[k] * ab_x + g_z_0_x_0_zzz_xyzzz[k];

                g_z_0_x_0_xzzz_zzzz[k] = -g_z_0_x_0_zzz_zzzz[k] * ab_x + g_z_0_x_0_zzz_xzzzz[k];
            }

            /// Set up 1500-1515 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyy_xxxx = cbuffer.data(gg_geom_1010_off + 1500 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_xxxy = cbuffer.data(gg_geom_1010_off + 1501 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_xxxz = cbuffer.data(gg_geom_1010_off + 1502 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_xxyy = cbuffer.data(gg_geom_1010_off + 1503 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_xxyz = cbuffer.data(gg_geom_1010_off + 1504 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_xxzz = cbuffer.data(gg_geom_1010_off + 1505 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_xyyy = cbuffer.data(gg_geom_1010_off + 1506 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_xyyz = cbuffer.data(gg_geom_1010_off + 1507 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_xyzz = cbuffer.data(gg_geom_1010_off + 1508 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_xzzz = cbuffer.data(gg_geom_1010_off + 1509 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1510 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1511 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1512 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1513 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1514 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyy_xxxx, g_z_0_x_0_yyy_xxxxy, g_z_0_x_0_yyy_xxxy, g_z_0_x_0_yyy_xxxyy, g_z_0_x_0_yyy_xxxyz, g_z_0_x_0_yyy_xxxz, g_z_0_x_0_yyy_xxyy, g_z_0_x_0_yyy_xxyyy, g_z_0_x_0_yyy_xxyyz, g_z_0_x_0_yyy_xxyz, g_z_0_x_0_yyy_xxyzz, g_z_0_x_0_yyy_xxzz, g_z_0_x_0_yyy_xyyy, g_z_0_x_0_yyy_xyyyy, g_z_0_x_0_yyy_xyyyz, g_z_0_x_0_yyy_xyyz, g_z_0_x_0_yyy_xyyzz, g_z_0_x_0_yyy_xyzz, g_z_0_x_0_yyy_xyzzz, g_z_0_x_0_yyy_xzzz, g_z_0_x_0_yyy_yyyy, g_z_0_x_0_yyy_yyyyy, g_z_0_x_0_yyy_yyyyz, g_z_0_x_0_yyy_yyyz, g_z_0_x_0_yyy_yyyzz, g_z_0_x_0_yyy_yyzz, g_z_0_x_0_yyy_yyzzz, g_z_0_x_0_yyy_yzzz, g_z_0_x_0_yyy_yzzzz, g_z_0_x_0_yyy_zzzz, g_z_0_x_0_yyyy_xxxx, g_z_0_x_0_yyyy_xxxy, g_z_0_x_0_yyyy_xxxz, g_z_0_x_0_yyyy_xxyy, g_z_0_x_0_yyyy_xxyz, g_z_0_x_0_yyyy_xxzz, g_z_0_x_0_yyyy_xyyy, g_z_0_x_0_yyyy_xyyz, g_z_0_x_0_yyyy_xyzz, g_z_0_x_0_yyyy_xzzz, g_z_0_x_0_yyyy_yyyy, g_z_0_x_0_yyyy_yyyz, g_z_0_x_0_yyyy_yyzz, g_z_0_x_0_yyyy_yzzz, g_z_0_x_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyy_xxxx[k] = -g_z_0_x_0_yyy_xxxx[k] * ab_y + g_z_0_x_0_yyy_xxxxy[k];

                g_z_0_x_0_yyyy_xxxy[k] = -g_z_0_x_0_yyy_xxxy[k] * ab_y + g_z_0_x_0_yyy_xxxyy[k];

                g_z_0_x_0_yyyy_xxxz[k] = -g_z_0_x_0_yyy_xxxz[k] * ab_y + g_z_0_x_0_yyy_xxxyz[k];

                g_z_0_x_0_yyyy_xxyy[k] = -g_z_0_x_0_yyy_xxyy[k] * ab_y + g_z_0_x_0_yyy_xxyyy[k];

                g_z_0_x_0_yyyy_xxyz[k] = -g_z_0_x_0_yyy_xxyz[k] * ab_y + g_z_0_x_0_yyy_xxyyz[k];

                g_z_0_x_0_yyyy_xxzz[k] = -g_z_0_x_0_yyy_xxzz[k] * ab_y + g_z_0_x_0_yyy_xxyzz[k];

                g_z_0_x_0_yyyy_xyyy[k] = -g_z_0_x_0_yyy_xyyy[k] * ab_y + g_z_0_x_0_yyy_xyyyy[k];

                g_z_0_x_0_yyyy_xyyz[k] = -g_z_0_x_0_yyy_xyyz[k] * ab_y + g_z_0_x_0_yyy_xyyyz[k];

                g_z_0_x_0_yyyy_xyzz[k] = -g_z_0_x_0_yyy_xyzz[k] * ab_y + g_z_0_x_0_yyy_xyyzz[k];

                g_z_0_x_0_yyyy_xzzz[k] = -g_z_0_x_0_yyy_xzzz[k] * ab_y + g_z_0_x_0_yyy_xyzzz[k];

                g_z_0_x_0_yyyy_yyyy[k] = -g_z_0_x_0_yyy_yyyy[k] * ab_y + g_z_0_x_0_yyy_yyyyy[k];

                g_z_0_x_0_yyyy_yyyz[k] = -g_z_0_x_0_yyy_yyyz[k] * ab_y + g_z_0_x_0_yyy_yyyyz[k];

                g_z_0_x_0_yyyy_yyzz[k] = -g_z_0_x_0_yyy_yyzz[k] * ab_y + g_z_0_x_0_yyy_yyyzz[k];

                g_z_0_x_0_yyyy_yzzz[k] = -g_z_0_x_0_yyy_yzzz[k] * ab_y + g_z_0_x_0_yyy_yyzzz[k];

                g_z_0_x_0_yyyy_zzzz[k] = -g_z_0_x_0_yyy_zzzz[k] * ab_y + g_z_0_x_0_yyy_yzzzz[k];
            }

            /// Set up 1515-1530 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1515 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1516 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1517 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1518 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1519 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1520 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1521 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1522 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1523 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1524 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1525 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1526 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1527 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1528 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1529 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyz_xxxx, g_z_0_x_0_yyyz_xxxy, g_z_0_x_0_yyyz_xxxz, g_z_0_x_0_yyyz_xxyy, g_z_0_x_0_yyyz_xxyz, g_z_0_x_0_yyyz_xxzz, g_z_0_x_0_yyyz_xyyy, g_z_0_x_0_yyyz_xyyz, g_z_0_x_0_yyyz_xyzz, g_z_0_x_0_yyyz_xzzz, g_z_0_x_0_yyyz_yyyy, g_z_0_x_0_yyyz_yyyz, g_z_0_x_0_yyyz_yyzz, g_z_0_x_0_yyyz_yzzz, g_z_0_x_0_yyyz_zzzz, g_z_0_x_0_yyz_xxxx, g_z_0_x_0_yyz_xxxxy, g_z_0_x_0_yyz_xxxy, g_z_0_x_0_yyz_xxxyy, g_z_0_x_0_yyz_xxxyz, g_z_0_x_0_yyz_xxxz, g_z_0_x_0_yyz_xxyy, g_z_0_x_0_yyz_xxyyy, g_z_0_x_0_yyz_xxyyz, g_z_0_x_0_yyz_xxyz, g_z_0_x_0_yyz_xxyzz, g_z_0_x_0_yyz_xxzz, g_z_0_x_0_yyz_xyyy, g_z_0_x_0_yyz_xyyyy, g_z_0_x_0_yyz_xyyyz, g_z_0_x_0_yyz_xyyz, g_z_0_x_0_yyz_xyyzz, g_z_0_x_0_yyz_xyzz, g_z_0_x_0_yyz_xyzzz, g_z_0_x_0_yyz_xzzz, g_z_0_x_0_yyz_yyyy, g_z_0_x_0_yyz_yyyyy, g_z_0_x_0_yyz_yyyyz, g_z_0_x_0_yyz_yyyz, g_z_0_x_0_yyz_yyyzz, g_z_0_x_0_yyz_yyzz, g_z_0_x_0_yyz_yyzzz, g_z_0_x_0_yyz_yzzz, g_z_0_x_0_yyz_yzzzz, g_z_0_x_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyz_xxxx[k] = -g_z_0_x_0_yyz_xxxx[k] * ab_y + g_z_0_x_0_yyz_xxxxy[k];

                g_z_0_x_0_yyyz_xxxy[k] = -g_z_0_x_0_yyz_xxxy[k] * ab_y + g_z_0_x_0_yyz_xxxyy[k];

                g_z_0_x_0_yyyz_xxxz[k] = -g_z_0_x_0_yyz_xxxz[k] * ab_y + g_z_0_x_0_yyz_xxxyz[k];

                g_z_0_x_0_yyyz_xxyy[k] = -g_z_0_x_0_yyz_xxyy[k] * ab_y + g_z_0_x_0_yyz_xxyyy[k];

                g_z_0_x_0_yyyz_xxyz[k] = -g_z_0_x_0_yyz_xxyz[k] * ab_y + g_z_0_x_0_yyz_xxyyz[k];

                g_z_0_x_0_yyyz_xxzz[k] = -g_z_0_x_0_yyz_xxzz[k] * ab_y + g_z_0_x_0_yyz_xxyzz[k];

                g_z_0_x_0_yyyz_xyyy[k] = -g_z_0_x_0_yyz_xyyy[k] * ab_y + g_z_0_x_0_yyz_xyyyy[k];

                g_z_0_x_0_yyyz_xyyz[k] = -g_z_0_x_0_yyz_xyyz[k] * ab_y + g_z_0_x_0_yyz_xyyyz[k];

                g_z_0_x_0_yyyz_xyzz[k] = -g_z_0_x_0_yyz_xyzz[k] * ab_y + g_z_0_x_0_yyz_xyyzz[k];

                g_z_0_x_0_yyyz_xzzz[k] = -g_z_0_x_0_yyz_xzzz[k] * ab_y + g_z_0_x_0_yyz_xyzzz[k];

                g_z_0_x_0_yyyz_yyyy[k] = -g_z_0_x_0_yyz_yyyy[k] * ab_y + g_z_0_x_0_yyz_yyyyy[k];

                g_z_0_x_0_yyyz_yyyz[k] = -g_z_0_x_0_yyz_yyyz[k] * ab_y + g_z_0_x_0_yyz_yyyyz[k];

                g_z_0_x_0_yyyz_yyzz[k] = -g_z_0_x_0_yyz_yyzz[k] * ab_y + g_z_0_x_0_yyz_yyyzz[k];

                g_z_0_x_0_yyyz_yzzz[k] = -g_z_0_x_0_yyz_yzzz[k] * ab_y + g_z_0_x_0_yyz_yyzzz[k];

                g_z_0_x_0_yyyz_zzzz[k] = -g_z_0_x_0_yyz_zzzz[k] * ab_y + g_z_0_x_0_yyz_yzzzz[k];
            }

            /// Set up 1530-1545 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1530 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1531 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1532 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1533 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1534 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1535 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1536 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1537 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1538 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1539 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1540 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1541 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1542 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1543 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1544 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyzz_xxxx, g_z_0_x_0_yyzz_xxxy, g_z_0_x_0_yyzz_xxxz, g_z_0_x_0_yyzz_xxyy, g_z_0_x_0_yyzz_xxyz, g_z_0_x_0_yyzz_xxzz, g_z_0_x_0_yyzz_xyyy, g_z_0_x_0_yyzz_xyyz, g_z_0_x_0_yyzz_xyzz, g_z_0_x_0_yyzz_xzzz, g_z_0_x_0_yyzz_yyyy, g_z_0_x_0_yyzz_yyyz, g_z_0_x_0_yyzz_yyzz, g_z_0_x_0_yyzz_yzzz, g_z_0_x_0_yyzz_zzzz, g_z_0_x_0_yzz_xxxx, g_z_0_x_0_yzz_xxxxy, g_z_0_x_0_yzz_xxxy, g_z_0_x_0_yzz_xxxyy, g_z_0_x_0_yzz_xxxyz, g_z_0_x_0_yzz_xxxz, g_z_0_x_0_yzz_xxyy, g_z_0_x_0_yzz_xxyyy, g_z_0_x_0_yzz_xxyyz, g_z_0_x_0_yzz_xxyz, g_z_0_x_0_yzz_xxyzz, g_z_0_x_0_yzz_xxzz, g_z_0_x_0_yzz_xyyy, g_z_0_x_0_yzz_xyyyy, g_z_0_x_0_yzz_xyyyz, g_z_0_x_0_yzz_xyyz, g_z_0_x_0_yzz_xyyzz, g_z_0_x_0_yzz_xyzz, g_z_0_x_0_yzz_xyzzz, g_z_0_x_0_yzz_xzzz, g_z_0_x_0_yzz_yyyy, g_z_0_x_0_yzz_yyyyy, g_z_0_x_0_yzz_yyyyz, g_z_0_x_0_yzz_yyyz, g_z_0_x_0_yzz_yyyzz, g_z_0_x_0_yzz_yyzz, g_z_0_x_0_yzz_yyzzz, g_z_0_x_0_yzz_yzzz, g_z_0_x_0_yzz_yzzzz, g_z_0_x_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyzz_xxxx[k] = -g_z_0_x_0_yzz_xxxx[k] * ab_y + g_z_0_x_0_yzz_xxxxy[k];

                g_z_0_x_0_yyzz_xxxy[k] = -g_z_0_x_0_yzz_xxxy[k] * ab_y + g_z_0_x_0_yzz_xxxyy[k];

                g_z_0_x_0_yyzz_xxxz[k] = -g_z_0_x_0_yzz_xxxz[k] * ab_y + g_z_0_x_0_yzz_xxxyz[k];

                g_z_0_x_0_yyzz_xxyy[k] = -g_z_0_x_0_yzz_xxyy[k] * ab_y + g_z_0_x_0_yzz_xxyyy[k];

                g_z_0_x_0_yyzz_xxyz[k] = -g_z_0_x_0_yzz_xxyz[k] * ab_y + g_z_0_x_0_yzz_xxyyz[k];

                g_z_0_x_0_yyzz_xxzz[k] = -g_z_0_x_0_yzz_xxzz[k] * ab_y + g_z_0_x_0_yzz_xxyzz[k];

                g_z_0_x_0_yyzz_xyyy[k] = -g_z_0_x_0_yzz_xyyy[k] * ab_y + g_z_0_x_0_yzz_xyyyy[k];

                g_z_0_x_0_yyzz_xyyz[k] = -g_z_0_x_0_yzz_xyyz[k] * ab_y + g_z_0_x_0_yzz_xyyyz[k];

                g_z_0_x_0_yyzz_xyzz[k] = -g_z_0_x_0_yzz_xyzz[k] * ab_y + g_z_0_x_0_yzz_xyyzz[k];

                g_z_0_x_0_yyzz_xzzz[k] = -g_z_0_x_0_yzz_xzzz[k] * ab_y + g_z_0_x_0_yzz_xyzzz[k];

                g_z_0_x_0_yyzz_yyyy[k] = -g_z_0_x_0_yzz_yyyy[k] * ab_y + g_z_0_x_0_yzz_yyyyy[k];

                g_z_0_x_0_yyzz_yyyz[k] = -g_z_0_x_0_yzz_yyyz[k] * ab_y + g_z_0_x_0_yzz_yyyyz[k];

                g_z_0_x_0_yyzz_yyzz[k] = -g_z_0_x_0_yzz_yyzz[k] * ab_y + g_z_0_x_0_yzz_yyyzz[k];

                g_z_0_x_0_yyzz_yzzz[k] = -g_z_0_x_0_yzz_yzzz[k] * ab_y + g_z_0_x_0_yzz_yyzzz[k];

                g_z_0_x_0_yyzz_zzzz[k] = -g_z_0_x_0_yzz_zzzz[k] * ab_y + g_z_0_x_0_yzz_yzzzz[k];
            }

            /// Set up 1545-1560 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1545 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1546 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1547 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1548 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1549 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1550 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1551 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1552 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1553 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1554 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1555 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1556 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1557 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1558 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1559 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yzzz_xxxx, g_z_0_x_0_yzzz_xxxy, g_z_0_x_0_yzzz_xxxz, g_z_0_x_0_yzzz_xxyy, g_z_0_x_0_yzzz_xxyz, g_z_0_x_0_yzzz_xxzz, g_z_0_x_0_yzzz_xyyy, g_z_0_x_0_yzzz_xyyz, g_z_0_x_0_yzzz_xyzz, g_z_0_x_0_yzzz_xzzz, g_z_0_x_0_yzzz_yyyy, g_z_0_x_0_yzzz_yyyz, g_z_0_x_0_yzzz_yyzz, g_z_0_x_0_yzzz_yzzz, g_z_0_x_0_yzzz_zzzz, g_z_0_x_0_zzz_xxxx, g_z_0_x_0_zzz_xxxxy, g_z_0_x_0_zzz_xxxy, g_z_0_x_0_zzz_xxxyy, g_z_0_x_0_zzz_xxxyz, g_z_0_x_0_zzz_xxxz, g_z_0_x_0_zzz_xxyy, g_z_0_x_0_zzz_xxyyy, g_z_0_x_0_zzz_xxyyz, g_z_0_x_0_zzz_xxyz, g_z_0_x_0_zzz_xxyzz, g_z_0_x_0_zzz_xxzz, g_z_0_x_0_zzz_xyyy, g_z_0_x_0_zzz_xyyyy, g_z_0_x_0_zzz_xyyyz, g_z_0_x_0_zzz_xyyz, g_z_0_x_0_zzz_xyyzz, g_z_0_x_0_zzz_xyzz, g_z_0_x_0_zzz_xyzzz, g_z_0_x_0_zzz_xzzz, g_z_0_x_0_zzz_yyyy, g_z_0_x_0_zzz_yyyyy, g_z_0_x_0_zzz_yyyyz, g_z_0_x_0_zzz_yyyz, g_z_0_x_0_zzz_yyyzz, g_z_0_x_0_zzz_yyzz, g_z_0_x_0_zzz_yyzzz, g_z_0_x_0_zzz_yzzz, g_z_0_x_0_zzz_yzzzz, g_z_0_x_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yzzz_xxxx[k] = -g_z_0_x_0_zzz_xxxx[k] * ab_y + g_z_0_x_0_zzz_xxxxy[k];

                g_z_0_x_0_yzzz_xxxy[k] = -g_z_0_x_0_zzz_xxxy[k] * ab_y + g_z_0_x_0_zzz_xxxyy[k];

                g_z_0_x_0_yzzz_xxxz[k] = -g_z_0_x_0_zzz_xxxz[k] * ab_y + g_z_0_x_0_zzz_xxxyz[k];

                g_z_0_x_0_yzzz_xxyy[k] = -g_z_0_x_0_zzz_xxyy[k] * ab_y + g_z_0_x_0_zzz_xxyyy[k];

                g_z_0_x_0_yzzz_xxyz[k] = -g_z_0_x_0_zzz_xxyz[k] * ab_y + g_z_0_x_0_zzz_xxyyz[k];

                g_z_0_x_0_yzzz_xxzz[k] = -g_z_0_x_0_zzz_xxzz[k] * ab_y + g_z_0_x_0_zzz_xxyzz[k];

                g_z_0_x_0_yzzz_xyyy[k] = -g_z_0_x_0_zzz_xyyy[k] * ab_y + g_z_0_x_0_zzz_xyyyy[k];

                g_z_0_x_0_yzzz_xyyz[k] = -g_z_0_x_0_zzz_xyyz[k] * ab_y + g_z_0_x_0_zzz_xyyyz[k];

                g_z_0_x_0_yzzz_xyzz[k] = -g_z_0_x_0_zzz_xyzz[k] * ab_y + g_z_0_x_0_zzz_xyyzz[k];

                g_z_0_x_0_yzzz_xzzz[k] = -g_z_0_x_0_zzz_xzzz[k] * ab_y + g_z_0_x_0_zzz_xyzzz[k];

                g_z_0_x_0_yzzz_yyyy[k] = -g_z_0_x_0_zzz_yyyy[k] * ab_y + g_z_0_x_0_zzz_yyyyy[k];

                g_z_0_x_0_yzzz_yyyz[k] = -g_z_0_x_0_zzz_yyyz[k] * ab_y + g_z_0_x_0_zzz_yyyyz[k];

                g_z_0_x_0_yzzz_yyzz[k] = -g_z_0_x_0_zzz_yyzz[k] * ab_y + g_z_0_x_0_zzz_yyyzz[k];

                g_z_0_x_0_yzzz_yzzz[k] = -g_z_0_x_0_zzz_yzzz[k] * ab_y + g_z_0_x_0_zzz_yyzzz[k];

                g_z_0_x_0_yzzz_zzzz[k] = -g_z_0_x_0_zzz_zzzz[k] * ab_y + g_z_0_x_0_zzz_yzzzz[k];
            }

            /// Set up 1560-1575 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_zzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1560 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1561 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1562 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1563 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1564 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1565 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1566 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1567 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1568 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1569 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1570 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1571 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1572 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1573 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1574 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_zzz_xxxx, g_0_0_x_0_zzz_xxxy, g_0_0_x_0_zzz_xxxz, g_0_0_x_0_zzz_xxyy, g_0_0_x_0_zzz_xxyz, g_0_0_x_0_zzz_xxzz, g_0_0_x_0_zzz_xyyy, g_0_0_x_0_zzz_xyyz, g_0_0_x_0_zzz_xyzz, g_0_0_x_0_zzz_xzzz, g_0_0_x_0_zzz_yyyy, g_0_0_x_0_zzz_yyyz, g_0_0_x_0_zzz_yyzz, g_0_0_x_0_zzz_yzzz, g_0_0_x_0_zzz_zzzz, g_z_0_x_0_zzz_xxxx, g_z_0_x_0_zzz_xxxxz, g_z_0_x_0_zzz_xxxy, g_z_0_x_0_zzz_xxxyz, g_z_0_x_0_zzz_xxxz, g_z_0_x_0_zzz_xxxzz, g_z_0_x_0_zzz_xxyy, g_z_0_x_0_zzz_xxyyz, g_z_0_x_0_zzz_xxyz, g_z_0_x_0_zzz_xxyzz, g_z_0_x_0_zzz_xxzz, g_z_0_x_0_zzz_xxzzz, g_z_0_x_0_zzz_xyyy, g_z_0_x_0_zzz_xyyyz, g_z_0_x_0_zzz_xyyz, g_z_0_x_0_zzz_xyyzz, g_z_0_x_0_zzz_xyzz, g_z_0_x_0_zzz_xyzzz, g_z_0_x_0_zzz_xzzz, g_z_0_x_0_zzz_xzzzz, g_z_0_x_0_zzz_yyyy, g_z_0_x_0_zzz_yyyyz, g_z_0_x_0_zzz_yyyz, g_z_0_x_0_zzz_yyyzz, g_z_0_x_0_zzz_yyzz, g_z_0_x_0_zzz_yyzzz, g_z_0_x_0_zzz_yzzz, g_z_0_x_0_zzz_yzzzz, g_z_0_x_0_zzz_zzzz, g_z_0_x_0_zzz_zzzzz, g_z_0_x_0_zzzz_xxxx, g_z_0_x_0_zzzz_xxxy, g_z_0_x_0_zzzz_xxxz, g_z_0_x_0_zzzz_xxyy, g_z_0_x_0_zzzz_xxyz, g_z_0_x_0_zzzz_xxzz, g_z_0_x_0_zzzz_xyyy, g_z_0_x_0_zzzz_xyyz, g_z_0_x_0_zzzz_xyzz, g_z_0_x_0_zzzz_xzzz, g_z_0_x_0_zzzz_yyyy, g_z_0_x_0_zzzz_yyyz, g_z_0_x_0_zzzz_yyzz, g_z_0_x_0_zzzz_yzzz, g_z_0_x_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_zzzz_xxxx[k] = -g_0_0_x_0_zzz_xxxx[k] - g_z_0_x_0_zzz_xxxx[k] * ab_z + g_z_0_x_0_zzz_xxxxz[k];

                g_z_0_x_0_zzzz_xxxy[k] = -g_0_0_x_0_zzz_xxxy[k] - g_z_0_x_0_zzz_xxxy[k] * ab_z + g_z_0_x_0_zzz_xxxyz[k];

                g_z_0_x_0_zzzz_xxxz[k] = -g_0_0_x_0_zzz_xxxz[k] - g_z_0_x_0_zzz_xxxz[k] * ab_z + g_z_0_x_0_zzz_xxxzz[k];

                g_z_0_x_0_zzzz_xxyy[k] = -g_0_0_x_0_zzz_xxyy[k] - g_z_0_x_0_zzz_xxyy[k] * ab_z + g_z_0_x_0_zzz_xxyyz[k];

                g_z_0_x_0_zzzz_xxyz[k] = -g_0_0_x_0_zzz_xxyz[k] - g_z_0_x_0_zzz_xxyz[k] * ab_z + g_z_0_x_0_zzz_xxyzz[k];

                g_z_0_x_0_zzzz_xxzz[k] = -g_0_0_x_0_zzz_xxzz[k] - g_z_0_x_0_zzz_xxzz[k] * ab_z + g_z_0_x_0_zzz_xxzzz[k];

                g_z_0_x_0_zzzz_xyyy[k] = -g_0_0_x_0_zzz_xyyy[k] - g_z_0_x_0_zzz_xyyy[k] * ab_z + g_z_0_x_0_zzz_xyyyz[k];

                g_z_0_x_0_zzzz_xyyz[k] = -g_0_0_x_0_zzz_xyyz[k] - g_z_0_x_0_zzz_xyyz[k] * ab_z + g_z_0_x_0_zzz_xyyzz[k];

                g_z_0_x_0_zzzz_xyzz[k] = -g_0_0_x_0_zzz_xyzz[k] - g_z_0_x_0_zzz_xyzz[k] * ab_z + g_z_0_x_0_zzz_xyzzz[k];

                g_z_0_x_0_zzzz_xzzz[k] = -g_0_0_x_0_zzz_xzzz[k] - g_z_0_x_0_zzz_xzzz[k] * ab_z + g_z_0_x_0_zzz_xzzzz[k];

                g_z_0_x_0_zzzz_yyyy[k] = -g_0_0_x_0_zzz_yyyy[k] - g_z_0_x_0_zzz_yyyy[k] * ab_z + g_z_0_x_0_zzz_yyyyz[k];

                g_z_0_x_0_zzzz_yyyz[k] = -g_0_0_x_0_zzz_yyyz[k] - g_z_0_x_0_zzz_yyyz[k] * ab_z + g_z_0_x_0_zzz_yyyzz[k];

                g_z_0_x_0_zzzz_yyzz[k] = -g_0_0_x_0_zzz_yyzz[k] - g_z_0_x_0_zzz_yyzz[k] * ab_z + g_z_0_x_0_zzz_yyzzz[k];

                g_z_0_x_0_zzzz_yzzz[k] = -g_0_0_x_0_zzz_yzzz[k] - g_z_0_x_0_zzz_yzzz[k] * ab_z + g_z_0_x_0_zzz_yzzzz[k];

                g_z_0_x_0_zzzz_zzzz[k] = -g_0_0_x_0_zzz_zzzz[k] - g_z_0_x_0_zzz_zzzz[k] * ab_z + g_z_0_x_0_zzz_zzzzz[k];
            }

            /// Set up 1575-1590 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxx_xxxx = cbuffer.data(gg_geom_1010_off + 1575 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_xxxy = cbuffer.data(gg_geom_1010_off + 1576 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_xxxz = cbuffer.data(gg_geom_1010_off + 1577 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_xxyy = cbuffer.data(gg_geom_1010_off + 1578 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_xxyz = cbuffer.data(gg_geom_1010_off + 1579 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_xxzz = cbuffer.data(gg_geom_1010_off + 1580 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_xyyy = cbuffer.data(gg_geom_1010_off + 1581 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_xyyz = cbuffer.data(gg_geom_1010_off + 1582 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_xyzz = cbuffer.data(gg_geom_1010_off + 1583 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_xzzz = cbuffer.data(gg_geom_1010_off + 1584 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_yyyy = cbuffer.data(gg_geom_1010_off + 1585 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_yyyz = cbuffer.data(gg_geom_1010_off + 1586 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_yyzz = cbuffer.data(gg_geom_1010_off + 1587 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_yzzz = cbuffer.data(gg_geom_1010_off + 1588 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_zzzz = cbuffer.data(gg_geom_1010_off + 1589 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxx_xxxx, g_z_0_y_0_xxx_xxxxx, g_z_0_y_0_xxx_xxxxy, g_z_0_y_0_xxx_xxxxz, g_z_0_y_0_xxx_xxxy, g_z_0_y_0_xxx_xxxyy, g_z_0_y_0_xxx_xxxyz, g_z_0_y_0_xxx_xxxz, g_z_0_y_0_xxx_xxxzz, g_z_0_y_0_xxx_xxyy, g_z_0_y_0_xxx_xxyyy, g_z_0_y_0_xxx_xxyyz, g_z_0_y_0_xxx_xxyz, g_z_0_y_0_xxx_xxyzz, g_z_0_y_0_xxx_xxzz, g_z_0_y_0_xxx_xxzzz, g_z_0_y_0_xxx_xyyy, g_z_0_y_0_xxx_xyyyy, g_z_0_y_0_xxx_xyyyz, g_z_0_y_0_xxx_xyyz, g_z_0_y_0_xxx_xyyzz, g_z_0_y_0_xxx_xyzz, g_z_0_y_0_xxx_xyzzz, g_z_0_y_0_xxx_xzzz, g_z_0_y_0_xxx_xzzzz, g_z_0_y_0_xxx_yyyy, g_z_0_y_0_xxx_yyyz, g_z_0_y_0_xxx_yyzz, g_z_0_y_0_xxx_yzzz, g_z_0_y_0_xxx_zzzz, g_z_0_y_0_xxxx_xxxx, g_z_0_y_0_xxxx_xxxy, g_z_0_y_0_xxxx_xxxz, g_z_0_y_0_xxxx_xxyy, g_z_0_y_0_xxxx_xxyz, g_z_0_y_0_xxxx_xxzz, g_z_0_y_0_xxxx_xyyy, g_z_0_y_0_xxxx_xyyz, g_z_0_y_0_xxxx_xyzz, g_z_0_y_0_xxxx_xzzz, g_z_0_y_0_xxxx_yyyy, g_z_0_y_0_xxxx_yyyz, g_z_0_y_0_xxxx_yyzz, g_z_0_y_0_xxxx_yzzz, g_z_0_y_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxx_xxxx[k] = -g_z_0_y_0_xxx_xxxx[k] * ab_x + g_z_0_y_0_xxx_xxxxx[k];

                g_z_0_y_0_xxxx_xxxy[k] = -g_z_0_y_0_xxx_xxxy[k] * ab_x + g_z_0_y_0_xxx_xxxxy[k];

                g_z_0_y_0_xxxx_xxxz[k] = -g_z_0_y_0_xxx_xxxz[k] * ab_x + g_z_0_y_0_xxx_xxxxz[k];

                g_z_0_y_0_xxxx_xxyy[k] = -g_z_0_y_0_xxx_xxyy[k] * ab_x + g_z_0_y_0_xxx_xxxyy[k];

                g_z_0_y_0_xxxx_xxyz[k] = -g_z_0_y_0_xxx_xxyz[k] * ab_x + g_z_0_y_0_xxx_xxxyz[k];

                g_z_0_y_0_xxxx_xxzz[k] = -g_z_0_y_0_xxx_xxzz[k] * ab_x + g_z_0_y_0_xxx_xxxzz[k];

                g_z_0_y_0_xxxx_xyyy[k] = -g_z_0_y_0_xxx_xyyy[k] * ab_x + g_z_0_y_0_xxx_xxyyy[k];

                g_z_0_y_0_xxxx_xyyz[k] = -g_z_0_y_0_xxx_xyyz[k] * ab_x + g_z_0_y_0_xxx_xxyyz[k];

                g_z_0_y_0_xxxx_xyzz[k] = -g_z_0_y_0_xxx_xyzz[k] * ab_x + g_z_0_y_0_xxx_xxyzz[k];

                g_z_0_y_0_xxxx_xzzz[k] = -g_z_0_y_0_xxx_xzzz[k] * ab_x + g_z_0_y_0_xxx_xxzzz[k];

                g_z_0_y_0_xxxx_yyyy[k] = -g_z_0_y_0_xxx_yyyy[k] * ab_x + g_z_0_y_0_xxx_xyyyy[k];

                g_z_0_y_0_xxxx_yyyz[k] = -g_z_0_y_0_xxx_yyyz[k] * ab_x + g_z_0_y_0_xxx_xyyyz[k];

                g_z_0_y_0_xxxx_yyzz[k] = -g_z_0_y_0_xxx_yyzz[k] * ab_x + g_z_0_y_0_xxx_xyyzz[k];

                g_z_0_y_0_xxxx_yzzz[k] = -g_z_0_y_0_xxx_yzzz[k] * ab_x + g_z_0_y_0_xxx_xyzzz[k];

                g_z_0_y_0_xxxx_zzzz[k] = -g_z_0_y_0_xxx_zzzz[k] * ab_x + g_z_0_y_0_xxx_xzzzz[k];
            }

            /// Set up 1590-1605 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxy_xxxx = cbuffer.data(gg_geom_1010_off + 1590 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_xxxy = cbuffer.data(gg_geom_1010_off + 1591 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_xxxz = cbuffer.data(gg_geom_1010_off + 1592 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_xxyy = cbuffer.data(gg_geom_1010_off + 1593 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_xxyz = cbuffer.data(gg_geom_1010_off + 1594 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_xxzz = cbuffer.data(gg_geom_1010_off + 1595 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_xyyy = cbuffer.data(gg_geom_1010_off + 1596 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_xyyz = cbuffer.data(gg_geom_1010_off + 1597 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_xyzz = cbuffer.data(gg_geom_1010_off + 1598 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_xzzz = cbuffer.data(gg_geom_1010_off + 1599 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_yyyy = cbuffer.data(gg_geom_1010_off + 1600 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_yyyz = cbuffer.data(gg_geom_1010_off + 1601 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_yyzz = cbuffer.data(gg_geom_1010_off + 1602 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_yzzz = cbuffer.data(gg_geom_1010_off + 1603 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_zzzz = cbuffer.data(gg_geom_1010_off + 1604 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxy_xxxx, g_z_0_y_0_xxxy_xxxy, g_z_0_y_0_xxxy_xxxz, g_z_0_y_0_xxxy_xxyy, g_z_0_y_0_xxxy_xxyz, g_z_0_y_0_xxxy_xxzz, g_z_0_y_0_xxxy_xyyy, g_z_0_y_0_xxxy_xyyz, g_z_0_y_0_xxxy_xyzz, g_z_0_y_0_xxxy_xzzz, g_z_0_y_0_xxxy_yyyy, g_z_0_y_0_xxxy_yyyz, g_z_0_y_0_xxxy_yyzz, g_z_0_y_0_xxxy_yzzz, g_z_0_y_0_xxxy_zzzz, g_z_0_y_0_xxy_xxxx, g_z_0_y_0_xxy_xxxxx, g_z_0_y_0_xxy_xxxxy, g_z_0_y_0_xxy_xxxxz, g_z_0_y_0_xxy_xxxy, g_z_0_y_0_xxy_xxxyy, g_z_0_y_0_xxy_xxxyz, g_z_0_y_0_xxy_xxxz, g_z_0_y_0_xxy_xxxzz, g_z_0_y_0_xxy_xxyy, g_z_0_y_0_xxy_xxyyy, g_z_0_y_0_xxy_xxyyz, g_z_0_y_0_xxy_xxyz, g_z_0_y_0_xxy_xxyzz, g_z_0_y_0_xxy_xxzz, g_z_0_y_0_xxy_xxzzz, g_z_0_y_0_xxy_xyyy, g_z_0_y_0_xxy_xyyyy, g_z_0_y_0_xxy_xyyyz, g_z_0_y_0_xxy_xyyz, g_z_0_y_0_xxy_xyyzz, g_z_0_y_0_xxy_xyzz, g_z_0_y_0_xxy_xyzzz, g_z_0_y_0_xxy_xzzz, g_z_0_y_0_xxy_xzzzz, g_z_0_y_0_xxy_yyyy, g_z_0_y_0_xxy_yyyz, g_z_0_y_0_xxy_yyzz, g_z_0_y_0_xxy_yzzz, g_z_0_y_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxy_xxxx[k] = -g_z_0_y_0_xxy_xxxx[k] * ab_x + g_z_0_y_0_xxy_xxxxx[k];

                g_z_0_y_0_xxxy_xxxy[k] = -g_z_0_y_0_xxy_xxxy[k] * ab_x + g_z_0_y_0_xxy_xxxxy[k];

                g_z_0_y_0_xxxy_xxxz[k] = -g_z_0_y_0_xxy_xxxz[k] * ab_x + g_z_0_y_0_xxy_xxxxz[k];

                g_z_0_y_0_xxxy_xxyy[k] = -g_z_0_y_0_xxy_xxyy[k] * ab_x + g_z_0_y_0_xxy_xxxyy[k];

                g_z_0_y_0_xxxy_xxyz[k] = -g_z_0_y_0_xxy_xxyz[k] * ab_x + g_z_0_y_0_xxy_xxxyz[k];

                g_z_0_y_0_xxxy_xxzz[k] = -g_z_0_y_0_xxy_xxzz[k] * ab_x + g_z_0_y_0_xxy_xxxzz[k];

                g_z_0_y_0_xxxy_xyyy[k] = -g_z_0_y_0_xxy_xyyy[k] * ab_x + g_z_0_y_0_xxy_xxyyy[k];

                g_z_0_y_0_xxxy_xyyz[k] = -g_z_0_y_0_xxy_xyyz[k] * ab_x + g_z_0_y_0_xxy_xxyyz[k];

                g_z_0_y_0_xxxy_xyzz[k] = -g_z_0_y_0_xxy_xyzz[k] * ab_x + g_z_0_y_0_xxy_xxyzz[k];

                g_z_0_y_0_xxxy_xzzz[k] = -g_z_0_y_0_xxy_xzzz[k] * ab_x + g_z_0_y_0_xxy_xxzzz[k];

                g_z_0_y_0_xxxy_yyyy[k] = -g_z_0_y_0_xxy_yyyy[k] * ab_x + g_z_0_y_0_xxy_xyyyy[k];

                g_z_0_y_0_xxxy_yyyz[k] = -g_z_0_y_0_xxy_yyyz[k] * ab_x + g_z_0_y_0_xxy_xyyyz[k];

                g_z_0_y_0_xxxy_yyzz[k] = -g_z_0_y_0_xxy_yyzz[k] * ab_x + g_z_0_y_0_xxy_xyyzz[k];

                g_z_0_y_0_xxxy_yzzz[k] = -g_z_0_y_0_xxy_yzzz[k] * ab_x + g_z_0_y_0_xxy_xyzzz[k];

                g_z_0_y_0_xxxy_zzzz[k] = -g_z_0_y_0_xxy_zzzz[k] * ab_x + g_z_0_y_0_xxy_xzzzz[k];
            }

            /// Set up 1605-1620 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxz_xxxx = cbuffer.data(gg_geom_1010_off + 1605 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_xxxy = cbuffer.data(gg_geom_1010_off + 1606 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_xxxz = cbuffer.data(gg_geom_1010_off + 1607 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_xxyy = cbuffer.data(gg_geom_1010_off + 1608 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_xxyz = cbuffer.data(gg_geom_1010_off + 1609 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_xxzz = cbuffer.data(gg_geom_1010_off + 1610 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_xyyy = cbuffer.data(gg_geom_1010_off + 1611 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_xyyz = cbuffer.data(gg_geom_1010_off + 1612 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_xyzz = cbuffer.data(gg_geom_1010_off + 1613 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_xzzz = cbuffer.data(gg_geom_1010_off + 1614 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_yyyy = cbuffer.data(gg_geom_1010_off + 1615 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_yyyz = cbuffer.data(gg_geom_1010_off + 1616 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_yyzz = cbuffer.data(gg_geom_1010_off + 1617 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_yzzz = cbuffer.data(gg_geom_1010_off + 1618 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_zzzz = cbuffer.data(gg_geom_1010_off + 1619 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxz_xxxx, g_z_0_y_0_xxxz_xxxy, g_z_0_y_0_xxxz_xxxz, g_z_0_y_0_xxxz_xxyy, g_z_0_y_0_xxxz_xxyz, g_z_0_y_0_xxxz_xxzz, g_z_0_y_0_xxxz_xyyy, g_z_0_y_0_xxxz_xyyz, g_z_0_y_0_xxxz_xyzz, g_z_0_y_0_xxxz_xzzz, g_z_0_y_0_xxxz_yyyy, g_z_0_y_0_xxxz_yyyz, g_z_0_y_0_xxxz_yyzz, g_z_0_y_0_xxxz_yzzz, g_z_0_y_0_xxxz_zzzz, g_z_0_y_0_xxz_xxxx, g_z_0_y_0_xxz_xxxxx, g_z_0_y_0_xxz_xxxxy, g_z_0_y_0_xxz_xxxxz, g_z_0_y_0_xxz_xxxy, g_z_0_y_0_xxz_xxxyy, g_z_0_y_0_xxz_xxxyz, g_z_0_y_0_xxz_xxxz, g_z_0_y_0_xxz_xxxzz, g_z_0_y_0_xxz_xxyy, g_z_0_y_0_xxz_xxyyy, g_z_0_y_0_xxz_xxyyz, g_z_0_y_0_xxz_xxyz, g_z_0_y_0_xxz_xxyzz, g_z_0_y_0_xxz_xxzz, g_z_0_y_0_xxz_xxzzz, g_z_0_y_0_xxz_xyyy, g_z_0_y_0_xxz_xyyyy, g_z_0_y_0_xxz_xyyyz, g_z_0_y_0_xxz_xyyz, g_z_0_y_0_xxz_xyyzz, g_z_0_y_0_xxz_xyzz, g_z_0_y_0_xxz_xyzzz, g_z_0_y_0_xxz_xzzz, g_z_0_y_0_xxz_xzzzz, g_z_0_y_0_xxz_yyyy, g_z_0_y_0_xxz_yyyz, g_z_0_y_0_xxz_yyzz, g_z_0_y_0_xxz_yzzz, g_z_0_y_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxz_xxxx[k] = -g_z_0_y_0_xxz_xxxx[k] * ab_x + g_z_0_y_0_xxz_xxxxx[k];

                g_z_0_y_0_xxxz_xxxy[k] = -g_z_0_y_0_xxz_xxxy[k] * ab_x + g_z_0_y_0_xxz_xxxxy[k];

                g_z_0_y_0_xxxz_xxxz[k] = -g_z_0_y_0_xxz_xxxz[k] * ab_x + g_z_0_y_0_xxz_xxxxz[k];

                g_z_0_y_0_xxxz_xxyy[k] = -g_z_0_y_0_xxz_xxyy[k] * ab_x + g_z_0_y_0_xxz_xxxyy[k];

                g_z_0_y_0_xxxz_xxyz[k] = -g_z_0_y_0_xxz_xxyz[k] * ab_x + g_z_0_y_0_xxz_xxxyz[k];

                g_z_0_y_0_xxxz_xxzz[k] = -g_z_0_y_0_xxz_xxzz[k] * ab_x + g_z_0_y_0_xxz_xxxzz[k];

                g_z_0_y_0_xxxz_xyyy[k] = -g_z_0_y_0_xxz_xyyy[k] * ab_x + g_z_0_y_0_xxz_xxyyy[k];

                g_z_0_y_0_xxxz_xyyz[k] = -g_z_0_y_0_xxz_xyyz[k] * ab_x + g_z_0_y_0_xxz_xxyyz[k];

                g_z_0_y_0_xxxz_xyzz[k] = -g_z_0_y_0_xxz_xyzz[k] * ab_x + g_z_0_y_0_xxz_xxyzz[k];

                g_z_0_y_0_xxxz_xzzz[k] = -g_z_0_y_0_xxz_xzzz[k] * ab_x + g_z_0_y_0_xxz_xxzzz[k];

                g_z_0_y_0_xxxz_yyyy[k] = -g_z_0_y_0_xxz_yyyy[k] * ab_x + g_z_0_y_0_xxz_xyyyy[k];

                g_z_0_y_0_xxxz_yyyz[k] = -g_z_0_y_0_xxz_yyyz[k] * ab_x + g_z_0_y_0_xxz_xyyyz[k];

                g_z_0_y_0_xxxz_yyzz[k] = -g_z_0_y_0_xxz_yyzz[k] * ab_x + g_z_0_y_0_xxz_xyyzz[k];

                g_z_0_y_0_xxxz_yzzz[k] = -g_z_0_y_0_xxz_yzzz[k] * ab_x + g_z_0_y_0_xxz_xyzzz[k];

                g_z_0_y_0_xxxz_zzzz[k] = -g_z_0_y_0_xxz_zzzz[k] * ab_x + g_z_0_y_0_xxz_xzzzz[k];
            }

            /// Set up 1620-1635 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyy_xxxx = cbuffer.data(gg_geom_1010_off + 1620 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_xxxy = cbuffer.data(gg_geom_1010_off + 1621 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_xxxz = cbuffer.data(gg_geom_1010_off + 1622 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_xxyy = cbuffer.data(gg_geom_1010_off + 1623 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_xxyz = cbuffer.data(gg_geom_1010_off + 1624 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_xxzz = cbuffer.data(gg_geom_1010_off + 1625 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_xyyy = cbuffer.data(gg_geom_1010_off + 1626 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_xyyz = cbuffer.data(gg_geom_1010_off + 1627 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_xyzz = cbuffer.data(gg_geom_1010_off + 1628 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_xzzz = cbuffer.data(gg_geom_1010_off + 1629 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_yyyy = cbuffer.data(gg_geom_1010_off + 1630 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_yyyz = cbuffer.data(gg_geom_1010_off + 1631 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_yyzz = cbuffer.data(gg_geom_1010_off + 1632 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_yzzz = cbuffer.data(gg_geom_1010_off + 1633 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_zzzz = cbuffer.data(gg_geom_1010_off + 1634 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyy_xxxx, g_z_0_y_0_xxyy_xxxy, g_z_0_y_0_xxyy_xxxz, g_z_0_y_0_xxyy_xxyy, g_z_0_y_0_xxyy_xxyz, g_z_0_y_0_xxyy_xxzz, g_z_0_y_0_xxyy_xyyy, g_z_0_y_0_xxyy_xyyz, g_z_0_y_0_xxyy_xyzz, g_z_0_y_0_xxyy_xzzz, g_z_0_y_0_xxyy_yyyy, g_z_0_y_0_xxyy_yyyz, g_z_0_y_0_xxyy_yyzz, g_z_0_y_0_xxyy_yzzz, g_z_0_y_0_xxyy_zzzz, g_z_0_y_0_xyy_xxxx, g_z_0_y_0_xyy_xxxxx, g_z_0_y_0_xyy_xxxxy, g_z_0_y_0_xyy_xxxxz, g_z_0_y_0_xyy_xxxy, g_z_0_y_0_xyy_xxxyy, g_z_0_y_0_xyy_xxxyz, g_z_0_y_0_xyy_xxxz, g_z_0_y_0_xyy_xxxzz, g_z_0_y_0_xyy_xxyy, g_z_0_y_0_xyy_xxyyy, g_z_0_y_0_xyy_xxyyz, g_z_0_y_0_xyy_xxyz, g_z_0_y_0_xyy_xxyzz, g_z_0_y_0_xyy_xxzz, g_z_0_y_0_xyy_xxzzz, g_z_0_y_0_xyy_xyyy, g_z_0_y_0_xyy_xyyyy, g_z_0_y_0_xyy_xyyyz, g_z_0_y_0_xyy_xyyz, g_z_0_y_0_xyy_xyyzz, g_z_0_y_0_xyy_xyzz, g_z_0_y_0_xyy_xyzzz, g_z_0_y_0_xyy_xzzz, g_z_0_y_0_xyy_xzzzz, g_z_0_y_0_xyy_yyyy, g_z_0_y_0_xyy_yyyz, g_z_0_y_0_xyy_yyzz, g_z_0_y_0_xyy_yzzz, g_z_0_y_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyy_xxxx[k] = -g_z_0_y_0_xyy_xxxx[k] * ab_x + g_z_0_y_0_xyy_xxxxx[k];

                g_z_0_y_0_xxyy_xxxy[k] = -g_z_0_y_0_xyy_xxxy[k] * ab_x + g_z_0_y_0_xyy_xxxxy[k];

                g_z_0_y_0_xxyy_xxxz[k] = -g_z_0_y_0_xyy_xxxz[k] * ab_x + g_z_0_y_0_xyy_xxxxz[k];

                g_z_0_y_0_xxyy_xxyy[k] = -g_z_0_y_0_xyy_xxyy[k] * ab_x + g_z_0_y_0_xyy_xxxyy[k];

                g_z_0_y_0_xxyy_xxyz[k] = -g_z_0_y_0_xyy_xxyz[k] * ab_x + g_z_0_y_0_xyy_xxxyz[k];

                g_z_0_y_0_xxyy_xxzz[k] = -g_z_0_y_0_xyy_xxzz[k] * ab_x + g_z_0_y_0_xyy_xxxzz[k];

                g_z_0_y_0_xxyy_xyyy[k] = -g_z_0_y_0_xyy_xyyy[k] * ab_x + g_z_0_y_0_xyy_xxyyy[k];

                g_z_0_y_0_xxyy_xyyz[k] = -g_z_0_y_0_xyy_xyyz[k] * ab_x + g_z_0_y_0_xyy_xxyyz[k];

                g_z_0_y_0_xxyy_xyzz[k] = -g_z_0_y_0_xyy_xyzz[k] * ab_x + g_z_0_y_0_xyy_xxyzz[k];

                g_z_0_y_0_xxyy_xzzz[k] = -g_z_0_y_0_xyy_xzzz[k] * ab_x + g_z_0_y_0_xyy_xxzzz[k];

                g_z_0_y_0_xxyy_yyyy[k] = -g_z_0_y_0_xyy_yyyy[k] * ab_x + g_z_0_y_0_xyy_xyyyy[k];

                g_z_0_y_0_xxyy_yyyz[k] = -g_z_0_y_0_xyy_yyyz[k] * ab_x + g_z_0_y_0_xyy_xyyyz[k];

                g_z_0_y_0_xxyy_yyzz[k] = -g_z_0_y_0_xyy_yyzz[k] * ab_x + g_z_0_y_0_xyy_xyyzz[k];

                g_z_0_y_0_xxyy_yzzz[k] = -g_z_0_y_0_xyy_yzzz[k] * ab_x + g_z_0_y_0_xyy_xyzzz[k];

                g_z_0_y_0_xxyy_zzzz[k] = -g_z_0_y_0_xyy_zzzz[k] * ab_x + g_z_0_y_0_xyy_xzzzz[k];
            }

            /// Set up 1635-1650 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyz_xxxx = cbuffer.data(gg_geom_1010_off + 1635 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_xxxy = cbuffer.data(gg_geom_1010_off + 1636 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_xxxz = cbuffer.data(gg_geom_1010_off + 1637 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_xxyy = cbuffer.data(gg_geom_1010_off + 1638 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_xxyz = cbuffer.data(gg_geom_1010_off + 1639 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_xxzz = cbuffer.data(gg_geom_1010_off + 1640 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_xyyy = cbuffer.data(gg_geom_1010_off + 1641 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_xyyz = cbuffer.data(gg_geom_1010_off + 1642 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_xyzz = cbuffer.data(gg_geom_1010_off + 1643 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_xzzz = cbuffer.data(gg_geom_1010_off + 1644 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_yyyy = cbuffer.data(gg_geom_1010_off + 1645 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_yyyz = cbuffer.data(gg_geom_1010_off + 1646 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_yyzz = cbuffer.data(gg_geom_1010_off + 1647 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_yzzz = cbuffer.data(gg_geom_1010_off + 1648 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_zzzz = cbuffer.data(gg_geom_1010_off + 1649 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyz_xxxx, g_z_0_y_0_xxyz_xxxy, g_z_0_y_0_xxyz_xxxz, g_z_0_y_0_xxyz_xxyy, g_z_0_y_0_xxyz_xxyz, g_z_0_y_0_xxyz_xxzz, g_z_0_y_0_xxyz_xyyy, g_z_0_y_0_xxyz_xyyz, g_z_0_y_0_xxyz_xyzz, g_z_0_y_0_xxyz_xzzz, g_z_0_y_0_xxyz_yyyy, g_z_0_y_0_xxyz_yyyz, g_z_0_y_0_xxyz_yyzz, g_z_0_y_0_xxyz_yzzz, g_z_0_y_0_xxyz_zzzz, g_z_0_y_0_xyz_xxxx, g_z_0_y_0_xyz_xxxxx, g_z_0_y_0_xyz_xxxxy, g_z_0_y_0_xyz_xxxxz, g_z_0_y_0_xyz_xxxy, g_z_0_y_0_xyz_xxxyy, g_z_0_y_0_xyz_xxxyz, g_z_0_y_0_xyz_xxxz, g_z_0_y_0_xyz_xxxzz, g_z_0_y_0_xyz_xxyy, g_z_0_y_0_xyz_xxyyy, g_z_0_y_0_xyz_xxyyz, g_z_0_y_0_xyz_xxyz, g_z_0_y_0_xyz_xxyzz, g_z_0_y_0_xyz_xxzz, g_z_0_y_0_xyz_xxzzz, g_z_0_y_0_xyz_xyyy, g_z_0_y_0_xyz_xyyyy, g_z_0_y_0_xyz_xyyyz, g_z_0_y_0_xyz_xyyz, g_z_0_y_0_xyz_xyyzz, g_z_0_y_0_xyz_xyzz, g_z_0_y_0_xyz_xyzzz, g_z_0_y_0_xyz_xzzz, g_z_0_y_0_xyz_xzzzz, g_z_0_y_0_xyz_yyyy, g_z_0_y_0_xyz_yyyz, g_z_0_y_0_xyz_yyzz, g_z_0_y_0_xyz_yzzz, g_z_0_y_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyz_xxxx[k] = -g_z_0_y_0_xyz_xxxx[k] * ab_x + g_z_0_y_0_xyz_xxxxx[k];

                g_z_0_y_0_xxyz_xxxy[k] = -g_z_0_y_0_xyz_xxxy[k] * ab_x + g_z_0_y_0_xyz_xxxxy[k];

                g_z_0_y_0_xxyz_xxxz[k] = -g_z_0_y_0_xyz_xxxz[k] * ab_x + g_z_0_y_0_xyz_xxxxz[k];

                g_z_0_y_0_xxyz_xxyy[k] = -g_z_0_y_0_xyz_xxyy[k] * ab_x + g_z_0_y_0_xyz_xxxyy[k];

                g_z_0_y_0_xxyz_xxyz[k] = -g_z_0_y_0_xyz_xxyz[k] * ab_x + g_z_0_y_0_xyz_xxxyz[k];

                g_z_0_y_0_xxyz_xxzz[k] = -g_z_0_y_0_xyz_xxzz[k] * ab_x + g_z_0_y_0_xyz_xxxzz[k];

                g_z_0_y_0_xxyz_xyyy[k] = -g_z_0_y_0_xyz_xyyy[k] * ab_x + g_z_0_y_0_xyz_xxyyy[k];

                g_z_0_y_0_xxyz_xyyz[k] = -g_z_0_y_0_xyz_xyyz[k] * ab_x + g_z_0_y_0_xyz_xxyyz[k];

                g_z_0_y_0_xxyz_xyzz[k] = -g_z_0_y_0_xyz_xyzz[k] * ab_x + g_z_0_y_0_xyz_xxyzz[k];

                g_z_0_y_0_xxyz_xzzz[k] = -g_z_0_y_0_xyz_xzzz[k] * ab_x + g_z_0_y_0_xyz_xxzzz[k];

                g_z_0_y_0_xxyz_yyyy[k] = -g_z_0_y_0_xyz_yyyy[k] * ab_x + g_z_0_y_0_xyz_xyyyy[k];

                g_z_0_y_0_xxyz_yyyz[k] = -g_z_0_y_0_xyz_yyyz[k] * ab_x + g_z_0_y_0_xyz_xyyyz[k];

                g_z_0_y_0_xxyz_yyzz[k] = -g_z_0_y_0_xyz_yyzz[k] * ab_x + g_z_0_y_0_xyz_xyyzz[k];

                g_z_0_y_0_xxyz_yzzz[k] = -g_z_0_y_0_xyz_yzzz[k] * ab_x + g_z_0_y_0_xyz_xyzzz[k];

                g_z_0_y_0_xxyz_zzzz[k] = -g_z_0_y_0_xyz_zzzz[k] * ab_x + g_z_0_y_0_xyz_xzzzz[k];
            }

            /// Set up 1650-1665 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxzz_xxxx = cbuffer.data(gg_geom_1010_off + 1650 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_xxxy = cbuffer.data(gg_geom_1010_off + 1651 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_xxxz = cbuffer.data(gg_geom_1010_off + 1652 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_xxyy = cbuffer.data(gg_geom_1010_off + 1653 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_xxyz = cbuffer.data(gg_geom_1010_off + 1654 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_xxzz = cbuffer.data(gg_geom_1010_off + 1655 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_xyyy = cbuffer.data(gg_geom_1010_off + 1656 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_xyyz = cbuffer.data(gg_geom_1010_off + 1657 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_xyzz = cbuffer.data(gg_geom_1010_off + 1658 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_xzzz = cbuffer.data(gg_geom_1010_off + 1659 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_yyyy = cbuffer.data(gg_geom_1010_off + 1660 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_yyyz = cbuffer.data(gg_geom_1010_off + 1661 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_yyzz = cbuffer.data(gg_geom_1010_off + 1662 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_yzzz = cbuffer.data(gg_geom_1010_off + 1663 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_zzzz = cbuffer.data(gg_geom_1010_off + 1664 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxzz_xxxx, g_z_0_y_0_xxzz_xxxy, g_z_0_y_0_xxzz_xxxz, g_z_0_y_0_xxzz_xxyy, g_z_0_y_0_xxzz_xxyz, g_z_0_y_0_xxzz_xxzz, g_z_0_y_0_xxzz_xyyy, g_z_0_y_0_xxzz_xyyz, g_z_0_y_0_xxzz_xyzz, g_z_0_y_0_xxzz_xzzz, g_z_0_y_0_xxzz_yyyy, g_z_0_y_0_xxzz_yyyz, g_z_0_y_0_xxzz_yyzz, g_z_0_y_0_xxzz_yzzz, g_z_0_y_0_xxzz_zzzz, g_z_0_y_0_xzz_xxxx, g_z_0_y_0_xzz_xxxxx, g_z_0_y_0_xzz_xxxxy, g_z_0_y_0_xzz_xxxxz, g_z_0_y_0_xzz_xxxy, g_z_0_y_0_xzz_xxxyy, g_z_0_y_0_xzz_xxxyz, g_z_0_y_0_xzz_xxxz, g_z_0_y_0_xzz_xxxzz, g_z_0_y_0_xzz_xxyy, g_z_0_y_0_xzz_xxyyy, g_z_0_y_0_xzz_xxyyz, g_z_0_y_0_xzz_xxyz, g_z_0_y_0_xzz_xxyzz, g_z_0_y_0_xzz_xxzz, g_z_0_y_0_xzz_xxzzz, g_z_0_y_0_xzz_xyyy, g_z_0_y_0_xzz_xyyyy, g_z_0_y_0_xzz_xyyyz, g_z_0_y_0_xzz_xyyz, g_z_0_y_0_xzz_xyyzz, g_z_0_y_0_xzz_xyzz, g_z_0_y_0_xzz_xyzzz, g_z_0_y_0_xzz_xzzz, g_z_0_y_0_xzz_xzzzz, g_z_0_y_0_xzz_yyyy, g_z_0_y_0_xzz_yyyz, g_z_0_y_0_xzz_yyzz, g_z_0_y_0_xzz_yzzz, g_z_0_y_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxzz_xxxx[k] = -g_z_0_y_0_xzz_xxxx[k] * ab_x + g_z_0_y_0_xzz_xxxxx[k];

                g_z_0_y_0_xxzz_xxxy[k] = -g_z_0_y_0_xzz_xxxy[k] * ab_x + g_z_0_y_0_xzz_xxxxy[k];

                g_z_0_y_0_xxzz_xxxz[k] = -g_z_0_y_0_xzz_xxxz[k] * ab_x + g_z_0_y_0_xzz_xxxxz[k];

                g_z_0_y_0_xxzz_xxyy[k] = -g_z_0_y_0_xzz_xxyy[k] * ab_x + g_z_0_y_0_xzz_xxxyy[k];

                g_z_0_y_0_xxzz_xxyz[k] = -g_z_0_y_0_xzz_xxyz[k] * ab_x + g_z_0_y_0_xzz_xxxyz[k];

                g_z_0_y_0_xxzz_xxzz[k] = -g_z_0_y_0_xzz_xxzz[k] * ab_x + g_z_0_y_0_xzz_xxxzz[k];

                g_z_0_y_0_xxzz_xyyy[k] = -g_z_0_y_0_xzz_xyyy[k] * ab_x + g_z_0_y_0_xzz_xxyyy[k];

                g_z_0_y_0_xxzz_xyyz[k] = -g_z_0_y_0_xzz_xyyz[k] * ab_x + g_z_0_y_0_xzz_xxyyz[k];

                g_z_0_y_0_xxzz_xyzz[k] = -g_z_0_y_0_xzz_xyzz[k] * ab_x + g_z_0_y_0_xzz_xxyzz[k];

                g_z_0_y_0_xxzz_xzzz[k] = -g_z_0_y_0_xzz_xzzz[k] * ab_x + g_z_0_y_0_xzz_xxzzz[k];

                g_z_0_y_0_xxzz_yyyy[k] = -g_z_0_y_0_xzz_yyyy[k] * ab_x + g_z_0_y_0_xzz_xyyyy[k];

                g_z_0_y_0_xxzz_yyyz[k] = -g_z_0_y_0_xzz_yyyz[k] * ab_x + g_z_0_y_0_xzz_xyyyz[k];

                g_z_0_y_0_xxzz_yyzz[k] = -g_z_0_y_0_xzz_yyzz[k] * ab_x + g_z_0_y_0_xzz_xyyzz[k];

                g_z_0_y_0_xxzz_yzzz[k] = -g_z_0_y_0_xzz_yzzz[k] * ab_x + g_z_0_y_0_xzz_xyzzz[k];

                g_z_0_y_0_xxzz_zzzz[k] = -g_z_0_y_0_xzz_zzzz[k] * ab_x + g_z_0_y_0_xzz_xzzzz[k];
            }

            /// Set up 1665-1680 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyy_xxxx = cbuffer.data(gg_geom_1010_off + 1665 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_xxxy = cbuffer.data(gg_geom_1010_off + 1666 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_xxxz = cbuffer.data(gg_geom_1010_off + 1667 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_xxyy = cbuffer.data(gg_geom_1010_off + 1668 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_xxyz = cbuffer.data(gg_geom_1010_off + 1669 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_xxzz = cbuffer.data(gg_geom_1010_off + 1670 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_xyyy = cbuffer.data(gg_geom_1010_off + 1671 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_xyyz = cbuffer.data(gg_geom_1010_off + 1672 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_xyzz = cbuffer.data(gg_geom_1010_off + 1673 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_xzzz = cbuffer.data(gg_geom_1010_off + 1674 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1675 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1676 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1677 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1678 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1679 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyy_xxxx, g_z_0_y_0_xyyy_xxxy, g_z_0_y_0_xyyy_xxxz, g_z_0_y_0_xyyy_xxyy, g_z_0_y_0_xyyy_xxyz, g_z_0_y_0_xyyy_xxzz, g_z_0_y_0_xyyy_xyyy, g_z_0_y_0_xyyy_xyyz, g_z_0_y_0_xyyy_xyzz, g_z_0_y_0_xyyy_xzzz, g_z_0_y_0_xyyy_yyyy, g_z_0_y_0_xyyy_yyyz, g_z_0_y_0_xyyy_yyzz, g_z_0_y_0_xyyy_yzzz, g_z_0_y_0_xyyy_zzzz, g_z_0_y_0_yyy_xxxx, g_z_0_y_0_yyy_xxxxx, g_z_0_y_0_yyy_xxxxy, g_z_0_y_0_yyy_xxxxz, g_z_0_y_0_yyy_xxxy, g_z_0_y_0_yyy_xxxyy, g_z_0_y_0_yyy_xxxyz, g_z_0_y_0_yyy_xxxz, g_z_0_y_0_yyy_xxxzz, g_z_0_y_0_yyy_xxyy, g_z_0_y_0_yyy_xxyyy, g_z_0_y_0_yyy_xxyyz, g_z_0_y_0_yyy_xxyz, g_z_0_y_0_yyy_xxyzz, g_z_0_y_0_yyy_xxzz, g_z_0_y_0_yyy_xxzzz, g_z_0_y_0_yyy_xyyy, g_z_0_y_0_yyy_xyyyy, g_z_0_y_0_yyy_xyyyz, g_z_0_y_0_yyy_xyyz, g_z_0_y_0_yyy_xyyzz, g_z_0_y_0_yyy_xyzz, g_z_0_y_0_yyy_xyzzz, g_z_0_y_0_yyy_xzzz, g_z_0_y_0_yyy_xzzzz, g_z_0_y_0_yyy_yyyy, g_z_0_y_0_yyy_yyyz, g_z_0_y_0_yyy_yyzz, g_z_0_y_0_yyy_yzzz, g_z_0_y_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyy_xxxx[k] = -g_z_0_y_0_yyy_xxxx[k] * ab_x + g_z_0_y_0_yyy_xxxxx[k];

                g_z_0_y_0_xyyy_xxxy[k] = -g_z_0_y_0_yyy_xxxy[k] * ab_x + g_z_0_y_0_yyy_xxxxy[k];

                g_z_0_y_0_xyyy_xxxz[k] = -g_z_0_y_0_yyy_xxxz[k] * ab_x + g_z_0_y_0_yyy_xxxxz[k];

                g_z_0_y_0_xyyy_xxyy[k] = -g_z_0_y_0_yyy_xxyy[k] * ab_x + g_z_0_y_0_yyy_xxxyy[k];

                g_z_0_y_0_xyyy_xxyz[k] = -g_z_0_y_0_yyy_xxyz[k] * ab_x + g_z_0_y_0_yyy_xxxyz[k];

                g_z_0_y_0_xyyy_xxzz[k] = -g_z_0_y_0_yyy_xxzz[k] * ab_x + g_z_0_y_0_yyy_xxxzz[k];

                g_z_0_y_0_xyyy_xyyy[k] = -g_z_0_y_0_yyy_xyyy[k] * ab_x + g_z_0_y_0_yyy_xxyyy[k];

                g_z_0_y_0_xyyy_xyyz[k] = -g_z_0_y_0_yyy_xyyz[k] * ab_x + g_z_0_y_0_yyy_xxyyz[k];

                g_z_0_y_0_xyyy_xyzz[k] = -g_z_0_y_0_yyy_xyzz[k] * ab_x + g_z_0_y_0_yyy_xxyzz[k];

                g_z_0_y_0_xyyy_xzzz[k] = -g_z_0_y_0_yyy_xzzz[k] * ab_x + g_z_0_y_0_yyy_xxzzz[k];

                g_z_0_y_0_xyyy_yyyy[k] = -g_z_0_y_0_yyy_yyyy[k] * ab_x + g_z_0_y_0_yyy_xyyyy[k];

                g_z_0_y_0_xyyy_yyyz[k] = -g_z_0_y_0_yyy_yyyz[k] * ab_x + g_z_0_y_0_yyy_xyyyz[k];

                g_z_0_y_0_xyyy_yyzz[k] = -g_z_0_y_0_yyy_yyzz[k] * ab_x + g_z_0_y_0_yyy_xyyzz[k];

                g_z_0_y_0_xyyy_yzzz[k] = -g_z_0_y_0_yyy_yzzz[k] * ab_x + g_z_0_y_0_yyy_xyzzz[k];

                g_z_0_y_0_xyyy_zzzz[k] = -g_z_0_y_0_yyy_zzzz[k] * ab_x + g_z_0_y_0_yyy_xzzzz[k];
            }

            /// Set up 1680-1695 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1680 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1681 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1682 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1683 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1684 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1685 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1686 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1687 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1688 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1689 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1690 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1691 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1692 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1693 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1694 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyz_xxxx, g_z_0_y_0_xyyz_xxxy, g_z_0_y_0_xyyz_xxxz, g_z_0_y_0_xyyz_xxyy, g_z_0_y_0_xyyz_xxyz, g_z_0_y_0_xyyz_xxzz, g_z_0_y_0_xyyz_xyyy, g_z_0_y_0_xyyz_xyyz, g_z_0_y_0_xyyz_xyzz, g_z_0_y_0_xyyz_xzzz, g_z_0_y_0_xyyz_yyyy, g_z_0_y_0_xyyz_yyyz, g_z_0_y_0_xyyz_yyzz, g_z_0_y_0_xyyz_yzzz, g_z_0_y_0_xyyz_zzzz, g_z_0_y_0_yyz_xxxx, g_z_0_y_0_yyz_xxxxx, g_z_0_y_0_yyz_xxxxy, g_z_0_y_0_yyz_xxxxz, g_z_0_y_0_yyz_xxxy, g_z_0_y_0_yyz_xxxyy, g_z_0_y_0_yyz_xxxyz, g_z_0_y_0_yyz_xxxz, g_z_0_y_0_yyz_xxxzz, g_z_0_y_0_yyz_xxyy, g_z_0_y_0_yyz_xxyyy, g_z_0_y_0_yyz_xxyyz, g_z_0_y_0_yyz_xxyz, g_z_0_y_0_yyz_xxyzz, g_z_0_y_0_yyz_xxzz, g_z_0_y_0_yyz_xxzzz, g_z_0_y_0_yyz_xyyy, g_z_0_y_0_yyz_xyyyy, g_z_0_y_0_yyz_xyyyz, g_z_0_y_0_yyz_xyyz, g_z_0_y_0_yyz_xyyzz, g_z_0_y_0_yyz_xyzz, g_z_0_y_0_yyz_xyzzz, g_z_0_y_0_yyz_xzzz, g_z_0_y_0_yyz_xzzzz, g_z_0_y_0_yyz_yyyy, g_z_0_y_0_yyz_yyyz, g_z_0_y_0_yyz_yyzz, g_z_0_y_0_yyz_yzzz, g_z_0_y_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyz_xxxx[k] = -g_z_0_y_0_yyz_xxxx[k] * ab_x + g_z_0_y_0_yyz_xxxxx[k];

                g_z_0_y_0_xyyz_xxxy[k] = -g_z_0_y_0_yyz_xxxy[k] * ab_x + g_z_0_y_0_yyz_xxxxy[k];

                g_z_0_y_0_xyyz_xxxz[k] = -g_z_0_y_0_yyz_xxxz[k] * ab_x + g_z_0_y_0_yyz_xxxxz[k];

                g_z_0_y_0_xyyz_xxyy[k] = -g_z_0_y_0_yyz_xxyy[k] * ab_x + g_z_0_y_0_yyz_xxxyy[k];

                g_z_0_y_0_xyyz_xxyz[k] = -g_z_0_y_0_yyz_xxyz[k] * ab_x + g_z_0_y_0_yyz_xxxyz[k];

                g_z_0_y_0_xyyz_xxzz[k] = -g_z_0_y_0_yyz_xxzz[k] * ab_x + g_z_0_y_0_yyz_xxxzz[k];

                g_z_0_y_0_xyyz_xyyy[k] = -g_z_0_y_0_yyz_xyyy[k] * ab_x + g_z_0_y_0_yyz_xxyyy[k];

                g_z_0_y_0_xyyz_xyyz[k] = -g_z_0_y_0_yyz_xyyz[k] * ab_x + g_z_0_y_0_yyz_xxyyz[k];

                g_z_0_y_0_xyyz_xyzz[k] = -g_z_0_y_0_yyz_xyzz[k] * ab_x + g_z_0_y_0_yyz_xxyzz[k];

                g_z_0_y_0_xyyz_xzzz[k] = -g_z_0_y_0_yyz_xzzz[k] * ab_x + g_z_0_y_0_yyz_xxzzz[k];

                g_z_0_y_0_xyyz_yyyy[k] = -g_z_0_y_0_yyz_yyyy[k] * ab_x + g_z_0_y_0_yyz_xyyyy[k];

                g_z_0_y_0_xyyz_yyyz[k] = -g_z_0_y_0_yyz_yyyz[k] * ab_x + g_z_0_y_0_yyz_xyyyz[k];

                g_z_0_y_0_xyyz_yyzz[k] = -g_z_0_y_0_yyz_yyzz[k] * ab_x + g_z_0_y_0_yyz_xyyzz[k];

                g_z_0_y_0_xyyz_yzzz[k] = -g_z_0_y_0_yyz_yzzz[k] * ab_x + g_z_0_y_0_yyz_xyzzz[k];

                g_z_0_y_0_xyyz_zzzz[k] = -g_z_0_y_0_yyz_zzzz[k] * ab_x + g_z_0_y_0_yyz_xzzzz[k];
            }

            /// Set up 1695-1710 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1695 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1696 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1697 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1698 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1699 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1700 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1701 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1702 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1703 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1704 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1705 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1706 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1707 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1708 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1709 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyzz_xxxx, g_z_0_y_0_xyzz_xxxy, g_z_0_y_0_xyzz_xxxz, g_z_0_y_0_xyzz_xxyy, g_z_0_y_0_xyzz_xxyz, g_z_0_y_0_xyzz_xxzz, g_z_0_y_0_xyzz_xyyy, g_z_0_y_0_xyzz_xyyz, g_z_0_y_0_xyzz_xyzz, g_z_0_y_0_xyzz_xzzz, g_z_0_y_0_xyzz_yyyy, g_z_0_y_0_xyzz_yyyz, g_z_0_y_0_xyzz_yyzz, g_z_0_y_0_xyzz_yzzz, g_z_0_y_0_xyzz_zzzz, g_z_0_y_0_yzz_xxxx, g_z_0_y_0_yzz_xxxxx, g_z_0_y_0_yzz_xxxxy, g_z_0_y_0_yzz_xxxxz, g_z_0_y_0_yzz_xxxy, g_z_0_y_0_yzz_xxxyy, g_z_0_y_0_yzz_xxxyz, g_z_0_y_0_yzz_xxxz, g_z_0_y_0_yzz_xxxzz, g_z_0_y_0_yzz_xxyy, g_z_0_y_0_yzz_xxyyy, g_z_0_y_0_yzz_xxyyz, g_z_0_y_0_yzz_xxyz, g_z_0_y_0_yzz_xxyzz, g_z_0_y_0_yzz_xxzz, g_z_0_y_0_yzz_xxzzz, g_z_0_y_0_yzz_xyyy, g_z_0_y_0_yzz_xyyyy, g_z_0_y_0_yzz_xyyyz, g_z_0_y_0_yzz_xyyz, g_z_0_y_0_yzz_xyyzz, g_z_0_y_0_yzz_xyzz, g_z_0_y_0_yzz_xyzzz, g_z_0_y_0_yzz_xzzz, g_z_0_y_0_yzz_xzzzz, g_z_0_y_0_yzz_yyyy, g_z_0_y_0_yzz_yyyz, g_z_0_y_0_yzz_yyzz, g_z_0_y_0_yzz_yzzz, g_z_0_y_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyzz_xxxx[k] = -g_z_0_y_0_yzz_xxxx[k] * ab_x + g_z_0_y_0_yzz_xxxxx[k];

                g_z_0_y_0_xyzz_xxxy[k] = -g_z_0_y_0_yzz_xxxy[k] * ab_x + g_z_0_y_0_yzz_xxxxy[k];

                g_z_0_y_0_xyzz_xxxz[k] = -g_z_0_y_0_yzz_xxxz[k] * ab_x + g_z_0_y_0_yzz_xxxxz[k];

                g_z_0_y_0_xyzz_xxyy[k] = -g_z_0_y_0_yzz_xxyy[k] * ab_x + g_z_0_y_0_yzz_xxxyy[k];

                g_z_0_y_0_xyzz_xxyz[k] = -g_z_0_y_0_yzz_xxyz[k] * ab_x + g_z_0_y_0_yzz_xxxyz[k];

                g_z_0_y_0_xyzz_xxzz[k] = -g_z_0_y_0_yzz_xxzz[k] * ab_x + g_z_0_y_0_yzz_xxxzz[k];

                g_z_0_y_0_xyzz_xyyy[k] = -g_z_0_y_0_yzz_xyyy[k] * ab_x + g_z_0_y_0_yzz_xxyyy[k];

                g_z_0_y_0_xyzz_xyyz[k] = -g_z_0_y_0_yzz_xyyz[k] * ab_x + g_z_0_y_0_yzz_xxyyz[k];

                g_z_0_y_0_xyzz_xyzz[k] = -g_z_0_y_0_yzz_xyzz[k] * ab_x + g_z_0_y_0_yzz_xxyzz[k];

                g_z_0_y_0_xyzz_xzzz[k] = -g_z_0_y_0_yzz_xzzz[k] * ab_x + g_z_0_y_0_yzz_xxzzz[k];

                g_z_0_y_0_xyzz_yyyy[k] = -g_z_0_y_0_yzz_yyyy[k] * ab_x + g_z_0_y_0_yzz_xyyyy[k];

                g_z_0_y_0_xyzz_yyyz[k] = -g_z_0_y_0_yzz_yyyz[k] * ab_x + g_z_0_y_0_yzz_xyyyz[k];

                g_z_0_y_0_xyzz_yyzz[k] = -g_z_0_y_0_yzz_yyzz[k] * ab_x + g_z_0_y_0_yzz_xyyzz[k];

                g_z_0_y_0_xyzz_yzzz[k] = -g_z_0_y_0_yzz_yzzz[k] * ab_x + g_z_0_y_0_yzz_xyzzz[k];

                g_z_0_y_0_xyzz_zzzz[k] = -g_z_0_y_0_yzz_zzzz[k] * ab_x + g_z_0_y_0_yzz_xzzzz[k];
            }

            /// Set up 1710-1725 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1710 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1711 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1712 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1713 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1714 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1715 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1716 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1717 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1718 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1719 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1720 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1721 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1722 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1723 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1724 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xzzz_xxxx, g_z_0_y_0_xzzz_xxxy, g_z_0_y_0_xzzz_xxxz, g_z_0_y_0_xzzz_xxyy, g_z_0_y_0_xzzz_xxyz, g_z_0_y_0_xzzz_xxzz, g_z_0_y_0_xzzz_xyyy, g_z_0_y_0_xzzz_xyyz, g_z_0_y_0_xzzz_xyzz, g_z_0_y_0_xzzz_xzzz, g_z_0_y_0_xzzz_yyyy, g_z_0_y_0_xzzz_yyyz, g_z_0_y_0_xzzz_yyzz, g_z_0_y_0_xzzz_yzzz, g_z_0_y_0_xzzz_zzzz, g_z_0_y_0_zzz_xxxx, g_z_0_y_0_zzz_xxxxx, g_z_0_y_0_zzz_xxxxy, g_z_0_y_0_zzz_xxxxz, g_z_0_y_0_zzz_xxxy, g_z_0_y_0_zzz_xxxyy, g_z_0_y_0_zzz_xxxyz, g_z_0_y_0_zzz_xxxz, g_z_0_y_0_zzz_xxxzz, g_z_0_y_0_zzz_xxyy, g_z_0_y_0_zzz_xxyyy, g_z_0_y_0_zzz_xxyyz, g_z_0_y_0_zzz_xxyz, g_z_0_y_0_zzz_xxyzz, g_z_0_y_0_zzz_xxzz, g_z_0_y_0_zzz_xxzzz, g_z_0_y_0_zzz_xyyy, g_z_0_y_0_zzz_xyyyy, g_z_0_y_0_zzz_xyyyz, g_z_0_y_0_zzz_xyyz, g_z_0_y_0_zzz_xyyzz, g_z_0_y_0_zzz_xyzz, g_z_0_y_0_zzz_xyzzz, g_z_0_y_0_zzz_xzzz, g_z_0_y_0_zzz_xzzzz, g_z_0_y_0_zzz_yyyy, g_z_0_y_0_zzz_yyyz, g_z_0_y_0_zzz_yyzz, g_z_0_y_0_zzz_yzzz, g_z_0_y_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xzzz_xxxx[k] = -g_z_0_y_0_zzz_xxxx[k] * ab_x + g_z_0_y_0_zzz_xxxxx[k];

                g_z_0_y_0_xzzz_xxxy[k] = -g_z_0_y_0_zzz_xxxy[k] * ab_x + g_z_0_y_0_zzz_xxxxy[k];

                g_z_0_y_0_xzzz_xxxz[k] = -g_z_0_y_0_zzz_xxxz[k] * ab_x + g_z_0_y_0_zzz_xxxxz[k];

                g_z_0_y_0_xzzz_xxyy[k] = -g_z_0_y_0_zzz_xxyy[k] * ab_x + g_z_0_y_0_zzz_xxxyy[k];

                g_z_0_y_0_xzzz_xxyz[k] = -g_z_0_y_0_zzz_xxyz[k] * ab_x + g_z_0_y_0_zzz_xxxyz[k];

                g_z_0_y_0_xzzz_xxzz[k] = -g_z_0_y_0_zzz_xxzz[k] * ab_x + g_z_0_y_0_zzz_xxxzz[k];

                g_z_0_y_0_xzzz_xyyy[k] = -g_z_0_y_0_zzz_xyyy[k] * ab_x + g_z_0_y_0_zzz_xxyyy[k];

                g_z_0_y_0_xzzz_xyyz[k] = -g_z_0_y_0_zzz_xyyz[k] * ab_x + g_z_0_y_0_zzz_xxyyz[k];

                g_z_0_y_0_xzzz_xyzz[k] = -g_z_0_y_0_zzz_xyzz[k] * ab_x + g_z_0_y_0_zzz_xxyzz[k];

                g_z_0_y_0_xzzz_xzzz[k] = -g_z_0_y_0_zzz_xzzz[k] * ab_x + g_z_0_y_0_zzz_xxzzz[k];

                g_z_0_y_0_xzzz_yyyy[k] = -g_z_0_y_0_zzz_yyyy[k] * ab_x + g_z_0_y_0_zzz_xyyyy[k];

                g_z_0_y_0_xzzz_yyyz[k] = -g_z_0_y_0_zzz_yyyz[k] * ab_x + g_z_0_y_0_zzz_xyyyz[k];

                g_z_0_y_0_xzzz_yyzz[k] = -g_z_0_y_0_zzz_yyzz[k] * ab_x + g_z_0_y_0_zzz_xyyzz[k];

                g_z_0_y_0_xzzz_yzzz[k] = -g_z_0_y_0_zzz_yzzz[k] * ab_x + g_z_0_y_0_zzz_xyzzz[k];

                g_z_0_y_0_xzzz_zzzz[k] = -g_z_0_y_0_zzz_zzzz[k] * ab_x + g_z_0_y_0_zzz_xzzzz[k];
            }

            /// Set up 1725-1740 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyy_xxxx = cbuffer.data(gg_geom_1010_off + 1725 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_xxxy = cbuffer.data(gg_geom_1010_off + 1726 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_xxxz = cbuffer.data(gg_geom_1010_off + 1727 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_xxyy = cbuffer.data(gg_geom_1010_off + 1728 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_xxyz = cbuffer.data(gg_geom_1010_off + 1729 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_xxzz = cbuffer.data(gg_geom_1010_off + 1730 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_xyyy = cbuffer.data(gg_geom_1010_off + 1731 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_xyyz = cbuffer.data(gg_geom_1010_off + 1732 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_xyzz = cbuffer.data(gg_geom_1010_off + 1733 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_xzzz = cbuffer.data(gg_geom_1010_off + 1734 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1735 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1736 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1737 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1738 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1739 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyy_xxxx, g_z_0_y_0_yyy_xxxxy, g_z_0_y_0_yyy_xxxy, g_z_0_y_0_yyy_xxxyy, g_z_0_y_0_yyy_xxxyz, g_z_0_y_0_yyy_xxxz, g_z_0_y_0_yyy_xxyy, g_z_0_y_0_yyy_xxyyy, g_z_0_y_0_yyy_xxyyz, g_z_0_y_0_yyy_xxyz, g_z_0_y_0_yyy_xxyzz, g_z_0_y_0_yyy_xxzz, g_z_0_y_0_yyy_xyyy, g_z_0_y_0_yyy_xyyyy, g_z_0_y_0_yyy_xyyyz, g_z_0_y_0_yyy_xyyz, g_z_0_y_0_yyy_xyyzz, g_z_0_y_0_yyy_xyzz, g_z_0_y_0_yyy_xyzzz, g_z_0_y_0_yyy_xzzz, g_z_0_y_0_yyy_yyyy, g_z_0_y_0_yyy_yyyyy, g_z_0_y_0_yyy_yyyyz, g_z_0_y_0_yyy_yyyz, g_z_0_y_0_yyy_yyyzz, g_z_0_y_0_yyy_yyzz, g_z_0_y_0_yyy_yyzzz, g_z_0_y_0_yyy_yzzz, g_z_0_y_0_yyy_yzzzz, g_z_0_y_0_yyy_zzzz, g_z_0_y_0_yyyy_xxxx, g_z_0_y_0_yyyy_xxxy, g_z_0_y_0_yyyy_xxxz, g_z_0_y_0_yyyy_xxyy, g_z_0_y_0_yyyy_xxyz, g_z_0_y_0_yyyy_xxzz, g_z_0_y_0_yyyy_xyyy, g_z_0_y_0_yyyy_xyyz, g_z_0_y_0_yyyy_xyzz, g_z_0_y_0_yyyy_xzzz, g_z_0_y_0_yyyy_yyyy, g_z_0_y_0_yyyy_yyyz, g_z_0_y_0_yyyy_yyzz, g_z_0_y_0_yyyy_yzzz, g_z_0_y_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyy_xxxx[k] = -g_z_0_y_0_yyy_xxxx[k] * ab_y + g_z_0_y_0_yyy_xxxxy[k];

                g_z_0_y_0_yyyy_xxxy[k] = -g_z_0_y_0_yyy_xxxy[k] * ab_y + g_z_0_y_0_yyy_xxxyy[k];

                g_z_0_y_0_yyyy_xxxz[k] = -g_z_0_y_0_yyy_xxxz[k] * ab_y + g_z_0_y_0_yyy_xxxyz[k];

                g_z_0_y_0_yyyy_xxyy[k] = -g_z_0_y_0_yyy_xxyy[k] * ab_y + g_z_0_y_0_yyy_xxyyy[k];

                g_z_0_y_0_yyyy_xxyz[k] = -g_z_0_y_0_yyy_xxyz[k] * ab_y + g_z_0_y_0_yyy_xxyyz[k];

                g_z_0_y_0_yyyy_xxzz[k] = -g_z_0_y_0_yyy_xxzz[k] * ab_y + g_z_0_y_0_yyy_xxyzz[k];

                g_z_0_y_0_yyyy_xyyy[k] = -g_z_0_y_0_yyy_xyyy[k] * ab_y + g_z_0_y_0_yyy_xyyyy[k];

                g_z_0_y_0_yyyy_xyyz[k] = -g_z_0_y_0_yyy_xyyz[k] * ab_y + g_z_0_y_0_yyy_xyyyz[k];

                g_z_0_y_0_yyyy_xyzz[k] = -g_z_0_y_0_yyy_xyzz[k] * ab_y + g_z_0_y_0_yyy_xyyzz[k];

                g_z_0_y_0_yyyy_xzzz[k] = -g_z_0_y_0_yyy_xzzz[k] * ab_y + g_z_0_y_0_yyy_xyzzz[k];

                g_z_0_y_0_yyyy_yyyy[k] = -g_z_0_y_0_yyy_yyyy[k] * ab_y + g_z_0_y_0_yyy_yyyyy[k];

                g_z_0_y_0_yyyy_yyyz[k] = -g_z_0_y_0_yyy_yyyz[k] * ab_y + g_z_0_y_0_yyy_yyyyz[k];

                g_z_0_y_0_yyyy_yyzz[k] = -g_z_0_y_0_yyy_yyzz[k] * ab_y + g_z_0_y_0_yyy_yyyzz[k];

                g_z_0_y_0_yyyy_yzzz[k] = -g_z_0_y_0_yyy_yzzz[k] * ab_y + g_z_0_y_0_yyy_yyzzz[k];

                g_z_0_y_0_yyyy_zzzz[k] = -g_z_0_y_0_yyy_zzzz[k] * ab_y + g_z_0_y_0_yyy_yzzzz[k];
            }

            /// Set up 1740-1755 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1740 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1741 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1742 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1743 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1744 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1745 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1746 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1747 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1748 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1749 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1750 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1751 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1752 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1753 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1754 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyz_xxxx, g_z_0_y_0_yyyz_xxxy, g_z_0_y_0_yyyz_xxxz, g_z_0_y_0_yyyz_xxyy, g_z_0_y_0_yyyz_xxyz, g_z_0_y_0_yyyz_xxzz, g_z_0_y_0_yyyz_xyyy, g_z_0_y_0_yyyz_xyyz, g_z_0_y_0_yyyz_xyzz, g_z_0_y_0_yyyz_xzzz, g_z_0_y_0_yyyz_yyyy, g_z_0_y_0_yyyz_yyyz, g_z_0_y_0_yyyz_yyzz, g_z_0_y_0_yyyz_yzzz, g_z_0_y_0_yyyz_zzzz, g_z_0_y_0_yyz_xxxx, g_z_0_y_0_yyz_xxxxy, g_z_0_y_0_yyz_xxxy, g_z_0_y_0_yyz_xxxyy, g_z_0_y_0_yyz_xxxyz, g_z_0_y_0_yyz_xxxz, g_z_0_y_0_yyz_xxyy, g_z_0_y_0_yyz_xxyyy, g_z_0_y_0_yyz_xxyyz, g_z_0_y_0_yyz_xxyz, g_z_0_y_0_yyz_xxyzz, g_z_0_y_0_yyz_xxzz, g_z_0_y_0_yyz_xyyy, g_z_0_y_0_yyz_xyyyy, g_z_0_y_0_yyz_xyyyz, g_z_0_y_0_yyz_xyyz, g_z_0_y_0_yyz_xyyzz, g_z_0_y_0_yyz_xyzz, g_z_0_y_0_yyz_xyzzz, g_z_0_y_0_yyz_xzzz, g_z_0_y_0_yyz_yyyy, g_z_0_y_0_yyz_yyyyy, g_z_0_y_0_yyz_yyyyz, g_z_0_y_0_yyz_yyyz, g_z_0_y_0_yyz_yyyzz, g_z_0_y_0_yyz_yyzz, g_z_0_y_0_yyz_yyzzz, g_z_0_y_0_yyz_yzzz, g_z_0_y_0_yyz_yzzzz, g_z_0_y_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyz_xxxx[k] = -g_z_0_y_0_yyz_xxxx[k] * ab_y + g_z_0_y_0_yyz_xxxxy[k];

                g_z_0_y_0_yyyz_xxxy[k] = -g_z_0_y_0_yyz_xxxy[k] * ab_y + g_z_0_y_0_yyz_xxxyy[k];

                g_z_0_y_0_yyyz_xxxz[k] = -g_z_0_y_0_yyz_xxxz[k] * ab_y + g_z_0_y_0_yyz_xxxyz[k];

                g_z_0_y_0_yyyz_xxyy[k] = -g_z_0_y_0_yyz_xxyy[k] * ab_y + g_z_0_y_0_yyz_xxyyy[k];

                g_z_0_y_0_yyyz_xxyz[k] = -g_z_0_y_0_yyz_xxyz[k] * ab_y + g_z_0_y_0_yyz_xxyyz[k];

                g_z_0_y_0_yyyz_xxzz[k] = -g_z_0_y_0_yyz_xxzz[k] * ab_y + g_z_0_y_0_yyz_xxyzz[k];

                g_z_0_y_0_yyyz_xyyy[k] = -g_z_0_y_0_yyz_xyyy[k] * ab_y + g_z_0_y_0_yyz_xyyyy[k];

                g_z_0_y_0_yyyz_xyyz[k] = -g_z_0_y_0_yyz_xyyz[k] * ab_y + g_z_0_y_0_yyz_xyyyz[k];

                g_z_0_y_0_yyyz_xyzz[k] = -g_z_0_y_0_yyz_xyzz[k] * ab_y + g_z_0_y_0_yyz_xyyzz[k];

                g_z_0_y_0_yyyz_xzzz[k] = -g_z_0_y_0_yyz_xzzz[k] * ab_y + g_z_0_y_0_yyz_xyzzz[k];

                g_z_0_y_0_yyyz_yyyy[k] = -g_z_0_y_0_yyz_yyyy[k] * ab_y + g_z_0_y_0_yyz_yyyyy[k];

                g_z_0_y_0_yyyz_yyyz[k] = -g_z_0_y_0_yyz_yyyz[k] * ab_y + g_z_0_y_0_yyz_yyyyz[k];

                g_z_0_y_0_yyyz_yyzz[k] = -g_z_0_y_0_yyz_yyzz[k] * ab_y + g_z_0_y_0_yyz_yyyzz[k];

                g_z_0_y_0_yyyz_yzzz[k] = -g_z_0_y_0_yyz_yzzz[k] * ab_y + g_z_0_y_0_yyz_yyzzz[k];

                g_z_0_y_0_yyyz_zzzz[k] = -g_z_0_y_0_yyz_zzzz[k] * ab_y + g_z_0_y_0_yyz_yzzzz[k];
            }

            /// Set up 1755-1770 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1755 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1756 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1757 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1758 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1759 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1760 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1761 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1762 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1763 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1764 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1765 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1766 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1767 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1768 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1769 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyzz_xxxx, g_z_0_y_0_yyzz_xxxy, g_z_0_y_0_yyzz_xxxz, g_z_0_y_0_yyzz_xxyy, g_z_0_y_0_yyzz_xxyz, g_z_0_y_0_yyzz_xxzz, g_z_0_y_0_yyzz_xyyy, g_z_0_y_0_yyzz_xyyz, g_z_0_y_0_yyzz_xyzz, g_z_0_y_0_yyzz_xzzz, g_z_0_y_0_yyzz_yyyy, g_z_0_y_0_yyzz_yyyz, g_z_0_y_0_yyzz_yyzz, g_z_0_y_0_yyzz_yzzz, g_z_0_y_0_yyzz_zzzz, g_z_0_y_0_yzz_xxxx, g_z_0_y_0_yzz_xxxxy, g_z_0_y_0_yzz_xxxy, g_z_0_y_0_yzz_xxxyy, g_z_0_y_0_yzz_xxxyz, g_z_0_y_0_yzz_xxxz, g_z_0_y_0_yzz_xxyy, g_z_0_y_0_yzz_xxyyy, g_z_0_y_0_yzz_xxyyz, g_z_0_y_0_yzz_xxyz, g_z_0_y_0_yzz_xxyzz, g_z_0_y_0_yzz_xxzz, g_z_0_y_0_yzz_xyyy, g_z_0_y_0_yzz_xyyyy, g_z_0_y_0_yzz_xyyyz, g_z_0_y_0_yzz_xyyz, g_z_0_y_0_yzz_xyyzz, g_z_0_y_0_yzz_xyzz, g_z_0_y_0_yzz_xyzzz, g_z_0_y_0_yzz_xzzz, g_z_0_y_0_yzz_yyyy, g_z_0_y_0_yzz_yyyyy, g_z_0_y_0_yzz_yyyyz, g_z_0_y_0_yzz_yyyz, g_z_0_y_0_yzz_yyyzz, g_z_0_y_0_yzz_yyzz, g_z_0_y_0_yzz_yyzzz, g_z_0_y_0_yzz_yzzz, g_z_0_y_0_yzz_yzzzz, g_z_0_y_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyzz_xxxx[k] = -g_z_0_y_0_yzz_xxxx[k] * ab_y + g_z_0_y_0_yzz_xxxxy[k];

                g_z_0_y_0_yyzz_xxxy[k] = -g_z_0_y_0_yzz_xxxy[k] * ab_y + g_z_0_y_0_yzz_xxxyy[k];

                g_z_0_y_0_yyzz_xxxz[k] = -g_z_0_y_0_yzz_xxxz[k] * ab_y + g_z_0_y_0_yzz_xxxyz[k];

                g_z_0_y_0_yyzz_xxyy[k] = -g_z_0_y_0_yzz_xxyy[k] * ab_y + g_z_0_y_0_yzz_xxyyy[k];

                g_z_0_y_0_yyzz_xxyz[k] = -g_z_0_y_0_yzz_xxyz[k] * ab_y + g_z_0_y_0_yzz_xxyyz[k];

                g_z_0_y_0_yyzz_xxzz[k] = -g_z_0_y_0_yzz_xxzz[k] * ab_y + g_z_0_y_0_yzz_xxyzz[k];

                g_z_0_y_0_yyzz_xyyy[k] = -g_z_0_y_0_yzz_xyyy[k] * ab_y + g_z_0_y_0_yzz_xyyyy[k];

                g_z_0_y_0_yyzz_xyyz[k] = -g_z_0_y_0_yzz_xyyz[k] * ab_y + g_z_0_y_0_yzz_xyyyz[k];

                g_z_0_y_0_yyzz_xyzz[k] = -g_z_0_y_0_yzz_xyzz[k] * ab_y + g_z_0_y_0_yzz_xyyzz[k];

                g_z_0_y_0_yyzz_xzzz[k] = -g_z_0_y_0_yzz_xzzz[k] * ab_y + g_z_0_y_0_yzz_xyzzz[k];

                g_z_0_y_0_yyzz_yyyy[k] = -g_z_0_y_0_yzz_yyyy[k] * ab_y + g_z_0_y_0_yzz_yyyyy[k];

                g_z_0_y_0_yyzz_yyyz[k] = -g_z_0_y_0_yzz_yyyz[k] * ab_y + g_z_0_y_0_yzz_yyyyz[k];

                g_z_0_y_0_yyzz_yyzz[k] = -g_z_0_y_0_yzz_yyzz[k] * ab_y + g_z_0_y_0_yzz_yyyzz[k];

                g_z_0_y_0_yyzz_yzzz[k] = -g_z_0_y_0_yzz_yzzz[k] * ab_y + g_z_0_y_0_yzz_yyzzz[k];

                g_z_0_y_0_yyzz_zzzz[k] = -g_z_0_y_0_yzz_zzzz[k] * ab_y + g_z_0_y_0_yzz_yzzzz[k];
            }

            /// Set up 1770-1785 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1770 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1771 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1772 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1773 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1774 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1775 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1776 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1777 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1778 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1779 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1780 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1781 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1782 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1783 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1784 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yzzz_xxxx, g_z_0_y_0_yzzz_xxxy, g_z_0_y_0_yzzz_xxxz, g_z_0_y_0_yzzz_xxyy, g_z_0_y_0_yzzz_xxyz, g_z_0_y_0_yzzz_xxzz, g_z_0_y_0_yzzz_xyyy, g_z_0_y_0_yzzz_xyyz, g_z_0_y_0_yzzz_xyzz, g_z_0_y_0_yzzz_xzzz, g_z_0_y_0_yzzz_yyyy, g_z_0_y_0_yzzz_yyyz, g_z_0_y_0_yzzz_yyzz, g_z_0_y_0_yzzz_yzzz, g_z_0_y_0_yzzz_zzzz, g_z_0_y_0_zzz_xxxx, g_z_0_y_0_zzz_xxxxy, g_z_0_y_0_zzz_xxxy, g_z_0_y_0_zzz_xxxyy, g_z_0_y_0_zzz_xxxyz, g_z_0_y_0_zzz_xxxz, g_z_0_y_0_zzz_xxyy, g_z_0_y_0_zzz_xxyyy, g_z_0_y_0_zzz_xxyyz, g_z_0_y_0_zzz_xxyz, g_z_0_y_0_zzz_xxyzz, g_z_0_y_0_zzz_xxzz, g_z_0_y_0_zzz_xyyy, g_z_0_y_0_zzz_xyyyy, g_z_0_y_0_zzz_xyyyz, g_z_0_y_0_zzz_xyyz, g_z_0_y_0_zzz_xyyzz, g_z_0_y_0_zzz_xyzz, g_z_0_y_0_zzz_xyzzz, g_z_0_y_0_zzz_xzzz, g_z_0_y_0_zzz_yyyy, g_z_0_y_0_zzz_yyyyy, g_z_0_y_0_zzz_yyyyz, g_z_0_y_0_zzz_yyyz, g_z_0_y_0_zzz_yyyzz, g_z_0_y_0_zzz_yyzz, g_z_0_y_0_zzz_yyzzz, g_z_0_y_0_zzz_yzzz, g_z_0_y_0_zzz_yzzzz, g_z_0_y_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yzzz_xxxx[k] = -g_z_0_y_0_zzz_xxxx[k] * ab_y + g_z_0_y_0_zzz_xxxxy[k];

                g_z_0_y_0_yzzz_xxxy[k] = -g_z_0_y_0_zzz_xxxy[k] * ab_y + g_z_0_y_0_zzz_xxxyy[k];

                g_z_0_y_0_yzzz_xxxz[k] = -g_z_0_y_0_zzz_xxxz[k] * ab_y + g_z_0_y_0_zzz_xxxyz[k];

                g_z_0_y_0_yzzz_xxyy[k] = -g_z_0_y_0_zzz_xxyy[k] * ab_y + g_z_0_y_0_zzz_xxyyy[k];

                g_z_0_y_0_yzzz_xxyz[k] = -g_z_0_y_0_zzz_xxyz[k] * ab_y + g_z_0_y_0_zzz_xxyyz[k];

                g_z_0_y_0_yzzz_xxzz[k] = -g_z_0_y_0_zzz_xxzz[k] * ab_y + g_z_0_y_0_zzz_xxyzz[k];

                g_z_0_y_0_yzzz_xyyy[k] = -g_z_0_y_0_zzz_xyyy[k] * ab_y + g_z_0_y_0_zzz_xyyyy[k];

                g_z_0_y_0_yzzz_xyyz[k] = -g_z_0_y_0_zzz_xyyz[k] * ab_y + g_z_0_y_0_zzz_xyyyz[k];

                g_z_0_y_0_yzzz_xyzz[k] = -g_z_0_y_0_zzz_xyzz[k] * ab_y + g_z_0_y_0_zzz_xyyzz[k];

                g_z_0_y_0_yzzz_xzzz[k] = -g_z_0_y_0_zzz_xzzz[k] * ab_y + g_z_0_y_0_zzz_xyzzz[k];

                g_z_0_y_0_yzzz_yyyy[k] = -g_z_0_y_0_zzz_yyyy[k] * ab_y + g_z_0_y_0_zzz_yyyyy[k];

                g_z_0_y_0_yzzz_yyyz[k] = -g_z_0_y_0_zzz_yyyz[k] * ab_y + g_z_0_y_0_zzz_yyyyz[k];

                g_z_0_y_0_yzzz_yyzz[k] = -g_z_0_y_0_zzz_yyzz[k] * ab_y + g_z_0_y_0_zzz_yyyzz[k];

                g_z_0_y_0_yzzz_yzzz[k] = -g_z_0_y_0_zzz_yzzz[k] * ab_y + g_z_0_y_0_zzz_yyzzz[k];

                g_z_0_y_0_yzzz_zzzz[k] = -g_z_0_y_0_zzz_zzzz[k] * ab_y + g_z_0_y_0_zzz_yzzzz[k];
            }

            /// Set up 1785-1800 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_zzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1785 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1786 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1787 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1788 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1789 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1790 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1791 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1792 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1793 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1794 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1795 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1796 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1797 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1798 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1799 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_zzz_xxxx, g_0_0_y_0_zzz_xxxy, g_0_0_y_0_zzz_xxxz, g_0_0_y_0_zzz_xxyy, g_0_0_y_0_zzz_xxyz, g_0_0_y_0_zzz_xxzz, g_0_0_y_0_zzz_xyyy, g_0_0_y_0_zzz_xyyz, g_0_0_y_0_zzz_xyzz, g_0_0_y_0_zzz_xzzz, g_0_0_y_0_zzz_yyyy, g_0_0_y_0_zzz_yyyz, g_0_0_y_0_zzz_yyzz, g_0_0_y_0_zzz_yzzz, g_0_0_y_0_zzz_zzzz, g_z_0_y_0_zzz_xxxx, g_z_0_y_0_zzz_xxxxz, g_z_0_y_0_zzz_xxxy, g_z_0_y_0_zzz_xxxyz, g_z_0_y_0_zzz_xxxz, g_z_0_y_0_zzz_xxxzz, g_z_0_y_0_zzz_xxyy, g_z_0_y_0_zzz_xxyyz, g_z_0_y_0_zzz_xxyz, g_z_0_y_0_zzz_xxyzz, g_z_0_y_0_zzz_xxzz, g_z_0_y_0_zzz_xxzzz, g_z_0_y_0_zzz_xyyy, g_z_0_y_0_zzz_xyyyz, g_z_0_y_0_zzz_xyyz, g_z_0_y_0_zzz_xyyzz, g_z_0_y_0_zzz_xyzz, g_z_0_y_0_zzz_xyzzz, g_z_0_y_0_zzz_xzzz, g_z_0_y_0_zzz_xzzzz, g_z_0_y_0_zzz_yyyy, g_z_0_y_0_zzz_yyyyz, g_z_0_y_0_zzz_yyyz, g_z_0_y_0_zzz_yyyzz, g_z_0_y_0_zzz_yyzz, g_z_0_y_0_zzz_yyzzz, g_z_0_y_0_zzz_yzzz, g_z_0_y_0_zzz_yzzzz, g_z_0_y_0_zzz_zzzz, g_z_0_y_0_zzz_zzzzz, g_z_0_y_0_zzzz_xxxx, g_z_0_y_0_zzzz_xxxy, g_z_0_y_0_zzzz_xxxz, g_z_0_y_0_zzzz_xxyy, g_z_0_y_0_zzzz_xxyz, g_z_0_y_0_zzzz_xxzz, g_z_0_y_0_zzzz_xyyy, g_z_0_y_0_zzzz_xyyz, g_z_0_y_0_zzzz_xyzz, g_z_0_y_0_zzzz_xzzz, g_z_0_y_0_zzzz_yyyy, g_z_0_y_0_zzzz_yyyz, g_z_0_y_0_zzzz_yyzz, g_z_0_y_0_zzzz_yzzz, g_z_0_y_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_zzzz_xxxx[k] = -g_0_0_y_0_zzz_xxxx[k] - g_z_0_y_0_zzz_xxxx[k] * ab_z + g_z_0_y_0_zzz_xxxxz[k];

                g_z_0_y_0_zzzz_xxxy[k] = -g_0_0_y_0_zzz_xxxy[k] - g_z_0_y_0_zzz_xxxy[k] * ab_z + g_z_0_y_0_zzz_xxxyz[k];

                g_z_0_y_0_zzzz_xxxz[k] = -g_0_0_y_0_zzz_xxxz[k] - g_z_0_y_0_zzz_xxxz[k] * ab_z + g_z_0_y_0_zzz_xxxzz[k];

                g_z_0_y_0_zzzz_xxyy[k] = -g_0_0_y_0_zzz_xxyy[k] - g_z_0_y_0_zzz_xxyy[k] * ab_z + g_z_0_y_0_zzz_xxyyz[k];

                g_z_0_y_0_zzzz_xxyz[k] = -g_0_0_y_0_zzz_xxyz[k] - g_z_0_y_0_zzz_xxyz[k] * ab_z + g_z_0_y_0_zzz_xxyzz[k];

                g_z_0_y_0_zzzz_xxzz[k] = -g_0_0_y_0_zzz_xxzz[k] - g_z_0_y_0_zzz_xxzz[k] * ab_z + g_z_0_y_0_zzz_xxzzz[k];

                g_z_0_y_0_zzzz_xyyy[k] = -g_0_0_y_0_zzz_xyyy[k] - g_z_0_y_0_zzz_xyyy[k] * ab_z + g_z_0_y_0_zzz_xyyyz[k];

                g_z_0_y_0_zzzz_xyyz[k] = -g_0_0_y_0_zzz_xyyz[k] - g_z_0_y_0_zzz_xyyz[k] * ab_z + g_z_0_y_0_zzz_xyyzz[k];

                g_z_0_y_0_zzzz_xyzz[k] = -g_0_0_y_0_zzz_xyzz[k] - g_z_0_y_0_zzz_xyzz[k] * ab_z + g_z_0_y_0_zzz_xyzzz[k];

                g_z_0_y_0_zzzz_xzzz[k] = -g_0_0_y_0_zzz_xzzz[k] - g_z_0_y_0_zzz_xzzz[k] * ab_z + g_z_0_y_0_zzz_xzzzz[k];

                g_z_0_y_0_zzzz_yyyy[k] = -g_0_0_y_0_zzz_yyyy[k] - g_z_0_y_0_zzz_yyyy[k] * ab_z + g_z_0_y_0_zzz_yyyyz[k];

                g_z_0_y_0_zzzz_yyyz[k] = -g_0_0_y_0_zzz_yyyz[k] - g_z_0_y_0_zzz_yyyz[k] * ab_z + g_z_0_y_0_zzz_yyyzz[k];

                g_z_0_y_0_zzzz_yyzz[k] = -g_0_0_y_0_zzz_yyzz[k] - g_z_0_y_0_zzz_yyzz[k] * ab_z + g_z_0_y_0_zzz_yyzzz[k];

                g_z_0_y_0_zzzz_yzzz[k] = -g_0_0_y_0_zzz_yzzz[k] - g_z_0_y_0_zzz_yzzz[k] * ab_z + g_z_0_y_0_zzz_yzzzz[k];

                g_z_0_y_0_zzzz_zzzz[k] = -g_0_0_y_0_zzz_zzzz[k] - g_z_0_y_0_zzz_zzzz[k] * ab_z + g_z_0_y_0_zzz_zzzzz[k];
            }

            /// Set up 1800-1815 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxx_xxxx = cbuffer.data(gg_geom_1010_off + 1800 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_xxxy = cbuffer.data(gg_geom_1010_off + 1801 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_xxxz = cbuffer.data(gg_geom_1010_off + 1802 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_xxyy = cbuffer.data(gg_geom_1010_off + 1803 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_xxyz = cbuffer.data(gg_geom_1010_off + 1804 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_xxzz = cbuffer.data(gg_geom_1010_off + 1805 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_xyyy = cbuffer.data(gg_geom_1010_off + 1806 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_xyyz = cbuffer.data(gg_geom_1010_off + 1807 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_xyzz = cbuffer.data(gg_geom_1010_off + 1808 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_xzzz = cbuffer.data(gg_geom_1010_off + 1809 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_yyyy = cbuffer.data(gg_geom_1010_off + 1810 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_yyyz = cbuffer.data(gg_geom_1010_off + 1811 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_yyzz = cbuffer.data(gg_geom_1010_off + 1812 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_yzzz = cbuffer.data(gg_geom_1010_off + 1813 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_zzzz = cbuffer.data(gg_geom_1010_off + 1814 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxx_xxxx, g_z_0_z_0_xxx_xxxxx, g_z_0_z_0_xxx_xxxxy, g_z_0_z_0_xxx_xxxxz, g_z_0_z_0_xxx_xxxy, g_z_0_z_0_xxx_xxxyy, g_z_0_z_0_xxx_xxxyz, g_z_0_z_0_xxx_xxxz, g_z_0_z_0_xxx_xxxzz, g_z_0_z_0_xxx_xxyy, g_z_0_z_0_xxx_xxyyy, g_z_0_z_0_xxx_xxyyz, g_z_0_z_0_xxx_xxyz, g_z_0_z_0_xxx_xxyzz, g_z_0_z_0_xxx_xxzz, g_z_0_z_0_xxx_xxzzz, g_z_0_z_0_xxx_xyyy, g_z_0_z_0_xxx_xyyyy, g_z_0_z_0_xxx_xyyyz, g_z_0_z_0_xxx_xyyz, g_z_0_z_0_xxx_xyyzz, g_z_0_z_0_xxx_xyzz, g_z_0_z_0_xxx_xyzzz, g_z_0_z_0_xxx_xzzz, g_z_0_z_0_xxx_xzzzz, g_z_0_z_0_xxx_yyyy, g_z_0_z_0_xxx_yyyz, g_z_0_z_0_xxx_yyzz, g_z_0_z_0_xxx_yzzz, g_z_0_z_0_xxx_zzzz, g_z_0_z_0_xxxx_xxxx, g_z_0_z_0_xxxx_xxxy, g_z_0_z_0_xxxx_xxxz, g_z_0_z_0_xxxx_xxyy, g_z_0_z_0_xxxx_xxyz, g_z_0_z_0_xxxx_xxzz, g_z_0_z_0_xxxx_xyyy, g_z_0_z_0_xxxx_xyyz, g_z_0_z_0_xxxx_xyzz, g_z_0_z_0_xxxx_xzzz, g_z_0_z_0_xxxx_yyyy, g_z_0_z_0_xxxx_yyyz, g_z_0_z_0_xxxx_yyzz, g_z_0_z_0_xxxx_yzzz, g_z_0_z_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxx_xxxx[k] = -g_z_0_z_0_xxx_xxxx[k] * ab_x + g_z_0_z_0_xxx_xxxxx[k];

                g_z_0_z_0_xxxx_xxxy[k] = -g_z_0_z_0_xxx_xxxy[k] * ab_x + g_z_0_z_0_xxx_xxxxy[k];

                g_z_0_z_0_xxxx_xxxz[k] = -g_z_0_z_0_xxx_xxxz[k] * ab_x + g_z_0_z_0_xxx_xxxxz[k];

                g_z_0_z_0_xxxx_xxyy[k] = -g_z_0_z_0_xxx_xxyy[k] * ab_x + g_z_0_z_0_xxx_xxxyy[k];

                g_z_0_z_0_xxxx_xxyz[k] = -g_z_0_z_0_xxx_xxyz[k] * ab_x + g_z_0_z_0_xxx_xxxyz[k];

                g_z_0_z_0_xxxx_xxzz[k] = -g_z_0_z_0_xxx_xxzz[k] * ab_x + g_z_0_z_0_xxx_xxxzz[k];

                g_z_0_z_0_xxxx_xyyy[k] = -g_z_0_z_0_xxx_xyyy[k] * ab_x + g_z_0_z_0_xxx_xxyyy[k];

                g_z_0_z_0_xxxx_xyyz[k] = -g_z_0_z_0_xxx_xyyz[k] * ab_x + g_z_0_z_0_xxx_xxyyz[k];

                g_z_0_z_0_xxxx_xyzz[k] = -g_z_0_z_0_xxx_xyzz[k] * ab_x + g_z_0_z_0_xxx_xxyzz[k];

                g_z_0_z_0_xxxx_xzzz[k] = -g_z_0_z_0_xxx_xzzz[k] * ab_x + g_z_0_z_0_xxx_xxzzz[k];

                g_z_0_z_0_xxxx_yyyy[k] = -g_z_0_z_0_xxx_yyyy[k] * ab_x + g_z_0_z_0_xxx_xyyyy[k];

                g_z_0_z_0_xxxx_yyyz[k] = -g_z_0_z_0_xxx_yyyz[k] * ab_x + g_z_0_z_0_xxx_xyyyz[k];

                g_z_0_z_0_xxxx_yyzz[k] = -g_z_0_z_0_xxx_yyzz[k] * ab_x + g_z_0_z_0_xxx_xyyzz[k];

                g_z_0_z_0_xxxx_yzzz[k] = -g_z_0_z_0_xxx_yzzz[k] * ab_x + g_z_0_z_0_xxx_xyzzz[k];

                g_z_0_z_0_xxxx_zzzz[k] = -g_z_0_z_0_xxx_zzzz[k] * ab_x + g_z_0_z_0_xxx_xzzzz[k];
            }

            /// Set up 1815-1830 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxy_xxxx = cbuffer.data(gg_geom_1010_off + 1815 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_xxxy = cbuffer.data(gg_geom_1010_off + 1816 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_xxxz = cbuffer.data(gg_geom_1010_off + 1817 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_xxyy = cbuffer.data(gg_geom_1010_off + 1818 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_xxyz = cbuffer.data(gg_geom_1010_off + 1819 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_xxzz = cbuffer.data(gg_geom_1010_off + 1820 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_xyyy = cbuffer.data(gg_geom_1010_off + 1821 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_xyyz = cbuffer.data(gg_geom_1010_off + 1822 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_xyzz = cbuffer.data(gg_geom_1010_off + 1823 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_xzzz = cbuffer.data(gg_geom_1010_off + 1824 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_yyyy = cbuffer.data(gg_geom_1010_off + 1825 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_yyyz = cbuffer.data(gg_geom_1010_off + 1826 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_yyzz = cbuffer.data(gg_geom_1010_off + 1827 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_yzzz = cbuffer.data(gg_geom_1010_off + 1828 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_zzzz = cbuffer.data(gg_geom_1010_off + 1829 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxy_xxxx, g_z_0_z_0_xxxy_xxxy, g_z_0_z_0_xxxy_xxxz, g_z_0_z_0_xxxy_xxyy, g_z_0_z_0_xxxy_xxyz, g_z_0_z_0_xxxy_xxzz, g_z_0_z_0_xxxy_xyyy, g_z_0_z_0_xxxy_xyyz, g_z_0_z_0_xxxy_xyzz, g_z_0_z_0_xxxy_xzzz, g_z_0_z_0_xxxy_yyyy, g_z_0_z_0_xxxy_yyyz, g_z_0_z_0_xxxy_yyzz, g_z_0_z_0_xxxy_yzzz, g_z_0_z_0_xxxy_zzzz, g_z_0_z_0_xxy_xxxx, g_z_0_z_0_xxy_xxxxx, g_z_0_z_0_xxy_xxxxy, g_z_0_z_0_xxy_xxxxz, g_z_0_z_0_xxy_xxxy, g_z_0_z_0_xxy_xxxyy, g_z_0_z_0_xxy_xxxyz, g_z_0_z_0_xxy_xxxz, g_z_0_z_0_xxy_xxxzz, g_z_0_z_0_xxy_xxyy, g_z_0_z_0_xxy_xxyyy, g_z_0_z_0_xxy_xxyyz, g_z_0_z_0_xxy_xxyz, g_z_0_z_0_xxy_xxyzz, g_z_0_z_0_xxy_xxzz, g_z_0_z_0_xxy_xxzzz, g_z_0_z_0_xxy_xyyy, g_z_0_z_0_xxy_xyyyy, g_z_0_z_0_xxy_xyyyz, g_z_0_z_0_xxy_xyyz, g_z_0_z_0_xxy_xyyzz, g_z_0_z_0_xxy_xyzz, g_z_0_z_0_xxy_xyzzz, g_z_0_z_0_xxy_xzzz, g_z_0_z_0_xxy_xzzzz, g_z_0_z_0_xxy_yyyy, g_z_0_z_0_xxy_yyyz, g_z_0_z_0_xxy_yyzz, g_z_0_z_0_xxy_yzzz, g_z_0_z_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxy_xxxx[k] = -g_z_0_z_0_xxy_xxxx[k] * ab_x + g_z_0_z_0_xxy_xxxxx[k];

                g_z_0_z_0_xxxy_xxxy[k] = -g_z_0_z_0_xxy_xxxy[k] * ab_x + g_z_0_z_0_xxy_xxxxy[k];

                g_z_0_z_0_xxxy_xxxz[k] = -g_z_0_z_0_xxy_xxxz[k] * ab_x + g_z_0_z_0_xxy_xxxxz[k];

                g_z_0_z_0_xxxy_xxyy[k] = -g_z_0_z_0_xxy_xxyy[k] * ab_x + g_z_0_z_0_xxy_xxxyy[k];

                g_z_0_z_0_xxxy_xxyz[k] = -g_z_0_z_0_xxy_xxyz[k] * ab_x + g_z_0_z_0_xxy_xxxyz[k];

                g_z_0_z_0_xxxy_xxzz[k] = -g_z_0_z_0_xxy_xxzz[k] * ab_x + g_z_0_z_0_xxy_xxxzz[k];

                g_z_0_z_0_xxxy_xyyy[k] = -g_z_0_z_0_xxy_xyyy[k] * ab_x + g_z_0_z_0_xxy_xxyyy[k];

                g_z_0_z_0_xxxy_xyyz[k] = -g_z_0_z_0_xxy_xyyz[k] * ab_x + g_z_0_z_0_xxy_xxyyz[k];

                g_z_0_z_0_xxxy_xyzz[k] = -g_z_0_z_0_xxy_xyzz[k] * ab_x + g_z_0_z_0_xxy_xxyzz[k];

                g_z_0_z_0_xxxy_xzzz[k] = -g_z_0_z_0_xxy_xzzz[k] * ab_x + g_z_0_z_0_xxy_xxzzz[k];

                g_z_0_z_0_xxxy_yyyy[k] = -g_z_0_z_0_xxy_yyyy[k] * ab_x + g_z_0_z_0_xxy_xyyyy[k];

                g_z_0_z_0_xxxy_yyyz[k] = -g_z_0_z_0_xxy_yyyz[k] * ab_x + g_z_0_z_0_xxy_xyyyz[k];

                g_z_0_z_0_xxxy_yyzz[k] = -g_z_0_z_0_xxy_yyzz[k] * ab_x + g_z_0_z_0_xxy_xyyzz[k];

                g_z_0_z_0_xxxy_yzzz[k] = -g_z_0_z_0_xxy_yzzz[k] * ab_x + g_z_0_z_0_xxy_xyzzz[k];

                g_z_0_z_0_xxxy_zzzz[k] = -g_z_0_z_0_xxy_zzzz[k] * ab_x + g_z_0_z_0_xxy_xzzzz[k];
            }

            /// Set up 1830-1845 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxz_xxxx = cbuffer.data(gg_geom_1010_off + 1830 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_xxxy = cbuffer.data(gg_geom_1010_off + 1831 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_xxxz = cbuffer.data(gg_geom_1010_off + 1832 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_xxyy = cbuffer.data(gg_geom_1010_off + 1833 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_xxyz = cbuffer.data(gg_geom_1010_off + 1834 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_xxzz = cbuffer.data(gg_geom_1010_off + 1835 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_xyyy = cbuffer.data(gg_geom_1010_off + 1836 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_xyyz = cbuffer.data(gg_geom_1010_off + 1837 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_xyzz = cbuffer.data(gg_geom_1010_off + 1838 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_xzzz = cbuffer.data(gg_geom_1010_off + 1839 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_yyyy = cbuffer.data(gg_geom_1010_off + 1840 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_yyyz = cbuffer.data(gg_geom_1010_off + 1841 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_yyzz = cbuffer.data(gg_geom_1010_off + 1842 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_yzzz = cbuffer.data(gg_geom_1010_off + 1843 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_zzzz = cbuffer.data(gg_geom_1010_off + 1844 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxz_xxxx, g_z_0_z_0_xxxz_xxxy, g_z_0_z_0_xxxz_xxxz, g_z_0_z_0_xxxz_xxyy, g_z_0_z_0_xxxz_xxyz, g_z_0_z_0_xxxz_xxzz, g_z_0_z_0_xxxz_xyyy, g_z_0_z_0_xxxz_xyyz, g_z_0_z_0_xxxz_xyzz, g_z_0_z_0_xxxz_xzzz, g_z_0_z_0_xxxz_yyyy, g_z_0_z_0_xxxz_yyyz, g_z_0_z_0_xxxz_yyzz, g_z_0_z_0_xxxz_yzzz, g_z_0_z_0_xxxz_zzzz, g_z_0_z_0_xxz_xxxx, g_z_0_z_0_xxz_xxxxx, g_z_0_z_0_xxz_xxxxy, g_z_0_z_0_xxz_xxxxz, g_z_0_z_0_xxz_xxxy, g_z_0_z_0_xxz_xxxyy, g_z_0_z_0_xxz_xxxyz, g_z_0_z_0_xxz_xxxz, g_z_0_z_0_xxz_xxxzz, g_z_0_z_0_xxz_xxyy, g_z_0_z_0_xxz_xxyyy, g_z_0_z_0_xxz_xxyyz, g_z_0_z_0_xxz_xxyz, g_z_0_z_0_xxz_xxyzz, g_z_0_z_0_xxz_xxzz, g_z_0_z_0_xxz_xxzzz, g_z_0_z_0_xxz_xyyy, g_z_0_z_0_xxz_xyyyy, g_z_0_z_0_xxz_xyyyz, g_z_0_z_0_xxz_xyyz, g_z_0_z_0_xxz_xyyzz, g_z_0_z_0_xxz_xyzz, g_z_0_z_0_xxz_xyzzz, g_z_0_z_0_xxz_xzzz, g_z_0_z_0_xxz_xzzzz, g_z_0_z_0_xxz_yyyy, g_z_0_z_0_xxz_yyyz, g_z_0_z_0_xxz_yyzz, g_z_0_z_0_xxz_yzzz, g_z_0_z_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxz_xxxx[k] = -g_z_0_z_0_xxz_xxxx[k] * ab_x + g_z_0_z_0_xxz_xxxxx[k];

                g_z_0_z_0_xxxz_xxxy[k] = -g_z_0_z_0_xxz_xxxy[k] * ab_x + g_z_0_z_0_xxz_xxxxy[k];

                g_z_0_z_0_xxxz_xxxz[k] = -g_z_0_z_0_xxz_xxxz[k] * ab_x + g_z_0_z_0_xxz_xxxxz[k];

                g_z_0_z_0_xxxz_xxyy[k] = -g_z_0_z_0_xxz_xxyy[k] * ab_x + g_z_0_z_0_xxz_xxxyy[k];

                g_z_0_z_0_xxxz_xxyz[k] = -g_z_0_z_0_xxz_xxyz[k] * ab_x + g_z_0_z_0_xxz_xxxyz[k];

                g_z_0_z_0_xxxz_xxzz[k] = -g_z_0_z_0_xxz_xxzz[k] * ab_x + g_z_0_z_0_xxz_xxxzz[k];

                g_z_0_z_0_xxxz_xyyy[k] = -g_z_0_z_0_xxz_xyyy[k] * ab_x + g_z_0_z_0_xxz_xxyyy[k];

                g_z_0_z_0_xxxz_xyyz[k] = -g_z_0_z_0_xxz_xyyz[k] * ab_x + g_z_0_z_0_xxz_xxyyz[k];

                g_z_0_z_0_xxxz_xyzz[k] = -g_z_0_z_0_xxz_xyzz[k] * ab_x + g_z_0_z_0_xxz_xxyzz[k];

                g_z_0_z_0_xxxz_xzzz[k] = -g_z_0_z_0_xxz_xzzz[k] * ab_x + g_z_0_z_0_xxz_xxzzz[k];

                g_z_0_z_0_xxxz_yyyy[k] = -g_z_0_z_0_xxz_yyyy[k] * ab_x + g_z_0_z_0_xxz_xyyyy[k];

                g_z_0_z_0_xxxz_yyyz[k] = -g_z_0_z_0_xxz_yyyz[k] * ab_x + g_z_0_z_0_xxz_xyyyz[k];

                g_z_0_z_0_xxxz_yyzz[k] = -g_z_0_z_0_xxz_yyzz[k] * ab_x + g_z_0_z_0_xxz_xyyzz[k];

                g_z_0_z_0_xxxz_yzzz[k] = -g_z_0_z_0_xxz_yzzz[k] * ab_x + g_z_0_z_0_xxz_xyzzz[k];

                g_z_0_z_0_xxxz_zzzz[k] = -g_z_0_z_0_xxz_zzzz[k] * ab_x + g_z_0_z_0_xxz_xzzzz[k];
            }

            /// Set up 1845-1860 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyy_xxxx = cbuffer.data(gg_geom_1010_off + 1845 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_xxxy = cbuffer.data(gg_geom_1010_off + 1846 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_xxxz = cbuffer.data(gg_geom_1010_off + 1847 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_xxyy = cbuffer.data(gg_geom_1010_off + 1848 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_xxyz = cbuffer.data(gg_geom_1010_off + 1849 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_xxzz = cbuffer.data(gg_geom_1010_off + 1850 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_xyyy = cbuffer.data(gg_geom_1010_off + 1851 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_xyyz = cbuffer.data(gg_geom_1010_off + 1852 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_xyzz = cbuffer.data(gg_geom_1010_off + 1853 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_xzzz = cbuffer.data(gg_geom_1010_off + 1854 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_yyyy = cbuffer.data(gg_geom_1010_off + 1855 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_yyyz = cbuffer.data(gg_geom_1010_off + 1856 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_yyzz = cbuffer.data(gg_geom_1010_off + 1857 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_yzzz = cbuffer.data(gg_geom_1010_off + 1858 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_zzzz = cbuffer.data(gg_geom_1010_off + 1859 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyy_xxxx, g_z_0_z_0_xxyy_xxxy, g_z_0_z_0_xxyy_xxxz, g_z_0_z_0_xxyy_xxyy, g_z_0_z_0_xxyy_xxyz, g_z_0_z_0_xxyy_xxzz, g_z_0_z_0_xxyy_xyyy, g_z_0_z_0_xxyy_xyyz, g_z_0_z_0_xxyy_xyzz, g_z_0_z_0_xxyy_xzzz, g_z_0_z_0_xxyy_yyyy, g_z_0_z_0_xxyy_yyyz, g_z_0_z_0_xxyy_yyzz, g_z_0_z_0_xxyy_yzzz, g_z_0_z_0_xxyy_zzzz, g_z_0_z_0_xyy_xxxx, g_z_0_z_0_xyy_xxxxx, g_z_0_z_0_xyy_xxxxy, g_z_0_z_0_xyy_xxxxz, g_z_0_z_0_xyy_xxxy, g_z_0_z_0_xyy_xxxyy, g_z_0_z_0_xyy_xxxyz, g_z_0_z_0_xyy_xxxz, g_z_0_z_0_xyy_xxxzz, g_z_0_z_0_xyy_xxyy, g_z_0_z_0_xyy_xxyyy, g_z_0_z_0_xyy_xxyyz, g_z_0_z_0_xyy_xxyz, g_z_0_z_0_xyy_xxyzz, g_z_0_z_0_xyy_xxzz, g_z_0_z_0_xyy_xxzzz, g_z_0_z_0_xyy_xyyy, g_z_0_z_0_xyy_xyyyy, g_z_0_z_0_xyy_xyyyz, g_z_0_z_0_xyy_xyyz, g_z_0_z_0_xyy_xyyzz, g_z_0_z_0_xyy_xyzz, g_z_0_z_0_xyy_xyzzz, g_z_0_z_0_xyy_xzzz, g_z_0_z_0_xyy_xzzzz, g_z_0_z_0_xyy_yyyy, g_z_0_z_0_xyy_yyyz, g_z_0_z_0_xyy_yyzz, g_z_0_z_0_xyy_yzzz, g_z_0_z_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyy_xxxx[k] = -g_z_0_z_0_xyy_xxxx[k] * ab_x + g_z_0_z_0_xyy_xxxxx[k];

                g_z_0_z_0_xxyy_xxxy[k] = -g_z_0_z_0_xyy_xxxy[k] * ab_x + g_z_0_z_0_xyy_xxxxy[k];

                g_z_0_z_0_xxyy_xxxz[k] = -g_z_0_z_0_xyy_xxxz[k] * ab_x + g_z_0_z_0_xyy_xxxxz[k];

                g_z_0_z_0_xxyy_xxyy[k] = -g_z_0_z_0_xyy_xxyy[k] * ab_x + g_z_0_z_0_xyy_xxxyy[k];

                g_z_0_z_0_xxyy_xxyz[k] = -g_z_0_z_0_xyy_xxyz[k] * ab_x + g_z_0_z_0_xyy_xxxyz[k];

                g_z_0_z_0_xxyy_xxzz[k] = -g_z_0_z_0_xyy_xxzz[k] * ab_x + g_z_0_z_0_xyy_xxxzz[k];

                g_z_0_z_0_xxyy_xyyy[k] = -g_z_0_z_0_xyy_xyyy[k] * ab_x + g_z_0_z_0_xyy_xxyyy[k];

                g_z_0_z_0_xxyy_xyyz[k] = -g_z_0_z_0_xyy_xyyz[k] * ab_x + g_z_0_z_0_xyy_xxyyz[k];

                g_z_0_z_0_xxyy_xyzz[k] = -g_z_0_z_0_xyy_xyzz[k] * ab_x + g_z_0_z_0_xyy_xxyzz[k];

                g_z_0_z_0_xxyy_xzzz[k] = -g_z_0_z_0_xyy_xzzz[k] * ab_x + g_z_0_z_0_xyy_xxzzz[k];

                g_z_0_z_0_xxyy_yyyy[k] = -g_z_0_z_0_xyy_yyyy[k] * ab_x + g_z_0_z_0_xyy_xyyyy[k];

                g_z_0_z_0_xxyy_yyyz[k] = -g_z_0_z_0_xyy_yyyz[k] * ab_x + g_z_0_z_0_xyy_xyyyz[k];

                g_z_0_z_0_xxyy_yyzz[k] = -g_z_0_z_0_xyy_yyzz[k] * ab_x + g_z_0_z_0_xyy_xyyzz[k];

                g_z_0_z_0_xxyy_yzzz[k] = -g_z_0_z_0_xyy_yzzz[k] * ab_x + g_z_0_z_0_xyy_xyzzz[k];

                g_z_0_z_0_xxyy_zzzz[k] = -g_z_0_z_0_xyy_zzzz[k] * ab_x + g_z_0_z_0_xyy_xzzzz[k];
            }

            /// Set up 1860-1875 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyz_xxxx = cbuffer.data(gg_geom_1010_off + 1860 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_xxxy = cbuffer.data(gg_geom_1010_off + 1861 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_xxxz = cbuffer.data(gg_geom_1010_off + 1862 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_xxyy = cbuffer.data(gg_geom_1010_off + 1863 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_xxyz = cbuffer.data(gg_geom_1010_off + 1864 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_xxzz = cbuffer.data(gg_geom_1010_off + 1865 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_xyyy = cbuffer.data(gg_geom_1010_off + 1866 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_xyyz = cbuffer.data(gg_geom_1010_off + 1867 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_xyzz = cbuffer.data(gg_geom_1010_off + 1868 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_xzzz = cbuffer.data(gg_geom_1010_off + 1869 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_yyyy = cbuffer.data(gg_geom_1010_off + 1870 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_yyyz = cbuffer.data(gg_geom_1010_off + 1871 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_yyzz = cbuffer.data(gg_geom_1010_off + 1872 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_yzzz = cbuffer.data(gg_geom_1010_off + 1873 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_zzzz = cbuffer.data(gg_geom_1010_off + 1874 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyz_xxxx, g_z_0_z_0_xxyz_xxxy, g_z_0_z_0_xxyz_xxxz, g_z_0_z_0_xxyz_xxyy, g_z_0_z_0_xxyz_xxyz, g_z_0_z_0_xxyz_xxzz, g_z_0_z_0_xxyz_xyyy, g_z_0_z_0_xxyz_xyyz, g_z_0_z_0_xxyz_xyzz, g_z_0_z_0_xxyz_xzzz, g_z_0_z_0_xxyz_yyyy, g_z_0_z_0_xxyz_yyyz, g_z_0_z_0_xxyz_yyzz, g_z_0_z_0_xxyz_yzzz, g_z_0_z_0_xxyz_zzzz, g_z_0_z_0_xyz_xxxx, g_z_0_z_0_xyz_xxxxx, g_z_0_z_0_xyz_xxxxy, g_z_0_z_0_xyz_xxxxz, g_z_0_z_0_xyz_xxxy, g_z_0_z_0_xyz_xxxyy, g_z_0_z_0_xyz_xxxyz, g_z_0_z_0_xyz_xxxz, g_z_0_z_0_xyz_xxxzz, g_z_0_z_0_xyz_xxyy, g_z_0_z_0_xyz_xxyyy, g_z_0_z_0_xyz_xxyyz, g_z_0_z_0_xyz_xxyz, g_z_0_z_0_xyz_xxyzz, g_z_0_z_0_xyz_xxzz, g_z_0_z_0_xyz_xxzzz, g_z_0_z_0_xyz_xyyy, g_z_0_z_0_xyz_xyyyy, g_z_0_z_0_xyz_xyyyz, g_z_0_z_0_xyz_xyyz, g_z_0_z_0_xyz_xyyzz, g_z_0_z_0_xyz_xyzz, g_z_0_z_0_xyz_xyzzz, g_z_0_z_0_xyz_xzzz, g_z_0_z_0_xyz_xzzzz, g_z_0_z_0_xyz_yyyy, g_z_0_z_0_xyz_yyyz, g_z_0_z_0_xyz_yyzz, g_z_0_z_0_xyz_yzzz, g_z_0_z_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyz_xxxx[k] = -g_z_0_z_0_xyz_xxxx[k] * ab_x + g_z_0_z_0_xyz_xxxxx[k];

                g_z_0_z_0_xxyz_xxxy[k] = -g_z_0_z_0_xyz_xxxy[k] * ab_x + g_z_0_z_0_xyz_xxxxy[k];

                g_z_0_z_0_xxyz_xxxz[k] = -g_z_0_z_0_xyz_xxxz[k] * ab_x + g_z_0_z_0_xyz_xxxxz[k];

                g_z_0_z_0_xxyz_xxyy[k] = -g_z_0_z_0_xyz_xxyy[k] * ab_x + g_z_0_z_0_xyz_xxxyy[k];

                g_z_0_z_0_xxyz_xxyz[k] = -g_z_0_z_0_xyz_xxyz[k] * ab_x + g_z_0_z_0_xyz_xxxyz[k];

                g_z_0_z_0_xxyz_xxzz[k] = -g_z_0_z_0_xyz_xxzz[k] * ab_x + g_z_0_z_0_xyz_xxxzz[k];

                g_z_0_z_0_xxyz_xyyy[k] = -g_z_0_z_0_xyz_xyyy[k] * ab_x + g_z_0_z_0_xyz_xxyyy[k];

                g_z_0_z_0_xxyz_xyyz[k] = -g_z_0_z_0_xyz_xyyz[k] * ab_x + g_z_0_z_0_xyz_xxyyz[k];

                g_z_0_z_0_xxyz_xyzz[k] = -g_z_0_z_0_xyz_xyzz[k] * ab_x + g_z_0_z_0_xyz_xxyzz[k];

                g_z_0_z_0_xxyz_xzzz[k] = -g_z_0_z_0_xyz_xzzz[k] * ab_x + g_z_0_z_0_xyz_xxzzz[k];

                g_z_0_z_0_xxyz_yyyy[k] = -g_z_0_z_0_xyz_yyyy[k] * ab_x + g_z_0_z_0_xyz_xyyyy[k];

                g_z_0_z_0_xxyz_yyyz[k] = -g_z_0_z_0_xyz_yyyz[k] * ab_x + g_z_0_z_0_xyz_xyyyz[k];

                g_z_0_z_0_xxyz_yyzz[k] = -g_z_0_z_0_xyz_yyzz[k] * ab_x + g_z_0_z_0_xyz_xyyzz[k];

                g_z_0_z_0_xxyz_yzzz[k] = -g_z_0_z_0_xyz_yzzz[k] * ab_x + g_z_0_z_0_xyz_xyzzz[k];

                g_z_0_z_0_xxyz_zzzz[k] = -g_z_0_z_0_xyz_zzzz[k] * ab_x + g_z_0_z_0_xyz_xzzzz[k];
            }

            /// Set up 1875-1890 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxzz_xxxx = cbuffer.data(gg_geom_1010_off + 1875 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_xxxy = cbuffer.data(gg_geom_1010_off + 1876 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_xxxz = cbuffer.data(gg_geom_1010_off + 1877 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_xxyy = cbuffer.data(gg_geom_1010_off + 1878 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_xxyz = cbuffer.data(gg_geom_1010_off + 1879 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_xxzz = cbuffer.data(gg_geom_1010_off + 1880 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_xyyy = cbuffer.data(gg_geom_1010_off + 1881 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_xyyz = cbuffer.data(gg_geom_1010_off + 1882 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_xyzz = cbuffer.data(gg_geom_1010_off + 1883 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_xzzz = cbuffer.data(gg_geom_1010_off + 1884 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_yyyy = cbuffer.data(gg_geom_1010_off + 1885 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_yyyz = cbuffer.data(gg_geom_1010_off + 1886 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_yyzz = cbuffer.data(gg_geom_1010_off + 1887 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_yzzz = cbuffer.data(gg_geom_1010_off + 1888 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_zzzz = cbuffer.data(gg_geom_1010_off + 1889 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxzz_xxxx, g_z_0_z_0_xxzz_xxxy, g_z_0_z_0_xxzz_xxxz, g_z_0_z_0_xxzz_xxyy, g_z_0_z_0_xxzz_xxyz, g_z_0_z_0_xxzz_xxzz, g_z_0_z_0_xxzz_xyyy, g_z_0_z_0_xxzz_xyyz, g_z_0_z_0_xxzz_xyzz, g_z_0_z_0_xxzz_xzzz, g_z_0_z_0_xxzz_yyyy, g_z_0_z_0_xxzz_yyyz, g_z_0_z_0_xxzz_yyzz, g_z_0_z_0_xxzz_yzzz, g_z_0_z_0_xxzz_zzzz, g_z_0_z_0_xzz_xxxx, g_z_0_z_0_xzz_xxxxx, g_z_0_z_0_xzz_xxxxy, g_z_0_z_0_xzz_xxxxz, g_z_0_z_0_xzz_xxxy, g_z_0_z_0_xzz_xxxyy, g_z_0_z_0_xzz_xxxyz, g_z_0_z_0_xzz_xxxz, g_z_0_z_0_xzz_xxxzz, g_z_0_z_0_xzz_xxyy, g_z_0_z_0_xzz_xxyyy, g_z_0_z_0_xzz_xxyyz, g_z_0_z_0_xzz_xxyz, g_z_0_z_0_xzz_xxyzz, g_z_0_z_0_xzz_xxzz, g_z_0_z_0_xzz_xxzzz, g_z_0_z_0_xzz_xyyy, g_z_0_z_0_xzz_xyyyy, g_z_0_z_0_xzz_xyyyz, g_z_0_z_0_xzz_xyyz, g_z_0_z_0_xzz_xyyzz, g_z_0_z_0_xzz_xyzz, g_z_0_z_0_xzz_xyzzz, g_z_0_z_0_xzz_xzzz, g_z_0_z_0_xzz_xzzzz, g_z_0_z_0_xzz_yyyy, g_z_0_z_0_xzz_yyyz, g_z_0_z_0_xzz_yyzz, g_z_0_z_0_xzz_yzzz, g_z_0_z_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxzz_xxxx[k] = -g_z_0_z_0_xzz_xxxx[k] * ab_x + g_z_0_z_0_xzz_xxxxx[k];

                g_z_0_z_0_xxzz_xxxy[k] = -g_z_0_z_0_xzz_xxxy[k] * ab_x + g_z_0_z_0_xzz_xxxxy[k];

                g_z_0_z_0_xxzz_xxxz[k] = -g_z_0_z_0_xzz_xxxz[k] * ab_x + g_z_0_z_0_xzz_xxxxz[k];

                g_z_0_z_0_xxzz_xxyy[k] = -g_z_0_z_0_xzz_xxyy[k] * ab_x + g_z_0_z_0_xzz_xxxyy[k];

                g_z_0_z_0_xxzz_xxyz[k] = -g_z_0_z_0_xzz_xxyz[k] * ab_x + g_z_0_z_0_xzz_xxxyz[k];

                g_z_0_z_0_xxzz_xxzz[k] = -g_z_0_z_0_xzz_xxzz[k] * ab_x + g_z_0_z_0_xzz_xxxzz[k];

                g_z_0_z_0_xxzz_xyyy[k] = -g_z_0_z_0_xzz_xyyy[k] * ab_x + g_z_0_z_0_xzz_xxyyy[k];

                g_z_0_z_0_xxzz_xyyz[k] = -g_z_0_z_0_xzz_xyyz[k] * ab_x + g_z_0_z_0_xzz_xxyyz[k];

                g_z_0_z_0_xxzz_xyzz[k] = -g_z_0_z_0_xzz_xyzz[k] * ab_x + g_z_0_z_0_xzz_xxyzz[k];

                g_z_0_z_0_xxzz_xzzz[k] = -g_z_0_z_0_xzz_xzzz[k] * ab_x + g_z_0_z_0_xzz_xxzzz[k];

                g_z_0_z_0_xxzz_yyyy[k] = -g_z_0_z_0_xzz_yyyy[k] * ab_x + g_z_0_z_0_xzz_xyyyy[k];

                g_z_0_z_0_xxzz_yyyz[k] = -g_z_0_z_0_xzz_yyyz[k] * ab_x + g_z_0_z_0_xzz_xyyyz[k];

                g_z_0_z_0_xxzz_yyzz[k] = -g_z_0_z_0_xzz_yyzz[k] * ab_x + g_z_0_z_0_xzz_xyyzz[k];

                g_z_0_z_0_xxzz_yzzz[k] = -g_z_0_z_0_xzz_yzzz[k] * ab_x + g_z_0_z_0_xzz_xyzzz[k];

                g_z_0_z_0_xxzz_zzzz[k] = -g_z_0_z_0_xzz_zzzz[k] * ab_x + g_z_0_z_0_xzz_xzzzz[k];
            }

            /// Set up 1890-1905 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyy_xxxx = cbuffer.data(gg_geom_1010_off + 1890 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_xxxy = cbuffer.data(gg_geom_1010_off + 1891 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_xxxz = cbuffer.data(gg_geom_1010_off + 1892 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_xxyy = cbuffer.data(gg_geom_1010_off + 1893 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_xxyz = cbuffer.data(gg_geom_1010_off + 1894 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_xxzz = cbuffer.data(gg_geom_1010_off + 1895 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_xyyy = cbuffer.data(gg_geom_1010_off + 1896 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_xyyz = cbuffer.data(gg_geom_1010_off + 1897 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_xyzz = cbuffer.data(gg_geom_1010_off + 1898 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_xzzz = cbuffer.data(gg_geom_1010_off + 1899 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1900 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1901 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1902 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1903 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1904 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyy_xxxx, g_z_0_z_0_xyyy_xxxy, g_z_0_z_0_xyyy_xxxz, g_z_0_z_0_xyyy_xxyy, g_z_0_z_0_xyyy_xxyz, g_z_0_z_0_xyyy_xxzz, g_z_0_z_0_xyyy_xyyy, g_z_0_z_0_xyyy_xyyz, g_z_0_z_0_xyyy_xyzz, g_z_0_z_0_xyyy_xzzz, g_z_0_z_0_xyyy_yyyy, g_z_0_z_0_xyyy_yyyz, g_z_0_z_0_xyyy_yyzz, g_z_0_z_0_xyyy_yzzz, g_z_0_z_0_xyyy_zzzz, g_z_0_z_0_yyy_xxxx, g_z_0_z_0_yyy_xxxxx, g_z_0_z_0_yyy_xxxxy, g_z_0_z_0_yyy_xxxxz, g_z_0_z_0_yyy_xxxy, g_z_0_z_0_yyy_xxxyy, g_z_0_z_0_yyy_xxxyz, g_z_0_z_0_yyy_xxxz, g_z_0_z_0_yyy_xxxzz, g_z_0_z_0_yyy_xxyy, g_z_0_z_0_yyy_xxyyy, g_z_0_z_0_yyy_xxyyz, g_z_0_z_0_yyy_xxyz, g_z_0_z_0_yyy_xxyzz, g_z_0_z_0_yyy_xxzz, g_z_0_z_0_yyy_xxzzz, g_z_0_z_0_yyy_xyyy, g_z_0_z_0_yyy_xyyyy, g_z_0_z_0_yyy_xyyyz, g_z_0_z_0_yyy_xyyz, g_z_0_z_0_yyy_xyyzz, g_z_0_z_0_yyy_xyzz, g_z_0_z_0_yyy_xyzzz, g_z_0_z_0_yyy_xzzz, g_z_0_z_0_yyy_xzzzz, g_z_0_z_0_yyy_yyyy, g_z_0_z_0_yyy_yyyz, g_z_0_z_0_yyy_yyzz, g_z_0_z_0_yyy_yzzz, g_z_0_z_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyy_xxxx[k] = -g_z_0_z_0_yyy_xxxx[k] * ab_x + g_z_0_z_0_yyy_xxxxx[k];

                g_z_0_z_0_xyyy_xxxy[k] = -g_z_0_z_0_yyy_xxxy[k] * ab_x + g_z_0_z_0_yyy_xxxxy[k];

                g_z_0_z_0_xyyy_xxxz[k] = -g_z_0_z_0_yyy_xxxz[k] * ab_x + g_z_0_z_0_yyy_xxxxz[k];

                g_z_0_z_0_xyyy_xxyy[k] = -g_z_0_z_0_yyy_xxyy[k] * ab_x + g_z_0_z_0_yyy_xxxyy[k];

                g_z_0_z_0_xyyy_xxyz[k] = -g_z_0_z_0_yyy_xxyz[k] * ab_x + g_z_0_z_0_yyy_xxxyz[k];

                g_z_0_z_0_xyyy_xxzz[k] = -g_z_0_z_0_yyy_xxzz[k] * ab_x + g_z_0_z_0_yyy_xxxzz[k];

                g_z_0_z_0_xyyy_xyyy[k] = -g_z_0_z_0_yyy_xyyy[k] * ab_x + g_z_0_z_0_yyy_xxyyy[k];

                g_z_0_z_0_xyyy_xyyz[k] = -g_z_0_z_0_yyy_xyyz[k] * ab_x + g_z_0_z_0_yyy_xxyyz[k];

                g_z_0_z_0_xyyy_xyzz[k] = -g_z_0_z_0_yyy_xyzz[k] * ab_x + g_z_0_z_0_yyy_xxyzz[k];

                g_z_0_z_0_xyyy_xzzz[k] = -g_z_0_z_0_yyy_xzzz[k] * ab_x + g_z_0_z_0_yyy_xxzzz[k];

                g_z_0_z_0_xyyy_yyyy[k] = -g_z_0_z_0_yyy_yyyy[k] * ab_x + g_z_0_z_0_yyy_xyyyy[k];

                g_z_0_z_0_xyyy_yyyz[k] = -g_z_0_z_0_yyy_yyyz[k] * ab_x + g_z_0_z_0_yyy_xyyyz[k];

                g_z_0_z_0_xyyy_yyzz[k] = -g_z_0_z_0_yyy_yyzz[k] * ab_x + g_z_0_z_0_yyy_xyyzz[k];

                g_z_0_z_0_xyyy_yzzz[k] = -g_z_0_z_0_yyy_yzzz[k] * ab_x + g_z_0_z_0_yyy_xyzzz[k];

                g_z_0_z_0_xyyy_zzzz[k] = -g_z_0_z_0_yyy_zzzz[k] * ab_x + g_z_0_z_0_yyy_xzzzz[k];
            }

            /// Set up 1905-1920 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1905 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1906 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1907 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1908 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1909 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1910 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1911 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1912 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1913 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1914 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1915 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1916 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1917 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1918 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1919 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyz_xxxx, g_z_0_z_0_xyyz_xxxy, g_z_0_z_0_xyyz_xxxz, g_z_0_z_0_xyyz_xxyy, g_z_0_z_0_xyyz_xxyz, g_z_0_z_0_xyyz_xxzz, g_z_0_z_0_xyyz_xyyy, g_z_0_z_0_xyyz_xyyz, g_z_0_z_0_xyyz_xyzz, g_z_0_z_0_xyyz_xzzz, g_z_0_z_0_xyyz_yyyy, g_z_0_z_0_xyyz_yyyz, g_z_0_z_0_xyyz_yyzz, g_z_0_z_0_xyyz_yzzz, g_z_0_z_0_xyyz_zzzz, g_z_0_z_0_yyz_xxxx, g_z_0_z_0_yyz_xxxxx, g_z_0_z_0_yyz_xxxxy, g_z_0_z_0_yyz_xxxxz, g_z_0_z_0_yyz_xxxy, g_z_0_z_0_yyz_xxxyy, g_z_0_z_0_yyz_xxxyz, g_z_0_z_0_yyz_xxxz, g_z_0_z_0_yyz_xxxzz, g_z_0_z_0_yyz_xxyy, g_z_0_z_0_yyz_xxyyy, g_z_0_z_0_yyz_xxyyz, g_z_0_z_0_yyz_xxyz, g_z_0_z_0_yyz_xxyzz, g_z_0_z_0_yyz_xxzz, g_z_0_z_0_yyz_xxzzz, g_z_0_z_0_yyz_xyyy, g_z_0_z_0_yyz_xyyyy, g_z_0_z_0_yyz_xyyyz, g_z_0_z_0_yyz_xyyz, g_z_0_z_0_yyz_xyyzz, g_z_0_z_0_yyz_xyzz, g_z_0_z_0_yyz_xyzzz, g_z_0_z_0_yyz_xzzz, g_z_0_z_0_yyz_xzzzz, g_z_0_z_0_yyz_yyyy, g_z_0_z_0_yyz_yyyz, g_z_0_z_0_yyz_yyzz, g_z_0_z_0_yyz_yzzz, g_z_0_z_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyz_xxxx[k] = -g_z_0_z_0_yyz_xxxx[k] * ab_x + g_z_0_z_0_yyz_xxxxx[k];

                g_z_0_z_0_xyyz_xxxy[k] = -g_z_0_z_0_yyz_xxxy[k] * ab_x + g_z_0_z_0_yyz_xxxxy[k];

                g_z_0_z_0_xyyz_xxxz[k] = -g_z_0_z_0_yyz_xxxz[k] * ab_x + g_z_0_z_0_yyz_xxxxz[k];

                g_z_0_z_0_xyyz_xxyy[k] = -g_z_0_z_0_yyz_xxyy[k] * ab_x + g_z_0_z_0_yyz_xxxyy[k];

                g_z_0_z_0_xyyz_xxyz[k] = -g_z_0_z_0_yyz_xxyz[k] * ab_x + g_z_0_z_0_yyz_xxxyz[k];

                g_z_0_z_0_xyyz_xxzz[k] = -g_z_0_z_0_yyz_xxzz[k] * ab_x + g_z_0_z_0_yyz_xxxzz[k];

                g_z_0_z_0_xyyz_xyyy[k] = -g_z_0_z_0_yyz_xyyy[k] * ab_x + g_z_0_z_0_yyz_xxyyy[k];

                g_z_0_z_0_xyyz_xyyz[k] = -g_z_0_z_0_yyz_xyyz[k] * ab_x + g_z_0_z_0_yyz_xxyyz[k];

                g_z_0_z_0_xyyz_xyzz[k] = -g_z_0_z_0_yyz_xyzz[k] * ab_x + g_z_0_z_0_yyz_xxyzz[k];

                g_z_0_z_0_xyyz_xzzz[k] = -g_z_0_z_0_yyz_xzzz[k] * ab_x + g_z_0_z_0_yyz_xxzzz[k];

                g_z_0_z_0_xyyz_yyyy[k] = -g_z_0_z_0_yyz_yyyy[k] * ab_x + g_z_0_z_0_yyz_xyyyy[k];

                g_z_0_z_0_xyyz_yyyz[k] = -g_z_0_z_0_yyz_yyyz[k] * ab_x + g_z_0_z_0_yyz_xyyyz[k];

                g_z_0_z_0_xyyz_yyzz[k] = -g_z_0_z_0_yyz_yyzz[k] * ab_x + g_z_0_z_0_yyz_xyyzz[k];

                g_z_0_z_0_xyyz_yzzz[k] = -g_z_0_z_0_yyz_yzzz[k] * ab_x + g_z_0_z_0_yyz_xyzzz[k];

                g_z_0_z_0_xyyz_zzzz[k] = -g_z_0_z_0_yyz_zzzz[k] * ab_x + g_z_0_z_0_yyz_xzzzz[k];
            }

            /// Set up 1920-1935 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1920 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1921 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1922 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1923 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1924 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1925 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1926 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1927 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1928 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1929 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1930 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1931 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1932 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1933 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1934 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyzz_xxxx, g_z_0_z_0_xyzz_xxxy, g_z_0_z_0_xyzz_xxxz, g_z_0_z_0_xyzz_xxyy, g_z_0_z_0_xyzz_xxyz, g_z_0_z_0_xyzz_xxzz, g_z_0_z_0_xyzz_xyyy, g_z_0_z_0_xyzz_xyyz, g_z_0_z_0_xyzz_xyzz, g_z_0_z_0_xyzz_xzzz, g_z_0_z_0_xyzz_yyyy, g_z_0_z_0_xyzz_yyyz, g_z_0_z_0_xyzz_yyzz, g_z_0_z_0_xyzz_yzzz, g_z_0_z_0_xyzz_zzzz, g_z_0_z_0_yzz_xxxx, g_z_0_z_0_yzz_xxxxx, g_z_0_z_0_yzz_xxxxy, g_z_0_z_0_yzz_xxxxz, g_z_0_z_0_yzz_xxxy, g_z_0_z_0_yzz_xxxyy, g_z_0_z_0_yzz_xxxyz, g_z_0_z_0_yzz_xxxz, g_z_0_z_0_yzz_xxxzz, g_z_0_z_0_yzz_xxyy, g_z_0_z_0_yzz_xxyyy, g_z_0_z_0_yzz_xxyyz, g_z_0_z_0_yzz_xxyz, g_z_0_z_0_yzz_xxyzz, g_z_0_z_0_yzz_xxzz, g_z_0_z_0_yzz_xxzzz, g_z_0_z_0_yzz_xyyy, g_z_0_z_0_yzz_xyyyy, g_z_0_z_0_yzz_xyyyz, g_z_0_z_0_yzz_xyyz, g_z_0_z_0_yzz_xyyzz, g_z_0_z_0_yzz_xyzz, g_z_0_z_0_yzz_xyzzz, g_z_0_z_0_yzz_xzzz, g_z_0_z_0_yzz_xzzzz, g_z_0_z_0_yzz_yyyy, g_z_0_z_0_yzz_yyyz, g_z_0_z_0_yzz_yyzz, g_z_0_z_0_yzz_yzzz, g_z_0_z_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyzz_xxxx[k] = -g_z_0_z_0_yzz_xxxx[k] * ab_x + g_z_0_z_0_yzz_xxxxx[k];

                g_z_0_z_0_xyzz_xxxy[k] = -g_z_0_z_0_yzz_xxxy[k] * ab_x + g_z_0_z_0_yzz_xxxxy[k];

                g_z_0_z_0_xyzz_xxxz[k] = -g_z_0_z_0_yzz_xxxz[k] * ab_x + g_z_0_z_0_yzz_xxxxz[k];

                g_z_0_z_0_xyzz_xxyy[k] = -g_z_0_z_0_yzz_xxyy[k] * ab_x + g_z_0_z_0_yzz_xxxyy[k];

                g_z_0_z_0_xyzz_xxyz[k] = -g_z_0_z_0_yzz_xxyz[k] * ab_x + g_z_0_z_0_yzz_xxxyz[k];

                g_z_0_z_0_xyzz_xxzz[k] = -g_z_0_z_0_yzz_xxzz[k] * ab_x + g_z_0_z_0_yzz_xxxzz[k];

                g_z_0_z_0_xyzz_xyyy[k] = -g_z_0_z_0_yzz_xyyy[k] * ab_x + g_z_0_z_0_yzz_xxyyy[k];

                g_z_0_z_0_xyzz_xyyz[k] = -g_z_0_z_0_yzz_xyyz[k] * ab_x + g_z_0_z_0_yzz_xxyyz[k];

                g_z_0_z_0_xyzz_xyzz[k] = -g_z_0_z_0_yzz_xyzz[k] * ab_x + g_z_0_z_0_yzz_xxyzz[k];

                g_z_0_z_0_xyzz_xzzz[k] = -g_z_0_z_0_yzz_xzzz[k] * ab_x + g_z_0_z_0_yzz_xxzzz[k];

                g_z_0_z_0_xyzz_yyyy[k] = -g_z_0_z_0_yzz_yyyy[k] * ab_x + g_z_0_z_0_yzz_xyyyy[k];

                g_z_0_z_0_xyzz_yyyz[k] = -g_z_0_z_0_yzz_yyyz[k] * ab_x + g_z_0_z_0_yzz_xyyyz[k];

                g_z_0_z_0_xyzz_yyzz[k] = -g_z_0_z_0_yzz_yyzz[k] * ab_x + g_z_0_z_0_yzz_xyyzz[k];

                g_z_0_z_0_xyzz_yzzz[k] = -g_z_0_z_0_yzz_yzzz[k] * ab_x + g_z_0_z_0_yzz_xyzzz[k];

                g_z_0_z_0_xyzz_zzzz[k] = -g_z_0_z_0_yzz_zzzz[k] * ab_x + g_z_0_z_0_yzz_xzzzz[k];
            }

            /// Set up 1935-1950 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1935 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1936 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1937 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1938 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1939 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_xxzz = cbuffer.data(gg_geom_1010_off + 1940 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_xyyy = cbuffer.data(gg_geom_1010_off + 1941 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_xyyz = cbuffer.data(gg_geom_1010_off + 1942 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_xyzz = cbuffer.data(gg_geom_1010_off + 1943 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_xzzz = cbuffer.data(gg_geom_1010_off + 1944 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_yyyy = cbuffer.data(gg_geom_1010_off + 1945 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_yyyz = cbuffer.data(gg_geom_1010_off + 1946 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_yyzz = cbuffer.data(gg_geom_1010_off + 1947 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_yzzz = cbuffer.data(gg_geom_1010_off + 1948 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_zzzz = cbuffer.data(gg_geom_1010_off + 1949 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xzzz_xxxx, g_z_0_z_0_xzzz_xxxy, g_z_0_z_0_xzzz_xxxz, g_z_0_z_0_xzzz_xxyy, g_z_0_z_0_xzzz_xxyz, g_z_0_z_0_xzzz_xxzz, g_z_0_z_0_xzzz_xyyy, g_z_0_z_0_xzzz_xyyz, g_z_0_z_0_xzzz_xyzz, g_z_0_z_0_xzzz_xzzz, g_z_0_z_0_xzzz_yyyy, g_z_0_z_0_xzzz_yyyz, g_z_0_z_0_xzzz_yyzz, g_z_0_z_0_xzzz_yzzz, g_z_0_z_0_xzzz_zzzz, g_z_0_z_0_zzz_xxxx, g_z_0_z_0_zzz_xxxxx, g_z_0_z_0_zzz_xxxxy, g_z_0_z_0_zzz_xxxxz, g_z_0_z_0_zzz_xxxy, g_z_0_z_0_zzz_xxxyy, g_z_0_z_0_zzz_xxxyz, g_z_0_z_0_zzz_xxxz, g_z_0_z_0_zzz_xxxzz, g_z_0_z_0_zzz_xxyy, g_z_0_z_0_zzz_xxyyy, g_z_0_z_0_zzz_xxyyz, g_z_0_z_0_zzz_xxyz, g_z_0_z_0_zzz_xxyzz, g_z_0_z_0_zzz_xxzz, g_z_0_z_0_zzz_xxzzz, g_z_0_z_0_zzz_xyyy, g_z_0_z_0_zzz_xyyyy, g_z_0_z_0_zzz_xyyyz, g_z_0_z_0_zzz_xyyz, g_z_0_z_0_zzz_xyyzz, g_z_0_z_0_zzz_xyzz, g_z_0_z_0_zzz_xyzzz, g_z_0_z_0_zzz_xzzz, g_z_0_z_0_zzz_xzzzz, g_z_0_z_0_zzz_yyyy, g_z_0_z_0_zzz_yyyz, g_z_0_z_0_zzz_yyzz, g_z_0_z_0_zzz_yzzz, g_z_0_z_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xzzz_xxxx[k] = -g_z_0_z_0_zzz_xxxx[k] * ab_x + g_z_0_z_0_zzz_xxxxx[k];

                g_z_0_z_0_xzzz_xxxy[k] = -g_z_0_z_0_zzz_xxxy[k] * ab_x + g_z_0_z_0_zzz_xxxxy[k];

                g_z_0_z_0_xzzz_xxxz[k] = -g_z_0_z_0_zzz_xxxz[k] * ab_x + g_z_0_z_0_zzz_xxxxz[k];

                g_z_0_z_0_xzzz_xxyy[k] = -g_z_0_z_0_zzz_xxyy[k] * ab_x + g_z_0_z_0_zzz_xxxyy[k];

                g_z_0_z_0_xzzz_xxyz[k] = -g_z_0_z_0_zzz_xxyz[k] * ab_x + g_z_0_z_0_zzz_xxxyz[k];

                g_z_0_z_0_xzzz_xxzz[k] = -g_z_0_z_0_zzz_xxzz[k] * ab_x + g_z_0_z_0_zzz_xxxzz[k];

                g_z_0_z_0_xzzz_xyyy[k] = -g_z_0_z_0_zzz_xyyy[k] * ab_x + g_z_0_z_0_zzz_xxyyy[k];

                g_z_0_z_0_xzzz_xyyz[k] = -g_z_0_z_0_zzz_xyyz[k] * ab_x + g_z_0_z_0_zzz_xxyyz[k];

                g_z_0_z_0_xzzz_xyzz[k] = -g_z_0_z_0_zzz_xyzz[k] * ab_x + g_z_0_z_0_zzz_xxyzz[k];

                g_z_0_z_0_xzzz_xzzz[k] = -g_z_0_z_0_zzz_xzzz[k] * ab_x + g_z_0_z_0_zzz_xxzzz[k];

                g_z_0_z_0_xzzz_yyyy[k] = -g_z_0_z_0_zzz_yyyy[k] * ab_x + g_z_0_z_0_zzz_xyyyy[k];

                g_z_0_z_0_xzzz_yyyz[k] = -g_z_0_z_0_zzz_yyyz[k] * ab_x + g_z_0_z_0_zzz_xyyyz[k];

                g_z_0_z_0_xzzz_yyzz[k] = -g_z_0_z_0_zzz_yyzz[k] * ab_x + g_z_0_z_0_zzz_xyyzz[k];

                g_z_0_z_0_xzzz_yzzz[k] = -g_z_0_z_0_zzz_yzzz[k] * ab_x + g_z_0_z_0_zzz_xyzzz[k];

                g_z_0_z_0_xzzz_zzzz[k] = -g_z_0_z_0_zzz_zzzz[k] * ab_x + g_z_0_z_0_zzz_xzzzz[k];
            }

            /// Set up 1950-1965 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyy_xxxx = cbuffer.data(gg_geom_1010_off + 1950 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_xxxy = cbuffer.data(gg_geom_1010_off + 1951 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_xxxz = cbuffer.data(gg_geom_1010_off + 1952 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_xxyy = cbuffer.data(gg_geom_1010_off + 1953 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_xxyz = cbuffer.data(gg_geom_1010_off + 1954 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_xxzz = cbuffer.data(gg_geom_1010_off + 1955 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_xyyy = cbuffer.data(gg_geom_1010_off + 1956 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_xyyz = cbuffer.data(gg_geom_1010_off + 1957 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_xyzz = cbuffer.data(gg_geom_1010_off + 1958 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_xzzz = cbuffer.data(gg_geom_1010_off + 1959 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_yyyy = cbuffer.data(gg_geom_1010_off + 1960 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_yyyz = cbuffer.data(gg_geom_1010_off + 1961 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_yyzz = cbuffer.data(gg_geom_1010_off + 1962 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_yzzz = cbuffer.data(gg_geom_1010_off + 1963 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_zzzz = cbuffer.data(gg_geom_1010_off + 1964 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyy_xxxx, g_z_0_z_0_yyy_xxxxy, g_z_0_z_0_yyy_xxxy, g_z_0_z_0_yyy_xxxyy, g_z_0_z_0_yyy_xxxyz, g_z_0_z_0_yyy_xxxz, g_z_0_z_0_yyy_xxyy, g_z_0_z_0_yyy_xxyyy, g_z_0_z_0_yyy_xxyyz, g_z_0_z_0_yyy_xxyz, g_z_0_z_0_yyy_xxyzz, g_z_0_z_0_yyy_xxzz, g_z_0_z_0_yyy_xyyy, g_z_0_z_0_yyy_xyyyy, g_z_0_z_0_yyy_xyyyz, g_z_0_z_0_yyy_xyyz, g_z_0_z_0_yyy_xyyzz, g_z_0_z_0_yyy_xyzz, g_z_0_z_0_yyy_xyzzz, g_z_0_z_0_yyy_xzzz, g_z_0_z_0_yyy_yyyy, g_z_0_z_0_yyy_yyyyy, g_z_0_z_0_yyy_yyyyz, g_z_0_z_0_yyy_yyyz, g_z_0_z_0_yyy_yyyzz, g_z_0_z_0_yyy_yyzz, g_z_0_z_0_yyy_yyzzz, g_z_0_z_0_yyy_yzzz, g_z_0_z_0_yyy_yzzzz, g_z_0_z_0_yyy_zzzz, g_z_0_z_0_yyyy_xxxx, g_z_0_z_0_yyyy_xxxy, g_z_0_z_0_yyyy_xxxz, g_z_0_z_0_yyyy_xxyy, g_z_0_z_0_yyyy_xxyz, g_z_0_z_0_yyyy_xxzz, g_z_0_z_0_yyyy_xyyy, g_z_0_z_0_yyyy_xyyz, g_z_0_z_0_yyyy_xyzz, g_z_0_z_0_yyyy_xzzz, g_z_0_z_0_yyyy_yyyy, g_z_0_z_0_yyyy_yyyz, g_z_0_z_0_yyyy_yyzz, g_z_0_z_0_yyyy_yzzz, g_z_0_z_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyy_xxxx[k] = -g_z_0_z_0_yyy_xxxx[k] * ab_y + g_z_0_z_0_yyy_xxxxy[k];

                g_z_0_z_0_yyyy_xxxy[k] = -g_z_0_z_0_yyy_xxxy[k] * ab_y + g_z_0_z_0_yyy_xxxyy[k];

                g_z_0_z_0_yyyy_xxxz[k] = -g_z_0_z_0_yyy_xxxz[k] * ab_y + g_z_0_z_0_yyy_xxxyz[k];

                g_z_0_z_0_yyyy_xxyy[k] = -g_z_0_z_0_yyy_xxyy[k] * ab_y + g_z_0_z_0_yyy_xxyyy[k];

                g_z_0_z_0_yyyy_xxyz[k] = -g_z_0_z_0_yyy_xxyz[k] * ab_y + g_z_0_z_0_yyy_xxyyz[k];

                g_z_0_z_0_yyyy_xxzz[k] = -g_z_0_z_0_yyy_xxzz[k] * ab_y + g_z_0_z_0_yyy_xxyzz[k];

                g_z_0_z_0_yyyy_xyyy[k] = -g_z_0_z_0_yyy_xyyy[k] * ab_y + g_z_0_z_0_yyy_xyyyy[k];

                g_z_0_z_0_yyyy_xyyz[k] = -g_z_0_z_0_yyy_xyyz[k] * ab_y + g_z_0_z_0_yyy_xyyyz[k];

                g_z_0_z_0_yyyy_xyzz[k] = -g_z_0_z_0_yyy_xyzz[k] * ab_y + g_z_0_z_0_yyy_xyyzz[k];

                g_z_0_z_0_yyyy_xzzz[k] = -g_z_0_z_0_yyy_xzzz[k] * ab_y + g_z_0_z_0_yyy_xyzzz[k];

                g_z_0_z_0_yyyy_yyyy[k] = -g_z_0_z_0_yyy_yyyy[k] * ab_y + g_z_0_z_0_yyy_yyyyy[k];

                g_z_0_z_0_yyyy_yyyz[k] = -g_z_0_z_0_yyy_yyyz[k] * ab_y + g_z_0_z_0_yyy_yyyyz[k];

                g_z_0_z_0_yyyy_yyzz[k] = -g_z_0_z_0_yyy_yyzz[k] * ab_y + g_z_0_z_0_yyy_yyyzz[k];

                g_z_0_z_0_yyyy_yzzz[k] = -g_z_0_z_0_yyy_yzzz[k] * ab_y + g_z_0_z_0_yyy_yyzzz[k];

                g_z_0_z_0_yyyy_zzzz[k] = -g_z_0_z_0_yyy_zzzz[k] * ab_y + g_z_0_z_0_yyy_yzzzz[k];
            }

            /// Set up 1965-1980 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyz_xxxx = cbuffer.data(gg_geom_1010_off + 1965 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_xxxy = cbuffer.data(gg_geom_1010_off + 1966 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_xxxz = cbuffer.data(gg_geom_1010_off + 1967 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_xxyy = cbuffer.data(gg_geom_1010_off + 1968 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_xxyz = cbuffer.data(gg_geom_1010_off + 1969 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_xxzz = cbuffer.data(gg_geom_1010_off + 1970 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_xyyy = cbuffer.data(gg_geom_1010_off + 1971 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_xyyz = cbuffer.data(gg_geom_1010_off + 1972 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_xyzz = cbuffer.data(gg_geom_1010_off + 1973 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_xzzz = cbuffer.data(gg_geom_1010_off + 1974 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_yyyy = cbuffer.data(gg_geom_1010_off + 1975 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_yyyz = cbuffer.data(gg_geom_1010_off + 1976 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_yyzz = cbuffer.data(gg_geom_1010_off + 1977 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_yzzz = cbuffer.data(gg_geom_1010_off + 1978 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_zzzz = cbuffer.data(gg_geom_1010_off + 1979 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyz_xxxx, g_z_0_z_0_yyyz_xxxy, g_z_0_z_0_yyyz_xxxz, g_z_0_z_0_yyyz_xxyy, g_z_0_z_0_yyyz_xxyz, g_z_0_z_0_yyyz_xxzz, g_z_0_z_0_yyyz_xyyy, g_z_0_z_0_yyyz_xyyz, g_z_0_z_0_yyyz_xyzz, g_z_0_z_0_yyyz_xzzz, g_z_0_z_0_yyyz_yyyy, g_z_0_z_0_yyyz_yyyz, g_z_0_z_0_yyyz_yyzz, g_z_0_z_0_yyyz_yzzz, g_z_0_z_0_yyyz_zzzz, g_z_0_z_0_yyz_xxxx, g_z_0_z_0_yyz_xxxxy, g_z_0_z_0_yyz_xxxy, g_z_0_z_0_yyz_xxxyy, g_z_0_z_0_yyz_xxxyz, g_z_0_z_0_yyz_xxxz, g_z_0_z_0_yyz_xxyy, g_z_0_z_0_yyz_xxyyy, g_z_0_z_0_yyz_xxyyz, g_z_0_z_0_yyz_xxyz, g_z_0_z_0_yyz_xxyzz, g_z_0_z_0_yyz_xxzz, g_z_0_z_0_yyz_xyyy, g_z_0_z_0_yyz_xyyyy, g_z_0_z_0_yyz_xyyyz, g_z_0_z_0_yyz_xyyz, g_z_0_z_0_yyz_xyyzz, g_z_0_z_0_yyz_xyzz, g_z_0_z_0_yyz_xyzzz, g_z_0_z_0_yyz_xzzz, g_z_0_z_0_yyz_yyyy, g_z_0_z_0_yyz_yyyyy, g_z_0_z_0_yyz_yyyyz, g_z_0_z_0_yyz_yyyz, g_z_0_z_0_yyz_yyyzz, g_z_0_z_0_yyz_yyzz, g_z_0_z_0_yyz_yyzzz, g_z_0_z_0_yyz_yzzz, g_z_0_z_0_yyz_yzzzz, g_z_0_z_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyz_xxxx[k] = -g_z_0_z_0_yyz_xxxx[k] * ab_y + g_z_0_z_0_yyz_xxxxy[k];

                g_z_0_z_0_yyyz_xxxy[k] = -g_z_0_z_0_yyz_xxxy[k] * ab_y + g_z_0_z_0_yyz_xxxyy[k];

                g_z_0_z_0_yyyz_xxxz[k] = -g_z_0_z_0_yyz_xxxz[k] * ab_y + g_z_0_z_0_yyz_xxxyz[k];

                g_z_0_z_0_yyyz_xxyy[k] = -g_z_0_z_0_yyz_xxyy[k] * ab_y + g_z_0_z_0_yyz_xxyyy[k];

                g_z_0_z_0_yyyz_xxyz[k] = -g_z_0_z_0_yyz_xxyz[k] * ab_y + g_z_0_z_0_yyz_xxyyz[k];

                g_z_0_z_0_yyyz_xxzz[k] = -g_z_0_z_0_yyz_xxzz[k] * ab_y + g_z_0_z_0_yyz_xxyzz[k];

                g_z_0_z_0_yyyz_xyyy[k] = -g_z_0_z_0_yyz_xyyy[k] * ab_y + g_z_0_z_0_yyz_xyyyy[k];

                g_z_0_z_0_yyyz_xyyz[k] = -g_z_0_z_0_yyz_xyyz[k] * ab_y + g_z_0_z_0_yyz_xyyyz[k];

                g_z_0_z_0_yyyz_xyzz[k] = -g_z_0_z_0_yyz_xyzz[k] * ab_y + g_z_0_z_0_yyz_xyyzz[k];

                g_z_0_z_0_yyyz_xzzz[k] = -g_z_0_z_0_yyz_xzzz[k] * ab_y + g_z_0_z_0_yyz_xyzzz[k];

                g_z_0_z_0_yyyz_yyyy[k] = -g_z_0_z_0_yyz_yyyy[k] * ab_y + g_z_0_z_0_yyz_yyyyy[k];

                g_z_0_z_0_yyyz_yyyz[k] = -g_z_0_z_0_yyz_yyyz[k] * ab_y + g_z_0_z_0_yyz_yyyyz[k];

                g_z_0_z_0_yyyz_yyzz[k] = -g_z_0_z_0_yyz_yyzz[k] * ab_y + g_z_0_z_0_yyz_yyyzz[k];

                g_z_0_z_0_yyyz_yzzz[k] = -g_z_0_z_0_yyz_yzzz[k] * ab_y + g_z_0_z_0_yyz_yyzzz[k];

                g_z_0_z_0_yyyz_zzzz[k] = -g_z_0_z_0_yyz_zzzz[k] * ab_y + g_z_0_z_0_yyz_yzzzz[k];
            }

            /// Set up 1980-1995 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyzz_xxxx = cbuffer.data(gg_geom_1010_off + 1980 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_xxxy = cbuffer.data(gg_geom_1010_off + 1981 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_xxxz = cbuffer.data(gg_geom_1010_off + 1982 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_xxyy = cbuffer.data(gg_geom_1010_off + 1983 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_xxyz = cbuffer.data(gg_geom_1010_off + 1984 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_xxzz = cbuffer.data(gg_geom_1010_off + 1985 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_xyyy = cbuffer.data(gg_geom_1010_off + 1986 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_xyyz = cbuffer.data(gg_geom_1010_off + 1987 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_xyzz = cbuffer.data(gg_geom_1010_off + 1988 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_xzzz = cbuffer.data(gg_geom_1010_off + 1989 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_yyyy = cbuffer.data(gg_geom_1010_off + 1990 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_yyyz = cbuffer.data(gg_geom_1010_off + 1991 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_yyzz = cbuffer.data(gg_geom_1010_off + 1992 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_yzzz = cbuffer.data(gg_geom_1010_off + 1993 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_zzzz = cbuffer.data(gg_geom_1010_off + 1994 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyzz_xxxx, g_z_0_z_0_yyzz_xxxy, g_z_0_z_0_yyzz_xxxz, g_z_0_z_0_yyzz_xxyy, g_z_0_z_0_yyzz_xxyz, g_z_0_z_0_yyzz_xxzz, g_z_0_z_0_yyzz_xyyy, g_z_0_z_0_yyzz_xyyz, g_z_0_z_0_yyzz_xyzz, g_z_0_z_0_yyzz_xzzz, g_z_0_z_0_yyzz_yyyy, g_z_0_z_0_yyzz_yyyz, g_z_0_z_0_yyzz_yyzz, g_z_0_z_0_yyzz_yzzz, g_z_0_z_0_yyzz_zzzz, g_z_0_z_0_yzz_xxxx, g_z_0_z_0_yzz_xxxxy, g_z_0_z_0_yzz_xxxy, g_z_0_z_0_yzz_xxxyy, g_z_0_z_0_yzz_xxxyz, g_z_0_z_0_yzz_xxxz, g_z_0_z_0_yzz_xxyy, g_z_0_z_0_yzz_xxyyy, g_z_0_z_0_yzz_xxyyz, g_z_0_z_0_yzz_xxyz, g_z_0_z_0_yzz_xxyzz, g_z_0_z_0_yzz_xxzz, g_z_0_z_0_yzz_xyyy, g_z_0_z_0_yzz_xyyyy, g_z_0_z_0_yzz_xyyyz, g_z_0_z_0_yzz_xyyz, g_z_0_z_0_yzz_xyyzz, g_z_0_z_0_yzz_xyzz, g_z_0_z_0_yzz_xyzzz, g_z_0_z_0_yzz_xzzz, g_z_0_z_0_yzz_yyyy, g_z_0_z_0_yzz_yyyyy, g_z_0_z_0_yzz_yyyyz, g_z_0_z_0_yzz_yyyz, g_z_0_z_0_yzz_yyyzz, g_z_0_z_0_yzz_yyzz, g_z_0_z_0_yzz_yyzzz, g_z_0_z_0_yzz_yzzz, g_z_0_z_0_yzz_yzzzz, g_z_0_z_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyzz_xxxx[k] = -g_z_0_z_0_yzz_xxxx[k] * ab_y + g_z_0_z_0_yzz_xxxxy[k];

                g_z_0_z_0_yyzz_xxxy[k] = -g_z_0_z_0_yzz_xxxy[k] * ab_y + g_z_0_z_0_yzz_xxxyy[k];

                g_z_0_z_0_yyzz_xxxz[k] = -g_z_0_z_0_yzz_xxxz[k] * ab_y + g_z_0_z_0_yzz_xxxyz[k];

                g_z_0_z_0_yyzz_xxyy[k] = -g_z_0_z_0_yzz_xxyy[k] * ab_y + g_z_0_z_0_yzz_xxyyy[k];

                g_z_0_z_0_yyzz_xxyz[k] = -g_z_0_z_0_yzz_xxyz[k] * ab_y + g_z_0_z_0_yzz_xxyyz[k];

                g_z_0_z_0_yyzz_xxzz[k] = -g_z_0_z_0_yzz_xxzz[k] * ab_y + g_z_0_z_0_yzz_xxyzz[k];

                g_z_0_z_0_yyzz_xyyy[k] = -g_z_0_z_0_yzz_xyyy[k] * ab_y + g_z_0_z_0_yzz_xyyyy[k];

                g_z_0_z_0_yyzz_xyyz[k] = -g_z_0_z_0_yzz_xyyz[k] * ab_y + g_z_0_z_0_yzz_xyyyz[k];

                g_z_0_z_0_yyzz_xyzz[k] = -g_z_0_z_0_yzz_xyzz[k] * ab_y + g_z_0_z_0_yzz_xyyzz[k];

                g_z_0_z_0_yyzz_xzzz[k] = -g_z_0_z_0_yzz_xzzz[k] * ab_y + g_z_0_z_0_yzz_xyzzz[k];

                g_z_0_z_0_yyzz_yyyy[k] = -g_z_0_z_0_yzz_yyyy[k] * ab_y + g_z_0_z_0_yzz_yyyyy[k];

                g_z_0_z_0_yyzz_yyyz[k] = -g_z_0_z_0_yzz_yyyz[k] * ab_y + g_z_0_z_0_yzz_yyyyz[k];

                g_z_0_z_0_yyzz_yyzz[k] = -g_z_0_z_0_yzz_yyzz[k] * ab_y + g_z_0_z_0_yzz_yyyzz[k];

                g_z_0_z_0_yyzz_yzzz[k] = -g_z_0_z_0_yzz_yzzz[k] * ab_y + g_z_0_z_0_yzz_yyzzz[k];

                g_z_0_z_0_yyzz_zzzz[k] = -g_z_0_z_0_yzz_zzzz[k] * ab_y + g_z_0_z_0_yzz_yzzzz[k];
            }

            /// Set up 1995-2010 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yzzz_xxxx = cbuffer.data(gg_geom_1010_off + 1995 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_xxxy = cbuffer.data(gg_geom_1010_off + 1996 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_xxxz = cbuffer.data(gg_geom_1010_off + 1997 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_xxyy = cbuffer.data(gg_geom_1010_off + 1998 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_xxyz = cbuffer.data(gg_geom_1010_off + 1999 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_xxzz = cbuffer.data(gg_geom_1010_off + 2000 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_xyyy = cbuffer.data(gg_geom_1010_off + 2001 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_xyyz = cbuffer.data(gg_geom_1010_off + 2002 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_xyzz = cbuffer.data(gg_geom_1010_off + 2003 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_xzzz = cbuffer.data(gg_geom_1010_off + 2004 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_yyyy = cbuffer.data(gg_geom_1010_off + 2005 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_yyyz = cbuffer.data(gg_geom_1010_off + 2006 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_yyzz = cbuffer.data(gg_geom_1010_off + 2007 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_yzzz = cbuffer.data(gg_geom_1010_off + 2008 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_zzzz = cbuffer.data(gg_geom_1010_off + 2009 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yzzz_xxxx, g_z_0_z_0_yzzz_xxxy, g_z_0_z_0_yzzz_xxxz, g_z_0_z_0_yzzz_xxyy, g_z_0_z_0_yzzz_xxyz, g_z_0_z_0_yzzz_xxzz, g_z_0_z_0_yzzz_xyyy, g_z_0_z_0_yzzz_xyyz, g_z_0_z_0_yzzz_xyzz, g_z_0_z_0_yzzz_xzzz, g_z_0_z_0_yzzz_yyyy, g_z_0_z_0_yzzz_yyyz, g_z_0_z_0_yzzz_yyzz, g_z_0_z_0_yzzz_yzzz, g_z_0_z_0_yzzz_zzzz, g_z_0_z_0_zzz_xxxx, g_z_0_z_0_zzz_xxxxy, g_z_0_z_0_zzz_xxxy, g_z_0_z_0_zzz_xxxyy, g_z_0_z_0_zzz_xxxyz, g_z_0_z_0_zzz_xxxz, g_z_0_z_0_zzz_xxyy, g_z_0_z_0_zzz_xxyyy, g_z_0_z_0_zzz_xxyyz, g_z_0_z_0_zzz_xxyz, g_z_0_z_0_zzz_xxyzz, g_z_0_z_0_zzz_xxzz, g_z_0_z_0_zzz_xyyy, g_z_0_z_0_zzz_xyyyy, g_z_0_z_0_zzz_xyyyz, g_z_0_z_0_zzz_xyyz, g_z_0_z_0_zzz_xyyzz, g_z_0_z_0_zzz_xyzz, g_z_0_z_0_zzz_xyzzz, g_z_0_z_0_zzz_xzzz, g_z_0_z_0_zzz_yyyy, g_z_0_z_0_zzz_yyyyy, g_z_0_z_0_zzz_yyyyz, g_z_0_z_0_zzz_yyyz, g_z_0_z_0_zzz_yyyzz, g_z_0_z_0_zzz_yyzz, g_z_0_z_0_zzz_yyzzz, g_z_0_z_0_zzz_yzzz, g_z_0_z_0_zzz_yzzzz, g_z_0_z_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yzzz_xxxx[k] = -g_z_0_z_0_zzz_xxxx[k] * ab_y + g_z_0_z_0_zzz_xxxxy[k];

                g_z_0_z_0_yzzz_xxxy[k] = -g_z_0_z_0_zzz_xxxy[k] * ab_y + g_z_0_z_0_zzz_xxxyy[k];

                g_z_0_z_0_yzzz_xxxz[k] = -g_z_0_z_0_zzz_xxxz[k] * ab_y + g_z_0_z_0_zzz_xxxyz[k];

                g_z_0_z_0_yzzz_xxyy[k] = -g_z_0_z_0_zzz_xxyy[k] * ab_y + g_z_0_z_0_zzz_xxyyy[k];

                g_z_0_z_0_yzzz_xxyz[k] = -g_z_0_z_0_zzz_xxyz[k] * ab_y + g_z_0_z_0_zzz_xxyyz[k];

                g_z_0_z_0_yzzz_xxzz[k] = -g_z_0_z_0_zzz_xxzz[k] * ab_y + g_z_0_z_0_zzz_xxyzz[k];

                g_z_0_z_0_yzzz_xyyy[k] = -g_z_0_z_0_zzz_xyyy[k] * ab_y + g_z_0_z_0_zzz_xyyyy[k];

                g_z_0_z_0_yzzz_xyyz[k] = -g_z_0_z_0_zzz_xyyz[k] * ab_y + g_z_0_z_0_zzz_xyyyz[k];

                g_z_0_z_0_yzzz_xyzz[k] = -g_z_0_z_0_zzz_xyzz[k] * ab_y + g_z_0_z_0_zzz_xyyzz[k];

                g_z_0_z_0_yzzz_xzzz[k] = -g_z_0_z_0_zzz_xzzz[k] * ab_y + g_z_0_z_0_zzz_xyzzz[k];

                g_z_0_z_0_yzzz_yyyy[k] = -g_z_0_z_0_zzz_yyyy[k] * ab_y + g_z_0_z_0_zzz_yyyyy[k];

                g_z_0_z_0_yzzz_yyyz[k] = -g_z_0_z_0_zzz_yyyz[k] * ab_y + g_z_0_z_0_zzz_yyyyz[k];

                g_z_0_z_0_yzzz_yyzz[k] = -g_z_0_z_0_zzz_yyzz[k] * ab_y + g_z_0_z_0_zzz_yyyzz[k];

                g_z_0_z_0_yzzz_yzzz[k] = -g_z_0_z_0_zzz_yzzz[k] * ab_y + g_z_0_z_0_zzz_yyzzz[k];

                g_z_0_z_0_yzzz_zzzz[k] = -g_z_0_z_0_zzz_zzzz[k] * ab_y + g_z_0_z_0_zzz_yzzzz[k];
            }

            /// Set up 2010-2025 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_zzzz_xxxx = cbuffer.data(gg_geom_1010_off + 2010 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_xxxy = cbuffer.data(gg_geom_1010_off + 2011 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_xxxz = cbuffer.data(gg_geom_1010_off + 2012 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_xxyy = cbuffer.data(gg_geom_1010_off + 2013 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_xxyz = cbuffer.data(gg_geom_1010_off + 2014 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_xxzz = cbuffer.data(gg_geom_1010_off + 2015 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_xyyy = cbuffer.data(gg_geom_1010_off + 2016 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_xyyz = cbuffer.data(gg_geom_1010_off + 2017 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_xyzz = cbuffer.data(gg_geom_1010_off + 2018 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_xzzz = cbuffer.data(gg_geom_1010_off + 2019 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_yyyy = cbuffer.data(gg_geom_1010_off + 2020 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_yyyz = cbuffer.data(gg_geom_1010_off + 2021 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_yyzz = cbuffer.data(gg_geom_1010_off + 2022 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_yzzz = cbuffer.data(gg_geom_1010_off + 2023 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_zzzz = cbuffer.data(gg_geom_1010_off + 2024 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_zzz_xxxx, g_0_0_z_0_zzz_xxxy, g_0_0_z_0_zzz_xxxz, g_0_0_z_0_zzz_xxyy, g_0_0_z_0_zzz_xxyz, g_0_0_z_0_zzz_xxzz, g_0_0_z_0_zzz_xyyy, g_0_0_z_0_zzz_xyyz, g_0_0_z_0_zzz_xyzz, g_0_0_z_0_zzz_xzzz, g_0_0_z_0_zzz_yyyy, g_0_0_z_0_zzz_yyyz, g_0_0_z_0_zzz_yyzz, g_0_0_z_0_zzz_yzzz, g_0_0_z_0_zzz_zzzz, g_z_0_z_0_zzz_xxxx, g_z_0_z_0_zzz_xxxxz, g_z_0_z_0_zzz_xxxy, g_z_0_z_0_zzz_xxxyz, g_z_0_z_0_zzz_xxxz, g_z_0_z_0_zzz_xxxzz, g_z_0_z_0_zzz_xxyy, g_z_0_z_0_zzz_xxyyz, g_z_0_z_0_zzz_xxyz, g_z_0_z_0_zzz_xxyzz, g_z_0_z_0_zzz_xxzz, g_z_0_z_0_zzz_xxzzz, g_z_0_z_0_zzz_xyyy, g_z_0_z_0_zzz_xyyyz, g_z_0_z_0_zzz_xyyz, g_z_0_z_0_zzz_xyyzz, g_z_0_z_0_zzz_xyzz, g_z_0_z_0_zzz_xyzzz, g_z_0_z_0_zzz_xzzz, g_z_0_z_0_zzz_xzzzz, g_z_0_z_0_zzz_yyyy, g_z_0_z_0_zzz_yyyyz, g_z_0_z_0_zzz_yyyz, g_z_0_z_0_zzz_yyyzz, g_z_0_z_0_zzz_yyzz, g_z_0_z_0_zzz_yyzzz, g_z_0_z_0_zzz_yzzz, g_z_0_z_0_zzz_yzzzz, g_z_0_z_0_zzz_zzzz, g_z_0_z_0_zzz_zzzzz, g_z_0_z_0_zzzz_xxxx, g_z_0_z_0_zzzz_xxxy, g_z_0_z_0_zzzz_xxxz, g_z_0_z_0_zzzz_xxyy, g_z_0_z_0_zzzz_xxyz, g_z_0_z_0_zzzz_xxzz, g_z_0_z_0_zzzz_xyyy, g_z_0_z_0_zzzz_xyyz, g_z_0_z_0_zzzz_xyzz, g_z_0_z_0_zzzz_xzzz, g_z_0_z_0_zzzz_yyyy, g_z_0_z_0_zzzz_yyyz, g_z_0_z_0_zzzz_yyzz, g_z_0_z_0_zzzz_yzzz, g_z_0_z_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_zzzz_xxxx[k] = -g_0_0_z_0_zzz_xxxx[k] - g_z_0_z_0_zzz_xxxx[k] * ab_z + g_z_0_z_0_zzz_xxxxz[k];

                g_z_0_z_0_zzzz_xxxy[k] = -g_0_0_z_0_zzz_xxxy[k] - g_z_0_z_0_zzz_xxxy[k] * ab_z + g_z_0_z_0_zzz_xxxyz[k];

                g_z_0_z_0_zzzz_xxxz[k] = -g_0_0_z_0_zzz_xxxz[k] - g_z_0_z_0_zzz_xxxz[k] * ab_z + g_z_0_z_0_zzz_xxxzz[k];

                g_z_0_z_0_zzzz_xxyy[k] = -g_0_0_z_0_zzz_xxyy[k] - g_z_0_z_0_zzz_xxyy[k] * ab_z + g_z_0_z_0_zzz_xxyyz[k];

                g_z_0_z_0_zzzz_xxyz[k] = -g_0_0_z_0_zzz_xxyz[k] - g_z_0_z_0_zzz_xxyz[k] * ab_z + g_z_0_z_0_zzz_xxyzz[k];

                g_z_0_z_0_zzzz_xxzz[k] = -g_0_0_z_0_zzz_xxzz[k] - g_z_0_z_0_zzz_xxzz[k] * ab_z + g_z_0_z_0_zzz_xxzzz[k];

                g_z_0_z_0_zzzz_xyyy[k] = -g_0_0_z_0_zzz_xyyy[k] - g_z_0_z_0_zzz_xyyy[k] * ab_z + g_z_0_z_0_zzz_xyyyz[k];

                g_z_0_z_0_zzzz_xyyz[k] = -g_0_0_z_0_zzz_xyyz[k] - g_z_0_z_0_zzz_xyyz[k] * ab_z + g_z_0_z_0_zzz_xyyzz[k];

                g_z_0_z_0_zzzz_xyzz[k] = -g_0_0_z_0_zzz_xyzz[k] - g_z_0_z_0_zzz_xyzz[k] * ab_z + g_z_0_z_0_zzz_xyzzz[k];

                g_z_0_z_0_zzzz_xzzz[k] = -g_0_0_z_0_zzz_xzzz[k] - g_z_0_z_0_zzz_xzzz[k] * ab_z + g_z_0_z_0_zzz_xzzzz[k];

                g_z_0_z_0_zzzz_yyyy[k] = -g_0_0_z_0_zzz_yyyy[k] - g_z_0_z_0_zzz_yyyy[k] * ab_z + g_z_0_z_0_zzz_yyyyz[k];

                g_z_0_z_0_zzzz_yyyz[k] = -g_0_0_z_0_zzz_yyyz[k] - g_z_0_z_0_zzz_yyyz[k] * ab_z + g_z_0_z_0_zzz_yyyzz[k];

                g_z_0_z_0_zzzz_yyzz[k] = -g_0_0_z_0_zzz_yyzz[k] - g_z_0_z_0_zzz_yyzz[k] * ab_z + g_z_0_z_0_zzz_yyzzz[k];

                g_z_0_z_0_zzzz_yzzz[k] = -g_0_0_z_0_zzz_yzzz[k] - g_z_0_z_0_zzz_yzzz[k] * ab_z + g_z_0_z_0_zzz_yzzzz[k];

                g_z_0_z_0_zzzz_zzzz[k] = -g_0_0_z_0_zzz_zzzz[k] - g_z_0_z_0_zzz_zzzz[k] * ab_z + g_z_0_z_0_zzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

