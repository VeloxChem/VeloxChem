#include "ElectronRepulsionGeom1000ContrRecGGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_ggxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_ggxx,
                                            const size_t idx_fgxx,
                                            const size_t idx_geom_10_fgxx,
                                            const size_t idx_geom_10_fhxx,
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

            const auto fg_off = idx_fgxx + i * dcomps + j;

            auto g_xxx_xxxx = cbuffer.data(fg_off + 0 * ccomps * dcomps);

            auto g_xxx_xxxy = cbuffer.data(fg_off + 1 * ccomps * dcomps);

            auto g_xxx_xxxz = cbuffer.data(fg_off + 2 * ccomps * dcomps);

            auto g_xxx_xxyy = cbuffer.data(fg_off + 3 * ccomps * dcomps);

            auto g_xxx_xxyz = cbuffer.data(fg_off + 4 * ccomps * dcomps);

            auto g_xxx_xxzz = cbuffer.data(fg_off + 5 * ccomps * dcomps);

            auto g_xxx_xyyy = cbuffer.data(fg_off + 6 * ccomps * dcomps);

            auto g_xxx_xyyz = cbuffer.data(fg_off + 7 * ccomps * dcomps);

            auto g_xxx_xyzz = cbuffer.data(fg_off + 8 * ccomps * dcomps);

            auto g_xxx_xzzz = cbuffer.data(fg_off + 9 * ccomps * dcomps);

            auto g_xxx_yyyy = cbuffer.data(fg_off + 10 * ccomps * dcomps);

            auto g_xxx_yyyz = cbuffer.data(fg_off + 11 * ccomps * dcomps);

            auto g_xxx_yyzz = cbuffer.data(fg_off + 12 * ccomps * dcomps);

            auto g_xxx_yzzz = cbuffer.data(fg_off + 13 * ccomps * dcomps);

            auto g_xxx_zzzz = cbuffer.data(fg_off + 14 * ccomps * dcomps);

            auto g_xxy_xxxx = cbuffer.data(fg_off + 15 * ccomps * dcomps);

            auto g_xxy_xxxy = cbuffer.data(fg_off + 16 * ccomps * dcomps);

            auto g_xxy_xxxz = cbuffer.data(fg_off + 17 * ccomps * dcomps);

            auto g_xxy_xxyy = cbuffer.data(fg_off + 18 * ccomps * dcomps);

            auto g_xxy_xxyz = cbuffer.data(fg_off + 19 * ccomps * dcomps);

            auto g_xxy_xxzz = cbuffer.data(fg_off + 20 * ccomps * dcomps);

            auto g_xxy_xyyy = cbuffer.data(fg_off + 21 * ccomps * dcomps);

            auto g_xxy_xyyz = cbuffer.data(fg_off + 22 * ccomps * dcomps);

            auto g_xxy_xyzz = cbuffer.data(fg_off + 23 * ccomps * dcomps);

            auto g_xxy_xzzz = cbuffer.data(fg_off + 24 * ccomps * dcomps);

            auto g_xxy_yyyy = cbuffer.data(fg_off + 25 * ccomps * dcomps);

            auto g_xxy_yyyz = cbuffer.data(fg_off + 26 * ccomps * dcomps);

            auto g_xxy_yyzz = cbuffer.data(fg_off + 27 * ccomps * dcomps);

            auto g_xxy_yzzz = cbuffer.data(fg_off + 28 * ccomps * dcomps);

            auto g_xxy_zzzz = cbuffer.data(fg_off + 29 * ccomps * dcomps);

            auto g_xxz_xxxx = cbuffer.data(fg_off + 30 * ccomps * dcomps);

            auto g_xxz_xxxy = cbuffer.data(fg_off + 31 * ccomps * dcomps);

            auto g_xxz_xxxz = cbuffer.data(fg_off + 32 * ccomps * dcomps);

            auto g_xxz_xxyy = cbuffer.data(fg_off + 33 * ccomps * dcomps);

            auto g_xxz_xxyz = cbuffer.data(fg_off + 34 * ccomps * dcomps);

            auto g_xxz_xxzz = cbuffer.data(fg_off + 35 * ccomps * dcomps);

            auto g_xxz_xyyy = cbuffer.data(fg_off + 36 * ccomps * dcomps);

            auto g_xxz_xyyz = cbuffer.data(fg_off + 37 * ccomps * dcomps);

            auto g_xxz_xyzz = cbuffer.data(fg_off + 38 * ccomps * dcomps);

            auto g_xxz_xzzz = cbuffer.data(fg_off + 39 * ccomps * dcomps);

            auto g_xxz_yyyy = cbuffer.data(fg_off + 40 * ccomps * dcomps);

            auto g_xxz_yyyz = cbuffer.data(fg_off + 41 * ccomps * dcomps);

            auto g_xxz_yyzz = cbuffer.data(fg_off + 42 * ccomps * dcomps);

            auto g_xxz_yzzz = cbuffer.data(fg_off + 43 * ccomps * dcomps);

            auto g_xxz_zzzz = cbuffer.data(fg_off + 44 * ccomps * dcomps);

            auto g_xyy_xxxx = cbuffer.data(fg_off + 45 * ccomps * dcomps);

            auto g_xyy_xxxy = cbuffer.data(fg_off + 46 * ccomps * dcomps);

            auto g_xyy_xxxz = cbuffer.data(fg_off + 47 * ccomps * dcomps);

            auto g_xyy_xxyy = cbuffer.data(fg_off + 48 * ccomps * dcomps);

            auto g_xyy_xxyz = cbuffer.data(fg_off + 49 * ccomps * dcomps);

            auto g_xyy_xxzz = cbuffer.data(fg_off + 50 * ccomps * dcomps);

            auto g_xyy_xyyy = cbuffer.data(fg_off + 51 * ccomps * dcomps);

            auto g_xyy_xyyz = cbuffer.data(fg_off + 52 * ccomps * dcomps);

            auto g_xyy_xyzz = cbuffer.data(fg_off + 53 * ccomps * dcomps);

            auto g_xyy_xzzz = cbuffer.data(fg_off + 54 * ccomps * dcomps);

            auto g_xyy_yyyy = cbuffer.data(fg_off + 55 * ccomps * dcomps);

            auto g_xyy_yyyz = cbuffer.data(fg_off + 56 * ccomps * dcomps);

            auto g_xyy_yyzz = cbuffer.data(fg_off + 57 * ccomps * dcomps);

            auto g_xyy_yzzz = cbuffer.data(fg_off + 58 * ccomps * dcomps);

            auto g_xyy_zzzz = cbuffer.data(fg_off + 59 * ccomps * dcomps);

            auto g_xyz_xxxx = cbuffer.data(fg_off + 60 * ccomps * dcomps);

            auto g_xyz_xxxy = cbuffer.data(fg_off + 61 * ccomps * dcomps);

            auto g_xyz_xxxz = cbuffer.data(fg_off + 62 * ccomps * dcomps);

            auto g_xyz_xxyy = cbuffer.data(fg_off + 63 * ccomps * dcomps);

            auto g_xyz_xxyz = cbuffer.data(fg_off + 64 * ccomps * dcomps);

            auto g_xyz_xxzz = cbuffer.data(fg_off + 65 * ccomps * dcomps);

            auto g_xyz_xyyy = cbuffer.data(fg_off + 66 * ccomps * dcomps);

            auto g_xyz_xyyz = cbuffer.data(fg_off + 67 * ccomps * dcomps);

            auto g_xyz_xyzz = cbuffer.data(fg_off + 68 * ccomps * dcomps);

            auto g_xyz_xzzz = cbuffer.data(fg_off + 69 * ccomps * dcomps);

            auto g_xyz_yyyy = cbuffer.data(fg_off + 70 * ccomps * dcomps);

            auto g_xyz_yyyz = cbuffer.data(fg_off + 71 * ccomps * dcomps);

            auto g_xyz_yyzz = cbuffer.data(fg_off + 72 * ccomps * dcomps);

            auto g_xyz_yzzz = cbuffer.data(fg_off + 73 * ccomps * dcomps);

            auto g_xyz_zzzz = cbuffer.data(fg_off + 74 * ccomps * dcomps);

            auto g_xzz_xxxx = cbuffer.data(fg_off + 75 * ccomps * dcomps);

            auto g_xzz_xxxy = cbuffer.data(fg_off + 76 * ccomps * dcomps);

            auto g_xzz_xxxz = cbuffer.data(fg_off + 77 * ccomps * dcomps);

            auto g_xzz_xxyy = cbuffer.data(fg_off + 78 * ccomps * dcomps);

            auto g_xzz_xxyz = cbuffer.data(fg_off + 79 * ccomps * dcomps);

            auto g_xzz_xxzz = cbuffer.data(fg_off + 80 * ccomps * dcomps);

            auto g_xzz_xyyy = cbuffer.data(fg_off + 81 * ccomps * dcomps);

            auto g_xzz_xyyz = cbuffer.data(fg_off + 82 * ccomps * dcomps);

            auto g_xzz_xyzz = cbuffer.data(fg_off + 83 * ccomps * dcomps);

            auto g_xzz_xzzz = cbuffer.data(fg_off + 84 * ccomps * dcomps);

            auto g_xzz_yyyy = cbuffer.data(fg_off + 85 * ccomps * dcomps);

            auto g_xzz_yyyz = cbuffer.data(fg_off + 86 * ccomps * dcomps);

            auto g_xzz_yyzz = cbuffer.data(fg_off + 87 * ccomps * dcomps);

            auto g_xzz_yzzz = cbuffer.data(fg_off + 88 * ccomps * dcomps);

            auto g_xzz_zzzz = cbuffer.data(fg_off + 89 * ccomps * dcomps);

            auto g_yyy_xxxx = cbuffer.data(fg_off + 90 * ccomps * dcomps);

            auto g_yyy_xxxy = cbuffer.data(fg_off + 91 * ccomps * dcomps);

            auto g_yyy_xxxz = cbuffer.data(fg_off + 92 * ccomps * dcomps);

            auto g_yyy_xxyy = cbuffer.data(fg_off + 93 * ccomps * dcomps);

            auto g_yyy_xxyz = cbuffer.data(fg_off + 94 * ccomps * dcomps);

            auto g_yyy_xxzz = cbuffer.data(fg_off + 95 * ccomps * dcomps);

            auto g_yyy_xyyy = cbuffer.data(fg_off + 96 * ccomps * dcomps);

            auto g_yyy_xyyz = cbuffer.data(fg_off + 97 * ccomps * dcomps);

            auto g_yyy_xyzz = cbuffer.data(fg_off + 98 * ccomps * dcomps);

            auto g_yyy_xzzz = cbuffer.data(fg_off + 99 * ccomps * dcomps);

            auto g_yyy_yyyy = cbuffer.data(fg_off + 100 * ccomps * dcomps);

            auto g_yyy_yyyz = cbuffer.data(fg_off + 101 * ccomps * dcomps);

            auto g_yyy_yyzz = cbuffer.data(fg_off + 102 * ccomps * dcomps);

            auto g_yyy_yzzz = cbuffer.data(fg_off + 103 * ccomps * dcomps);

            auto g_yyy_zzzz = cbuffer.data(fg_off + 104 * ccomps * dcomps);

            auto g_yyz_xxxx = cbuffer.data(fg_off + 105 * ccomps * dcomps);

            auto g_yyz_xxxy = cbuffer.data(fg_off + 106 * ccomps * dcomps);

            auto g_yyz_xxxz = cbuffer.data(fg_off + 107 * ccomps * dcomps);

            auto g_yyz_xxyy = cbuffer.data(fg_off + 108 * ccomps * dcomps);

            auto g_yyz_xxyz = cbuffer.data(fg_off + 109 * ccomps * dcomps);

            auto g_yyz_xxzz = cbuffer.data(fg_off + 110 * ccomps * dcomps);

            auto g_yyz_xyyy = cbuffer.data(fg_off + 111 * ccomps * dcomps);

            auto g_yyz_xyyz = cbuffer.data(fg_off + 112 * ccomps * dcomps);

            auto g_yyz_xyzz = cbuffer.data(fg_off + 113 * ccomps * dcomps);

            auto g_yyz_xzzz = cbuffer.data(fg_off + 114 * ccomps * dcomps);

            auto g_yyz_yyyy = cbuffer.data(fg_off + 115 * ccomps * dcomps);

            auto g_yyz_yyyz = cbuffer.data(fg_off + 116 * ccomps * dcomps);

            auto g_yyz_yyzz = cbuffer.data(fg_off + 117 * ccomps * dcomps);

            auto g_yyz_yzzz = cbuffer.data(fg_off + 118 * ccomps * dcomps);

            auto g_yyz_zzzz = cbuffer.data(fg_off + 119 * ccomps * dcomps);

            auto g_yzz_xxxx = cbuffer.data(fg_off + 120 * ccomps * dcomps);

            auto g_yzz_xxxy = cbuffer.data(fg_off + 121 * ccomps * dcomps);

            auto g_yzz_xxxz = cbuffer.data(fg_off + 122 * ccomps * dcomps);

            auto g_yzz_xxyy = cbuffer.data(fg_off + 123 * ccomps * dcomps);

            auto g_yzz_xxyz = cbuffer.data(fg_off + 124 * ccomps * dcomps);

            auto g_yzz_xxzz = cbuffer.data(fg_off + 125 * ccomps * dcomps);

            auto g_yzz_xyyy = cbuffer.data(fg_off + 126 * ccomps * dcomps);

            auto g_yzz_xyyz = cbuffer.data(fg_off + 127 * ccomps * dcomps);

            auto g_yzz_xyzz = cbuffer.data(fg_off + 128 * ccomps * dcomps);

            auto g_yzz_xzzz = cbuffer.data(fg_off + 129 * ccomps * dcomps);

            auto g_yzz_yyyy = cbuffer.data(fg_off + 130 * ccomps * dcomps);

            auto g_yzz_yyyz = cbuffer.data(fg_off + 131 * ccomps * dcomps);

            auto g_yzz_yyzz = cbuffer.data(fg_off + 132 * ccomps * dcomps);

            auto g_yzz_yzzz = cbuffer.data(fg_off + 133 * ccomps * dcomps);

            auto g_yzz_zzzz = cbuffer.data(fg_off + 134 * ccomps * dcomps);

            auto g_zzz_xxxx = cbuffer.data(fg_off + 135 * ccomps * dcomps);

            auto g_zzz_xxxy = cbuffer.data(fg_off + 136 * ccomps * dcomps);

            auto g_zzz_xxxz = cbuffer.data(fg_off + 137 * ccomps * dcomps);

            auto g_zzz_xxyy = cbuffer.data(fg_off + 138 * ccomps * dcomps);

            auto g_zzz_xxyz = cbuffer.data(fg_off + 139 * ccomps * dcomps);

            auto g_zzz_xxzz = cbuffer.data(fg_off + 140 * ccomps * dcomps);

            auto g_zzz_xyyy = cbuffer.data(fg_off + 141 * ccomps * dcomps);

            auto g_zzz_xyyz = cbuffer.data(fg_off + 142 * ccomps * dcomps);

            auto g_zzz_xyzz = cbuffer.data(fg_off + 143 * ccomps * dcomps);

            auto g_zzz_xzzz = cbuffer.data(fg_off + 144 * ccomps * dcomps);

            auto g_zzz_yyyy = cbuffer.data(fg_off + 145 * ccomps * dcomps);

            auto g_zzz_yyyz = cbuffer.data(fg_off + 146 * ccomps * dcomps);

            auto g_zzz_yyzz = cbuffer.data(fg_off + 147 * ccomps * dcomps);

            auto g_zzz_yzzz = cbuffer.data(fg_off + 148 * ccomps * dcomps);

            auto g_zzz_zzzz = cbuffer.data(fg_off + 149 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FGSS

            const auto fg_geom_10_off = idx_geom_10_fgxx + i * dcomps + j;

            auto g_x_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 149 * ccomps * dcomps);

            auto g_y_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 164 * ccomps * dcomps);

            auto g_y_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 179 * ccomps * dcomps);

            auto g_y_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 180 * ccomps * dcomps);

            auto g_y_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 181 * ccomps * dcomps);

            auto g_y_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 182 * ccomps * dcomps);

            auto g_y_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 183 * ccomps * dcomps);

            auto g_y_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 184 * ccomps * dcomps);

            auto g_y_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 185 * ccomps * dcomps);

            auto g_y_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 186 * ccomps * dcomps);

            auto g_y_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 187 * ccomps * dcomps);

            auto g_y_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 188 * ccomps * dcomps);

            auto g_y_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 189 * ccomps * dcomps);

            auto g_y_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 190 * ccomps * dcomps);

            auto g_y_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 191 * ccomps * dcomps);

            auto g_y_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 192 * ccomps * dcomps);

            auto g_y_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 193 * ccomps * dcomps);

            auto g_y_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 194 * ccomps * dcomps);

            auto g_y_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 195 * ccomps * dcomps);

            auto g_y_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 196 * ccomps * dcomps);

            auto g_y_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 197 * ccomps * dcomps);

            auto g_y_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 198 * ccomps * dcomps);

            auto g_y_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 199 * ccomps * dcomps);

            auto g_y_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 200 * ccomps * dcomps);

            auto g_y_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 201 * ccomps * dcomps);

            auto g_y_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 202 * ccomps * dcomps);

            auto g_y_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 203 * ccomps * dcomps);

            auto g_y_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 204 * ccomps * dcomps);

            auto g_y_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 205 * ccomps * dcomps);

            auto g_y_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 206 * ccomps * dcomps);

            auto g_y_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 207 * ccomps * dcomps);

            auto g_y_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 208 * ccomps * dcomps);

            auto g_y_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 209 * ccomps * dcomps);

            auto g_y_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 219 * ccomps * dcomps);

            auto g_y_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 223 * ccomps * dcomps);

            auto g_y_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 251 * ccomps * dcomps);

            auto g_y_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 269 * ccomps * dcomps);

            auto g_y_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 299 * ccomps * dcomps);

            auto g_z_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 300 * ccomps * dcomps);

            auto g_z_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 301 * ccomps * dcomps);

            auto g_z_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 302 * ccomps * dcomps);

            auto g_z_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 303 * ccomps * dcomps);

            auto g_z_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 304 * ccomps * dcomps);

            auto g_z_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 305 * ccomps * dcomps);

            auto g_z_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 306 * ccomps * dcomps);

            auto g_z_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 307 * ccomps * dcomps);

            auto g_z_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 308 * ccomps * dcomps);

            auto g_z_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 309 * ccomps * dcomps);

            auto g_z_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 310 * ccomps * dcomps);

            auto g_z_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 311 * ccomps * dcomps);

            auto g_z_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 312 * ccomps * dcomps);

            auto g_z_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 313 * ccomps * dcomps);

            auto g_z_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 314 * ccomps * dcomps);

            auto g_z_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 315 * ccomps * dcomps);

            auto g_z_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 316 * ccomps * dcomps);

            auto g_z_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 317 * ccomps * dcomps);

            auto g_z_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 318 * ccomps * dcomps);

            auto g_z_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 319 * ccomps * dcomps);

            auto g_z_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 320 * ccomps * dcomps);

            auto g_z_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 321 * ccomps * dcomps);

            auto g_z_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 322 * ccomps * dcomps);

            auto g_z_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 323 * ccomps * dcomps);

            auto g_z_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 324 * ccomps * dcomps);

            auto g_z_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 325 * ccomps * dcomps);

            auto g_z_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 326 * ccomps * dcomps);

            auto g_z_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 327 * ccomps * dcomps);

            auto g_z_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 328 * ccomps * dcomps);

            auto g_z_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 329 * ccomps * dcomps);

            auto g_z_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 330 * ccomps * dcomps);

            auto g_z_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 331 * ccomps * dcomps);

            auto g_z_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 332 * ccomps * dcomps);

            auto g_z_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 333 * ccomps * dcomps);

            auto g_z_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 334 * ccomps * dcomps);

            auto g_z_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 335 * ccomps * dcomps);

            auto g_z_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 336 * ccomps * dcomps);

            auto g_z_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 337 * ccomps * dcomps);

            auto g_z_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 338 * ccomps * dcomps);

            auto g_z_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 339 * ccomps * dcomps);

            auto g_z_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 340 * ccomps * dcomps);

            auto g_z_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 341 * ccomps * dcomps);

            auto g_z_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 342 * ccomps * dcomps);

            auto g_z_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 343 * ccomps * dcomps);

            auto g_z_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 344 * ccomps * dcomps);

            auto g_z_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 345 * ccomps * dcomps);

            auto g_z_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 346 * ccomps * dcomps);

            auto g_z_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 347 * ccomps * dcomps);

            auto g_z_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 348 * ccomps * dcomps);

            auto g_z_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 349 * ccomps * dcomps);

            auto g_z_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 350 * ccomps * dcomps);

            auto g_z_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 351 * ccomps * dcomps);

            auto g_z_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 352 * ccomps * dcomps);

            auto g_z_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 353 * ccomps * dcomps);

            auto g_z_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 354 * ccomps * dcomps);

            auto g_z_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 355 * ccomps * dcomps);

            auto g_z_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 356 * ccomps * dcomps);

            auto g_z_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 357 * ccomps * dcomps);

            auto g_z_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 358 * ccomps * dcomps);

            auto g_z_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 359 * ccomps * dcomps);

            auto g_z_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 360 * ccomps * dcomps);

            auto g_z_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 361 * ccomps * dcomps);

            auto g_z_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 362 * ccomps * dcomps);

            auto g_z_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 363 * ccomps * dcomps);

            auto g_z_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 364 * ccomps * dcomps);

            auto g_z_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 365 * ccomps * dcomps);

            auto g_z_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 366 * ccomps * dcomps);

            auto g_z_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 367 * ccomps * dcomps);

            auto g_z_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 368 * ccomps * dcomps);

            auto g_z_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 369 * ccomps * dcomps);

            auto g_z_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 370 * ccomps * dcomps);

            auto g_z_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 371 * ccomps * dcomps);

            auto g_z_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 372 * ccomps * dcomps);

            auto g_z_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 373 * ccomps * dcomps);

            auto g_z_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 374 * ccomps * dcomps);

            auto g_z_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 375 * ccomps * dcomps);

            auto g_z_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 376 * ccomps * dcomps);

            auto g_z_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 377 * ccomps * dcomps);

            auto g_z_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 378 * ccomps * dcomps);

            auto g_z_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 379 * ccomps * dcomps);

            auto g_z_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 380 * ccomps * dcomps);

            auto g_z_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 381 * ccomps * dcomps);

            auto g_z_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 382 * ccomps * dcomps);

            auto g_z_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 383 * ccomps * dcomps);

            auto g_z_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 384 * ccomps * dcomps);

            auto g_z_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 385 * ccomps * dcomps);

            auto g_z_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 386 * ccomps * dcomps);

            auto g_z_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 387 * ccomps * dcomps);

            auto g_z_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 388 * ccomps * dcomps);

            auto g_z_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 389 * ccomps * dcomps);

            auto g_z_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 390 * ccomps * dcomps);

            auto g_z_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 391 * ccomps * dcomps);

            auto g_z_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 392 * ccomps * dcomps);

            auto g_z_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 393 * ccomps * dcomps);

            auto g_z_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 394 * ccomps * dcomps);

            auto g_z_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 395 * ccomps * dcomps);

            auto g_z_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 396 * ccomps * dcomps);

            auto g_z_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 397 * ccomps * dcomps);

            auto g_z_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 398 * ccomps * dcomps);

            auto g_z_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 399 * ccomps * dcomps);

            auto g_z_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 400 * ccomps * dcomps);

            auto g_z_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 401 * ccomps * dcomps);

            auto g_z_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 402 * ccomps * dcomps);

            auto g_z_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 403 * ccomps * dcomps);

            auto g_z_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 404 * ccomps * dcomps);

            auto g_z_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 405 * ccomps * dcomps);

            auto g_z_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 406 * ccomps * dcomps);

            auto g_z_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 407 * ccomps * dcomps);

            auto g_z_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 408 * ccomps * dcomps);

            auto g_z_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 409 * ccomps * dcomps);

            auto g_z_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 410 * ccomps * dcomps);

            auto g_z_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 411 * ccomps * dcomps);

            auto g_z_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 412 * ccomps * dcomps);

            auto g_z_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 413 * ccomps * dcomps);

            auto g_z_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 414 * ccomps * dcomps);

            auto g_z_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 415 * ccomps * dcomps);

            auto g_z_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 416 * ccomps * dcomps);

            auto g_z_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 417 * ccomps * dcomps);

            auto g_z_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 418 * ccomps * dcomps);

            auto g_z_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 419 * ccomps * dcomps);

            auto g_z_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 420 * ccomps * dcomps);

            auto g_z_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 421 * ccomps * dcomps);

            auto g_z_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 422 * ccomps * dcomps);

            auto g_z_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 423 * ccomps * dcomps);

            auto g_z_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 424 * ccomps * dcomps);

            auto g_z_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 425 * ccomps * dcomps);

            auto g_z_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 426 * ccomps * dcomps);

            auto g_z_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 427 * ccomps * dcomps);

            auto g_z_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 428 * ccomps * dcomps);

            auto g_z_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 429 * ccomps * dcomps);

            auto g_z_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 430 * ccomps * dcomps);

            auto g_z_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 431 * ccomps * dcomps);

            auto g_z_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 432 * ccomps * dcomps);

            auto g_z_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 433 * ccomps * dcomps);

            auto g_z_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 434 * ccomps * dcomps);

            auto g_z_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 435 * ccomps * dcomps);

            auto g_z_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 436 * ccomps * dcomps);

            auto g_z_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 437 * ccomps * dcomps);

            auto g_z_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 438 * ccomps * dcomps);

            auto g_z_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 439 * ccomps * dcomps);

            auto g_z_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 440 * ccomps * dcomps);

            auto g_z_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 441 * ccomps * dcomps);

            auto g_z_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 442 * ccomps * dcomps);

            auto g_z_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 443 * ccomps * dcomps);

            auto g_z_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 444 * ccomps * dcomps);

            auto g_z_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 445 * ccomps * dcomps);

            auto g_z_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 446 * ccomps * dcomps);

            auto g_z_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 447 * ccomps * dcomps);

            auto g_z_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 448 * ccomps * dcomps);

            auto g_z_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 449 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FHSS

            const auto fh_geom_10_off = idx_geom_10_fhxx + i * dcomps + j;

            auto g_x_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 209 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 219 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 223 * ccomps * dcomps);

            auto g_y_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 251 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 269 * ccomps * dcomps);

            auto g_y_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 299 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 419 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 420 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 421 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 422 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 423 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 424 * ccomps * dcomps);

            auto g_z_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 425 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 426 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 427 * ccomps * dcomps);

            auto g_z_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 428 * ccomps * dcomps);

            auto g_z_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 429 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 430 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 431 * ccomps * dcomps);

            auto g_z_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 432 * ccomps * dcomps);

            auto g_z_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 433 * ccomps * dcomps);

            auto g_z_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 434 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 435 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 436 * ccomps * dcomps);

            auto g_z_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 437 * ccomps * dcomps);

            auto g_z_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 438 * ccomps * dcomps);

            auto g_z_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 439 * ccomps * dcomps);

            auto g_z_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 440 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 441 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 442 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 443 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 444 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 445 * ccomps * dcomps);

            auto g_z_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 446 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 447 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 448 * ccomps * dcomps);

            auto g_z_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 449 * ccomps * dcomps);

            auto g_z_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 450 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 451 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 452 * ccomps * dcomps);

            auto g_z_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 453 * ccomps * dcomps);

            auto g_z_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 454 * ccomps * dcomps);

            auto g_z_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 455 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 456 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 457 * ccomps * dcomps);

            auto g_z_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 458 * ccomps * dcomps);

            auto g_z_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 459 * ccomps * dcomps);

            auto g_z_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 460 * ccomps * dcomps);

            auto g_z_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 461 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 462 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 463 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 464 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 465 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 466 * ccomps * dcomps);

            auto g_z_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 467 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 468 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 469 * ccomps * dcomps);

            auto g_z_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 470 * ccomps * dcomps);

            auto g_z_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 471 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 472 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 473 * ccomps * dcomps);

            auto g_z_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 474 * ccomps * dcomps);

            auto g_z_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 475 * ccomps * dcomps);

            auto g_z_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 476 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 477 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 478 * ccomps * dcomps);

            auto g_z_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 479 * ccomps * dcomps);

            auto g_z_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 480 * ccomps * dcomps);

            auto g_z_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 481 * ccomps * dcomps);

            auto g_z_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 482 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 483 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 484 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 485 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 486 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 487 * ccomps * dcomps);

            auto g_z_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 488 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 489 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 490 * ccomps * dcomps);

            auto g_z_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 491 * ccomps * dcomps);

            auto g_z_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 492 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 493 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 494 * ccomps * dcomps);

            auto g_z_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 495 * ccomps * dcomps);

            auto g_z_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 496 * ccomps * dcomps);

            auto g_z_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 497 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 498 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 499 * ccomps * dcomps);

            auto g_z_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 500 * ccomps * dcomps);

            auto g_z_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 501 * ccomps * dcomps);

            auto g_z_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 502 * ccomps * dcomps);

            auto g_z_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 503 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 504 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 505 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 506 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 507 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 508 * ccomps * dcomps);

            auto g_z_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 509 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 510 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 511 * ccomps * dcomps);

            auto g_z_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 512 * ccomps * dcomps);

            auto g_z_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 513 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 514 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 515 * ccomps * dcomps);

            auto g_z_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 516 * ccomps * dcomps);

            auto g_z_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 517 * ccomps * dcomps);

            auto g_z_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 518 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 519 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 520 * ccomps * dcomps);

            auto g_z_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 521 * ccomps * dcomps);

            auto g_z_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 522 * ccomps * dcomps);

            auto g_z_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 523 * ccomps * dcomps);

            auto g_z_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 524 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 525 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 526 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 527 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 528 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 529 * ccomps * dcomps);

            auto g_z_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 530 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 531 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 532 * ccomps * dcomps);

            auto g_z_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 533 * ccomps * dcomps);

            auto g_z_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 534 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 535 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 536 * ccomps * dcomps);

            auto g_z_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 537 * ccomps * dcomps);

            auto g_z_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 538 * ccomps * dcomps);

            auto g_z_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 539 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 540 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 541 * ccomps * dcomps);

            auto g_z_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 542 * ccomps * dcomps);

            auto g_z_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 543 * ccomps * dcomps);

            auto g_z_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 544 * ccomps * dcomps);

            auto g_z_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 545 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 546 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 547 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 548 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 549 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 550 * ccomps * dcomps);

            auto g_z_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 551 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 552 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 553 * ccomps * dcomps);

            auto g_z_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 554 * ccomps * dcomps);

            auto g_z_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 555 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 556 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 557 * ccomps * dcomps);

            auto g_z_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 558 * ccomps * dcomps);

            auto g_z_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 559 * ccomps * dcomps);

            auto g_z_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 566 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 569 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 575 * ccomps * dcomps);

            auto g_z_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 579 * ccomps * dcomps);

            auto g_z_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 584 * ccomps * dcomps);

            auto g_z_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 587 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 589 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 599 * ccomps * dcomps);

            auto g_z_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 608 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 609 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 611 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 614 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 615 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 619 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 629 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ggxx

            const auto gg_geom_10_off = idx_geom_10_ggxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_xxxx, g_x_0_xxx_xxxxx, g_x_0_xxx_xxxxy, g_x_0_xxx_xxxxz, g_x_0_xxx_xxxy, g_x_0_xxx_xxxyy, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxz, g_x_0_xxx_xxxzz, g_x_0_xxx_xxyy, g_x_0_xxx_xxyyy, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxzz, g_x_0_xxx_xxzzz, g_x_0_xxx_xyyy, g_x_0_xxx_xyyyy, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xzzz, g_x_0_xxx_xzzzz, g_x_0_xxx_yyyy, g_x_0_xxx_yyyz, g_x_0_xxx_yyzz, g_x_0_xxx_yzzz, g_x_0_xxx_zzzz, g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_zzzz, g_xxx_xxxx, g_xxx_xxxy, g_xxx_xxxz, g_xxx_xxyy, g_xxx_xxyz, g_xxx_xxzz, g_xxx_xyyy, g_xxx_xyyz, g_xxx_xyzz, g_xxx_xzzz, g_xxx_yyyy, g_xxx_yyyz, g_xxx_yyzz, g_xxx_yzzz, g_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxx_xxxx[k] = -g_xxx_xxxx[k] - g_x_0_xxx_xxxx[k] * ab_x + g_x_0_xxx_xxxxx[k];

                g_x_0_xxxx_xxxy[k] = -g_xxx_xxxy[k] - g_x_0_xxx_xxxy[k] * ab_x + g_x_0_xxx_xxxxy[k];

                g_x_0_xxxx_xxxz[k] = -g_xxx_xxxz[k] - g_x_0_xxx_xxxz[k] * ab_x + g_x_0_xxx_xxxxz[k];

                g_x_0_xxxx_xxyy[k] = -g_xxx_xxyy[k] - g_x_0_xxx_xxyy[k] * ab_x + g_x_0_xxx_xxxyy[k];

                g_x_0_xxxx_xxyz[k] = -g_xxx_xxyz[k] - g_x_0_xxx_xxyz[k] * ab_x + g_x_0_xxx_xxxyz[k];

                g_x_0_xxxx_xxzz[k] = -g_xxx_xxzz[k] - g_x_0_xxx_xxzz[k] * ab_x + g_x_0_xxx_xxxzz[k];

                g_x_0_xxxx_xyyy[k] = -g_xxx_xyyy[k] - g_x_0_xxx_xyyy[k] * ab_x + g_x_0_xxx_xxyyy[k];

                g_x_0_xxxx_xyyz[k] = -g_xxx_xyyz[k] - g_x_0_xxx_xyyz[k] * ab_x + g_x_0_xxx_xxyyz[k];

                g_x_0_xxxx_xyzz[k] = -g_xxx_xyzz[k] - g_x_0_xxx_xyzz[k] * ab_x + g_x_0_xxx_xxyzz[k];

                g_x_0_xxxx_xzzz[k] = -g_xxx_xzzz[k] - g_x_0_xxx_xzzz[k] * ab_x + g_x_0_xxx_xxzzz[k];

                g_x_0_xxxx_yyyy[k] = -g_xxx_yyyy[k] - g_x_0_xxx_yyyy[k] * ab_x + g_x_0_xxx_xyyyy[k];

                g_x_0_xxxx_yyyz[k] = -g_xxx_yyyz[k] - g_x_0_xxx_yyyz[k] * ab_x + g_x_0_xxx_xyyyz[k];

                g_x_0_xxxx_yyzz[k] = -g_xxx_yyzz[k] - g_x_0_xxx_yyzz[k] * ab_x + g_x_0_xxx_xyyzz[k];

                g_x_0_xxxx_yzzz[k] = -g_xxx_yzzz[k] - g_x_0_xxx_yzzz[k] * ab_x + g_x_0_xxx_xyzzz[k];

                g_x_0_xxxx_zzzz[k] = -g_xxx_zzzz[k] - g_x_0_xxx_zzzz[k] * ab_x + g_x_0_xxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_xxxx, g_x_0_xxx_xxxxy, g_x_0_xxx_xxxy, g_x_0_xxx_xxxyy, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxz, g_x_0_xxx_xxyy, g_x_0_xxx_xxyyy, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxzz, g_x_0_xxx_xyyy, g_x_0_xxx_xyyyy, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xzzz, g_x_0_xxx_yyyy, g_x_0_xxx_yyyyy, g_x_0_xxx_yyyyz, g_x_0_xxx_yyyz, g_x_0_xxx_yyyzz, g_x_0_xxx_yyzz, g_x_0_xxx_yyzzz, g_x_0_xxx_yzzz, g_x_0_xxx_yzzzz, g_x_0_xxx_zzzz, g_x_0_xxxy_xxxx, g_x_0_xxxy_xxxy, g_x_0_xxxy_xxxz, g_x_0_xxxy_xxyy, g_x_0_xxxy_xxyz, g_x_0_xxxy_xxzz, g_x_0_xxxy_xyyy, g_x_0_xxxy_xyyz, g_x_0_xxxy_xyzz, g_x_0_xxxy_xzzz, g_x_0_xxxy_yyyy, g_x_0_xxxy_yyyz, g_x_0_xxxy_yyzz, g_x_0_xxxy_yzzz, g_x_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxy_xxxx[k] = -g_x_0_xxx_xxxx[k] * ab_y + g_x_0_xxx_xxxxy[k];

                g_x_0_xxxy_xxxy[k] = -g_x_0_xxx_xxxy[k] * ab_y + g_x_0_xxx_xxxyy[k];

                g_x_0_xxxy_xxxz[k] = -g_x_0_xxx_xxxz[k] * ab_y + g_x_0_xxx_xxxyz[k];

                g_x_0_xxxy_xxyy[k] = -g_x_0_xxx_xxyy[k] * ab_y + g_x_0_xxx_xxyyy[k];

                g_x_0_xxxy_xxyz[k] = -g_x_0_xxx_xxyz[k] * ab_y + g_x_0_xxx_xxyyz[k];

                g_x_0_xxxy_xxzz[k] = -g_x_0_xxx_xxzz[k] * ab_y + g_x_0_xxx_xxyzz[k];

                g_x_0_xxxy_xyyy[k] = -g_x_0_xxx_xyyy[k] * ab_y + g_x_0_xxx_xyyyy[k];

                g_x_0_xxxy_xyyz[k] = -g_x_0_xxx_xyyz[k] * ab_y + g_x_0_xxx_xyyyz[k];

                g_x_0_xxxy_xyzz[k] = -g_x_0_xxx_xyzz[k] * ab_y + g_x_0_xxx_xyyzz[k];

                g_x_0_xxxy_xzzz[k] = -g_x_0_xxx_xzzz[k] * ab_y + g_x_0_xxx_xyzzz[k];

                g_x_0_xxxy_yyyy[k] = -g_x_0_xxx_yyyy[k] * ab_y + g_x_0_xxx_yyyyy[k];

                g_x_0_xxxy_yyyz[k] = -g_x_0_xxx_yyyz[k] * ab_y + g_x_0_xxx_yyyyz[k];

                g_x_0_xxxy_yyzz[k] = -g_x_0_xxx_yyzz[k] * ab_y + g_x_0_xxx_yyyzz[k];

                g_x_0_xxxy_yzzz[k] = -g_x_0_xxx_yzzz[k] * ab_y + g_x_0_xxx_yyzzz[k];

                g_x_0_xxxy_zzzz[k] = -g_x_0_xxx_zzzz[k] * ab_y + g_x_0_xxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_xxxx, g_x_0_xxx_xxxxz, g_x_0_xxx_xxxy, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxz, g_x_0_xxx_xxxzz, g_x_0_xxx_xxyy, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxzz, g_x_0_xxx_xxzzz, g_x_0_xxx_xyyy, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xzzz, g_x_0_xxx_xzzzz, g_x_0_xxx_yyyy, g_x_0_xxx_yyyyz, g_x_0_xxx_yyyz, g_x_0_xxx_yyyzz, g_x_0_xxx_yyzz, g_x_0_xxx_yyzzz, g_x_0_xxx_yzzz, g_x_0_xxx_yzzzz, g_x_0_xxx_zzzz, g_x_0_xxx_zzzzz, g_x_0_xxxz_xxxx, g_x_0_xxxz_xxxy, g_x_0_xxxz_xxxz, g_x_0_xxxz_xxyy, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxzz, g_x_0_xxxz_xyyy, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xzzz, g_x_0_xxxz_yyyy, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxz_xxxx[k] = -g_x_0_xxx_xxxx[k] * ab_z + g_x_0_xxx_xxxxz[k];

                g_x_0_xxxz_xxxy[k] = -g_x_0_xxx_xxxy[k] * ab_z + g_x_0_xxx_xxxyz[k];

                g_x_0_xxxz_xxxz[k] = -g_x_0_xxx_xxxz[k] * ab_z + g_x_0_xxx_xxxzz[k];

                g_x_0_xxxz_xxyy[k] = -g_x_0_xxx_xxyy[k] * ab_z + g_x_0_xxx_xxyyz[k];

                g_x_0_xxxz_xxyz[k] = -g_x_0_xxx_xxyz[k] * ab_z + g_x_0_xxx_xxyzz[k];

                g_x_0_xxxz_xxzz[k] = -g_x_0_xxx_xxzz[k] * ab_z + g_x_0_xxx_xxzzz[k];

                g_x_0_xxxz_xyyy[k] = -g_x_0_xxx_xyyy[k] * ab_z + g_x_0_xxx_xyyyz[k];

                g_x_0_xxxz_xyyz[k] = -g_x_0_xxx_xyyz[k] * ab_z + g_x_0_xxx_xyyzz[k];

                g_x_0_xxxz_xyzz[k] = -g_x_0_xxx_xyzz[k] * ab_z + g_x_0_xxx_xyzzz[k];

                g_x_0_xxxz_xzzz[k] = -g_x_0_xxx_xzzz[k] * ab_z + g_x_0_xxx_xzzzz[k];

                g_x_0_xxxz_yyyy[k] = -g_x_0_xxx_yyyy[k] * ab_z + g_x_0_xxx_yyyyz[k];

                g_x_0_xxxz_yyyz[k] = -g_x_0_xxx_yyyz[k] * ab_z + g_x_0_xxx_yyyzz[k];

                g_x_0_xxxz_yyzz[k] = -g_x_0_xxx_yyzz[k] * ab_z + g_x_0_xxx_yyzzz[k];

                g_x_0_xxxz_yzzz[k] = -g_x_0_xxx_yzzz[k] * ab_z + g_x_0_xxx_yzzzz[k];

                g_x_0_xxxz_zzzz[k] = -g_x_0_xxx_zzzz[k] * ab_z + g_x_0_xxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxy_xxxx, g_x_0_xxy_xxxxy, g_x_0_xxy_xxxy, g_x_0_xxy_xxxyy, g_x_0_xxy_xxxyz, g_x_0_xxy_xxxz, g_x_0_xxy_xxyy, g_x_0_xxy_xxyyy, g_x_0_xxy_xxyyz, g_x_0_xxy_xxyz, g_x_0_xxy_xxyzz, g_x_0_xxy_xxzz, g_x_0_xxy_xyyy, g_x_0_xxy_xyyyy, g_x_0_xxy_xyyyz, g_x_0_xxy_xyyz, g_x_0_xxy_xyyzz, g_x_0_xxy_xyzz, g_x_0_xxy_xyzzz, g_x_0_xxy_xzzz, g_x_0_xxy_yyyy, g_x_0_xxy_yyyyy, g_x_0_xxy_yyyyz, g_x_0_xxy_yyyz, g_x_0_xxy_yyyzz, g_x_0_xxy_yyzz, g_x_0_xxy_yyzzz, g_x_0_xxy_yzzz, g_x_0_xxy_yzzzz, g_x_0_xxy_zzzz, g_x_0_xxyy_xxxx, g_x_0_xxyy_xxxy, g_x_0_xxyy_xxxz, g_x_0_xxyy_xxyy, g_x_0_xxyy_xxyz, g_x_0_xxyy_xxzz, g_x_0_xxyy_xyyy, g_x_0_xxyy_xyyz, g_x_0_xxyy_xyzz, g_x_0_xxyy_xzzz, g_x_0_xxyy_yyyy, g_x_0_xxyy_yyyz, g_x_0_xxyy_yyzz, g_x_0_xxyy_yzzz, g_x_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyy_xxxx[k] = -g_x_0_xxy_xxxx[k] * ab_y + g_x_0_xxy_xxxxy[k];

                g_x_0_xxyy_xxxy[k] = -g_x_0_xxy_xxxy[k] * ab_y + g_x_0_xxy_xxxyy[k];

                g_x_0_xxyy_xxxz[k] = -g_x_0_xxy_xxxz[k] * ab_y + g_x_0_xxy_xxxyz[k];

                g_x_0_xxyy_xxyy[k] = -g_x_0_xxy_xxyy[k] * ab_y + g_x_0_xxy_xxyyy[k];

                g_x_0_xxyy_xxyz[k] = -g_x_0_xxy_xxyz[k] * ab_y + g_x_0_xxy_xxyyz[k];

                g_x_0_xxyy_xxzz[k] = -g_x_0_xxy_xxzz[k] * ab_y + g_x_0_xxy_xxyzz[k];

                g_x_0_xxyy_xyyy[k] = -g_x_0_xxy_xyyy[k] * ab_y + g_x_0_xxy_xyyyy[k];

                g_x_0_xxyy_xyyz[k] = -g_x_0_xxy_xyyz[k] * ab_y + g_x_0_xxy_xyyyz[k];

                g_x_0_xxyy_xyzz[k] = -g_x_0_xxy_xyzz[k] * ab_y + g_x_0_xxy_xyyzz[k];

                g_x_0_xxyy_xzzz[k] = -g_x_0_xxy_xzzz[k] * ab_y + g_x_0_xxy_xyzzz[k];

                g_x_0_xxyy_yyyy[k] = -g_x_0_xxy_yyyy[k] * ab_y + g_x_0_xxy_yyyyy[k];

                g_x_0_xxyy_yyyz[k] = -g_x_0_xxy_yyyz[k] * ab_y + g_x_0_xxy_yyyyz[k];

                g_x_0_xxyy_yyzz[k] = -g_x_0_xxy_yyzz[k] * ab_y + g_x_0_xxy_yyyzz[k];

                g_x_0_xxyy_yzzz[k] = -g_x_0_xxy_yzzz[k] * ab_y + g_x_0_xxy_yyzzz[k];

                g_x_0_xxyy_zzzz[k] = -g_x_0_xxy_zzzz[k] * ab_y + g_x_0_xxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyz_xxxx, g_x_0_xxyz_xxxy, g_x_0_xxyz_xxxz, g_x_0_xxyz_xxyy, g_x_0_xxyz_xxyz, g_x_0_xxyz_xxzz, g_x_0_xxyz_xyyy, g_x_0_xxyz_xyyz, g_x_0_xxyz_xyzz, g_x_0_xxyz_xzzz, g_x_0_xxyz_yyyy, g_x_0_xxyz_yyyz, g_x_0_xxyz_yyzz, g_x_0_xxyz_yzzz, g_x_0_xxyz_zzzz, g_x_0_xxz_xxxx, g_x_0_xxz_xxxxy, g_x_0_xxz_xxxy, g_x_0_xxz_xxxyy, g_x_0_xxz_xxxyz, g_x_0_xxz_xxxz, g_x_0_xxz_xxyy, g_x_0_xxz_xxyyy, g_x_0_xxz_xxyyz, g_x_0_xxz_xxyz, g_x_0_xxz_xxyzz, g_x_0_xxz_xxzz, g_x_0_xxz_xyyy, g_x_0_xxz_xyyyy, g_x_0_xxz_xyyyz, g_x_0_xxz_xyyz, g_x_0_xxz_xyyzz, g_x_0_xxz_xyzz, g_x_0_xxz_xyzzz, g_x_0_xxz_xzzz, g_x_0_xxz_yyyy, g_x_0_xxz_yyyyy, g_x_0_xxz_yyyyz, g_x_0_xxz_yyyz, g_x_0_xxz_yyyzz, g_x_0_xxz_yyzz, g_x_0_xxz_yyzzz, g_x_0_xxz_yzzz, g_x_0_xxz_yzzzz, g_x_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyz_xxxx[k] = -g_x_0_xxz_xxxx[k] * ab_y + g_x_0_xxz_xxxxy[k];

                g_x_0_xxyz_xxxy[k] = -g_x_0_xxz_xxxy[k] * ab_y + g_x_0_xxz_xxxyy[k];

                g_x_0_xxyz_xxxz[k] = -g_x_0_xxz_xxxz[k] * ab_y + g_x_0_xxz_xxxyz[k];

                g_x_0_xxyz_xxyy[k] = -g_x_0_xxz_xxyy[k] * ab_y + g_x_0_xxz_xxyyy[k];

                g_x_0_xxyz_xxyz[k] = -g_x_0_xxz_xxyz[k] * ab_y + g_x_0_xxz_xxyyz[k];

                g_x_0_xxyz_xxzz[k] = -g_x_0_xxz_xxzz[k] * ab_y + g_x_0_xxz_xxyzz[k];

                g_x_0_xxyz_xyyy[k] = -g_x_0_xxz_xyyy[k] * ab_y + g_x_0_xxz_xyyyy[k];

                g_x_0_xxyz_xyyz[k] = -g_x_0_xxz_xyyz[k] * ab_y + g_x_0_xxz_xyyyz[k];

                g_x_0_xxyz_xyzz[k] = -g_x_0_xxz_xyzz[k] * ab_y + g_x_0_xxz_xyyzz[k];

                g_x_0_xxyz_xzzz[k] = -g_x_0_xxz_xzzz[k] * ab_y + g_x_0_xxz_xyzzz[k];

                g_x_0_xxyz_yyyy[k] = -g_x_0_xxz_yyyy[k] * ab_y + g_x_0_xxz_yyyyy[k];

                g_x_0_xxyz_yyyz[k] = -g_x_0_xxz_yyyz[k] * ab_y + g_x_0_xxz_yyyyz[k];

                g_x_0_xxyz_yyzz[k] = -g_x_0_xxz_yyzz[k] * ab_y + g_x_0_xxz_yyyzz[k];

                g_x_0_xxyz_yzzz[k] = -g_x_0_xxz_yzzz[k] * ab_y + g_x_0_xxz_yyzzz[k];

                g_x_0_xxyz_zzzz[k] = -g_x_0_xxz_zzzz[k] * ab_y + g_x_0_xxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxz_xxxx, g_x_0_xxz_xxxxz, g_x_0_xxz_xxxy, g_x_0_xxz_xxxyz, g_x_0_xxz_xxxz, g_x_0_xxz_xxxzz, g_x_0_xxz_xxyy, g_x_0_xxz_xxyyz, g_x_0_xxz_xxyz, g_x_0_xxz_xxyzz, g_x_0_xxz_xxzz, g_x_0_xxz_xxzzz, g_x_0_xxz_xyyy, g_x_0_xxz_xyyyz, g_x_0_xxz_xyyz, g_x_0_xxz_xyyzz, g_x_0_xxz_xyzz, g_x_0_xxz_xyzzz, g_x_0_xxz_xzzz, g_x_0_xxz_xzzzz, g_x_0_xxz_yyyy, g_x_0_xxz_yyyyz, g_x_0_xxz_yyyz, g_x_0_xxz_yyyzz, g_x_0_xxz_yyzz, g_x_0_xxz_yyzzz, g_x_0_xxz_yzzz, g_x_0_xxz_yzzzz, g_x_0_xxz_zzzz, g_x_0_xxz_zzzzz, g_x_0_xxzz_xxxx, g_x_0_xxzz_xxxy, g_x_0_xxzz_xxxz, g_x_0_xxzz_xxyy, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxzz, g_x_0_xxzz_xyyy, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xzzz, g_x_0_xxzz_yyyy, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzz_xxxx[k] = -g_x_0_xxz_xxxx[k] * ab_z + g_x_0_xxz_xxxxz[k];

                g_x_0_xxzz_xxxy[k] = -g_x_0_xxz_xxxy[k] * ab_z + g_x_0_xxz_xxxyz[k];

                g_x_0_xxzz_xxxz[k] = -g_x_0_xxz_xxxz[k] * ab_z + g_x_0_xxz_xxxzz[k];

                g_x_0_xxzz_xxyy[k] = -g_x_0_xxz_xxyy[k] * ab_z + g_x_0_xxz_xxyyz[k];

                g_x_0_xxzz_xxyz[k] = -g_x_0_xxz_xxyz[k] * ab_z + g_x_0_xxz_xxyzz[k];

                g_x_0_xxzz_xxzz[k] = -g_x_0_xxz_xxzz[k] * ab_z + g_x_0_xxz_xxzzz[k];

                g_x_0_xxzz_xyyy[k] = -g_x_0_xxz_xyyy[k] * ab_z + g_x_0_xxz_xyyyz[k];

                g_x_0_xxzz_xyyz[k] = -g_x_0_xxz_xyyz[k] * ab_z + g_x_0_xxz_xyyzz[k];

                g_x_0_xxzz_xyzz[k] = -g_x_0_xxz_xyzz[k] * ab_z + g_x_0_xxz_xyzzz[k];

                g_x_0_xxzz_xzzz[k] = -g_x_0_xxz_xzzz[k] * ab_z + g_x_0_xxz_xzzzz[k];

                g_x_0_xxzz_yyyy[k] = -g_x_0_xxz_yyyy[k] * ab_z + g_x_0_xxz_yyyyz[k];

                g_x_0_xxzz_yyyz[k] = -g_x_0_xxz_yyyz[k] * ab_z + g_x_0_xxz_yyyzz[k];

                g_x_0_xxzz_yyzz[k] = -g_x_0_xxz_yyzz[k] * ab_z + g_x_0_xxz_yyzzz[k];

                g_x_0_xxzz_yzzz[k] = -g_x_0_xxz_yzzz[k] * ab_z + g_x_0_xxz_yzzzz[k];

                g_x_0_xxzz_zzzz[k] = -g_x_0_xxz_zzzz[k] * ab_z + g_x_0_xxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyy_xxxx, g_x_0_xyy_xxxxy, g_x_0_xyy_xxxy, g_x_0_xyy_xxxyy, g_x_0_xyy_xxxyz, g_x_0_xyy_xxxz, g_x_0_xyy_xxyy, g_x_0_xyy_xxyyy, g_x_0_xyy_xxyyz, g_x_0_xyy_xxyz, g_x_0_xyy_xxyzz, g_x_0_xyy_xxzz, g_x_0_xyy_xyyy, g_x_0_xyy_xyyyy, g_x_0_xyy_xyyyz, g_x_0_xyy_xyyz, g_x_0_xyy_xyyzz, g_x_0_xyy_xyzz, g_x_0_xyy_xyzzz, g_x_0_xyy_xzzz, g_x_0_xyy_yyyy, g_x_0_xyy_yyyyy, g_x_0_xyy_yyyyz, g_x_0_xyy_yyyz, g_x_0_xyy_yyyzz, g_x_0_xyy_yyzz, g_x_0_xyy_yyzzz, g_x_0_xyy_yzzz, g_x_0_xyy_yzzzz, g_x_0_xyy_zzzz, g_x_0_xyyy_xxxx, g_x_0_xyyy_xxxy, g_x_0_xyyy_xxxz, g_x_0_xyyy_xxyy, g_x_0_xyyy_xxyz, g_x_0_xyyy_xxzz, g_x_0_xyyy_xyyy, g_x_0_xyyy_xyyz, g_x_0_xyyy_xyzz, g_x_0_xyyy_xzzz, g_x_0_xyyy_yyyy, g_x_0_xyyy_yyyz, g_x_0_xyyy_yyzz, g_x_0_xyyy_yzzz, g_x_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyy_xxxx[k] = -g_x_0_xyy_xxxx[k] * ab_y + g_x_0_xyy_xxxxy[k];

                g_x_0_xyyy_xxxy[k] = -g_x_0_xyy_xxxy[k] * ab_y + g_x_0_xyy_xxxyy[k];

                g_x_0_xyyy_xxxz[k] = -g_x_0_xyy_xxxz[k] * ab_y + g_x_0_xyy_xxxyz[k];

                g_x_0_xyyy_xxyy[k] = -g_x_0_xyy_xxyy[k] * ab_y + g_x_0_xyy_xxyyy[k];

                g_x_0_xyyy_xxyz[k] = -g_x_0_xyy_xxyz[k] * ab_y + g_x_0_xyy_xxyyz[k];

                g_x_0_xyyy_xxzz[k] = -g_x_0_xyy_xxzz[k] * ab_y + g_x_0_xyy_xxyzz[k];

                g_x_0_xyyy_xyyy[k] = -g_x_0_xyy_xyyy[k] * ab_y + g_x_0_xyy_xyyyy[k];

                g_x_0_xyyy_xyyz[k] = -g_x_0_xyy_xyyz[k] * ab_y + g_x_0_xyy_xyyyz[k];

                g_x_0_xyyy_xyzz[k] = -g_x_0_xyy_xyzz[k] * ab_y + g_x_0_xyy_xyyzz[k];

                g_x_0_xyyy_xzzz[k] = -g_x_0_xyy_xzzz[k] * ab_y + g_x_0_xyy_xyzzz[k];

                g_x_0_xyyy_yyyy[k] = -g_x_0_xyy_yyyy[k] * ab_y + g_x_0_xyy_yyyyy[k];

                g_x_0_xyyy_yyyz[k] = -g_x_0_xyy_yyyz[k] * ab_y + g_x_0_xyy_yyyyz[k];

                g_x_0_xyyy_yyzz[k] = -g_x_0_xyy_yyzz[k] * ab_y + g_x_0_xyy_yyyzz[k];

                g_x_0_xyyy_yzzz[k] = -g_x_0_xyy_yzzz[k] * ab_y + g_x_0_xyy_yyzzz[k];

                g_x_0_xyyy_zzzz[k] = -g_x_0_xyy_zzzz[k] * ab_y + g_x_0_xyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyz_xxxx, g_x_0_xyyz_xxxy, g_x_0_xyyz_xxxz, g_x_0_xyyz_xxyy, g_x_0_xyyz_xxyz, g_x_0_xyyz_xxzz, g_x_0_xyyz_xyyy, g_x_0_xyyz_xyyz, g_x_0_xyyz_xyzz, g_x_0_xyyz_xzzz, g_x_0_xyyz_yyyy, g_x_0_xyyz_yyyz, g_x_0_xyyz_yyzz, g_x_0_xyyz_yzzz, g_x_0_xyyz_zzzz, g_x_0_xyz_xxxx, g_x_0_xyz_xxxxy, g_x_0_xyz_xxxy, g_x_0_xyz_xxxyy, g_x_0_xyz_xxxyz, g_x_0_xyz_xxxz, g_x_0_xyz_xxyy, g_x_0_xyz_xxyyy, g_x_0_xyz_xxyyz, g_x_0_xyz_xxyz, g_x_0_xyz_xxyzz, g_x_0_xyz_xxzz, g_x_0_xyz_xyyy, g_x_0_xyz_xyyyy, g_x_0_xyz_xyyyz, g_x_0_xyz_xyyz, g_x_0_xyz_xyyzz, g_x_0_xyz_xyzz, g_x_0_xyz_xyzzz, g_x_0_xyz_xzzz, g_x_0_xyz_yyyy, g_x_0_xyz_yyyyy, g_x_0_xyz_yyyyz, g_x_0_xyz_yyyz, g_x_0_xyz_yyyzz, g_x_0_xyz_yyzz, g_x_0_xyz_yyzzz, g_x_0_xyz_yzzz, g_x_0_xyz_yzzzz, g_x_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyz_xxxx[k] = -g_x_0_xyz_xxxx[k] * ab_y + g_x_0_xyz_xxxxy[k];

                g_x_0_xyyz_xxxy[k] = -g_x_0_xyz_xxxy[k] * ab_y + g_x_0_xyz_xxxyy[k];

                g_x_0_xyyz_xxxz[k] = -g_x_0_xyz_xxxz[k] * ab_y + g_x_0_xyz_xxxyz[k];

                g_x_0_xyyz_xxyy[k] = -g_x_0_xyz_xxyy[k] * ab_y + g_x_0_xyz_xxyyy[k];

                g_x_0_xyyz_xxyz[k] = -g_x_0_xyz_xxyz[k] * ab_y + g_x_0_xyz_xxyyz[k];

                g_x_0_xyyz_xxzz[k] = -g_x_0_xyz_xxzz[k] * ab_y + g_x_0_xyz_xxyzz[k];

                g_x_0_xyyz_xyyy[k] = -g_x_0_xyz_xyyy[k] * ab_y + g_x_0_xyz_xyyyy[k];

                g_x_0_xyyz_xyyz[k] = -g_x_0_xyz_xyyz[k] * ab_y + g_x_0_xyz_xyyyz[k];

                g_x_0_xyyz_xyzz[k] = -g_x_0_xyz_xyzz[k] * ab_y + g_x_0_xyz_xyyzz[k];

                g_x_0_xyyz_xzzz[k] = -g_x_0_xyz_xzzz[k] * ab_y + g_x_0_xyz_xyzzz[k];

                g_x_0_xyyz_yyyy[k] = -g_x_0_xyz_yyyy[k] * ab_y + g_x_0_xyz_yyyyy[k];

                g_x_0_xyyz_yyyz[k] = -g_x_0_xyz_yyyz[k] * ab_y + g_x_0_xyz_yyyyz[k];

                g_x_0_xyyz_yyzz[k] = -g_x_0_xyz_yyzz[k] * ab_y + g_x_0_xyz_yyyzz[k];

                g_x_0_xyyz_yzzz[k] = -g_x_0_xyz_yzzz[k] * ab_y + g_x_0_xyz_yyzzz[k];

                g_x_0_xyyz_zzzz[k] = -g_x_0_xyz_zzzz[k] * ab_y + g_x_0_xyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzz_xxxx, g_x_0_xyzz_xxxy, g_x_0_xyzz_xxxz, g_x_0_xyzz_xxyy, g_x_0_xyzz_xxyz, g_x_0_xyzz_xxzz, g_x_0_xyzz_xyyy, g_x_0_xyzz_xyyz, g_x_0_xyzz_xyzz, g_x_0_xyzz_xzzz, g_x_0_xyzz_yyyy, g_x_0_xyzz_yyyz, g_x_0_xyzz_yyzz, g_x_0_xyzz_yzzz, g_x_0_xyzz_zzzz, g_x_0_xzz_xxxx, g_x_0_xzz_xxxxy, g_x_0_xzz_xxxy, g_x_0_xzz_xxxyy, g_x_0_xzz_xxxyz, g_x_0_xzz_xxxz, g_x_0_xzz_xxyy, g_x_0_xzz_xxyyy, g_x_0_xzz_xxyyz, g_x_0_xzz_xxyz, g_x_0_xzz_xxyzz, g_x_0_xzz_xxzz, g_x_0_xzz_xyyy, g_x_0_xzz_xyyyy, g_x_0_xzz_xyyyz, g_x_0_xzz_xyyz, g_x_0_xzz_xyyzz, g_x_0_xzz_xyzz, g_x_0_xzz_xyzzz, g_x_0_xzz_xzzz, g_x_0_xzz_yyyy, g_x_0_xzz_yyyyy, g_x_0_xzz_yyyyz, g_x_0_xzz_yyyz, g_x_0_xzz_yyyzz, g_x_0_xzz_yyzz, g_x_0_xzz_yyzzz, g_x_0_xzz_yzzz, g_x_0_xzz_yzzzz, g_x_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzz_xxxx[k] = -g_x_0_xzz_xxxx[k] * ab_y + g_x_0_xzz_xxxxy[k];

                g_x_0_xyzz_xxxy[k] = -g_x_0_xzz_xxxy[k] * ab_y + g_x_0_xzz_xxxyy[k];

                g_x_0_xyzz_xxxz[k] = -g_x_0_xzz_xxxz[k] * ab_y + g_x_0_xzz_xxxyz[k];

                g_x_0_xyzz_xxyy[k] = -g_x_0_xzz_xxyy[k] * ab_y + g_x_0_xzz_xxyyy[k];

                g_x_0_xyzz_xxyz[k] = -g_x_0_xzz_xxyz[k] * ab_y + g_x_0_xzz_xxyyz[k];

                g_x_0_xyzz_xxzz[k] = -g_x_0_xzz_xxzz[k] * ab_y + g_x_0_xzz_xxyzz[k];

                g_x_0_xyzz_xyyy[k] = -g_x_0_xzz_xyyy[k] * ab_y + g_x_0_xzz_xyyyy[k];

                g_x_0_xyzz_xyyz[k] = -g_x_0_xzz_xyyz[k] * ab_y + g_x_0_xzz_xyyyz[k];

                g_x_0_xyzz_xyzz[k] = -g_x_0_xzz_xyzz[k] * ab_y + g_x_0_xzz_xyyzz[k];

                g_x_0_xyzz_xzzz[k] = -g_x_0_xzz_xzzz[k] * ab_y + g_x_0_xzz_xyzzz[k];

                g_x_0_xyzz_yyyy[k] = -g_x_0_xzz_yyyy[k] * ab_y + g_x_0_xzz_yyyyy[k];

                g_x_0_xyzz_yyyz[k] = -g_x_0_xzz_yyyz[k] * ab_y + g_x_0_xzz_yyyyz[k];

                g_x_0_xyzz_yyzz[k] = -g_x_0_xzz_yyzz[k] * ab_y + g_x_0_xzz_yyyzz[k];

                g_x_0_xyzz_yzzz[k] = -g_x_0_xzz_yzzz[k] * ab_y + g_x_0_xzz_yyzzz[k];

                g_x_0_xyzz_zzzz[k] = -g_x_0_xzz_zzzz[k] * ab_y + g_x_0_xzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzz_xxxx, g_x_0_xzz_xxxxz, g_x_0_xzz_xxxy, g_x_0_xzz_xxxyz, g_x_0_xzz_xxxz, g_x_0_xzz_xxxzz, g_x_0_xzz_xxyy, g_x_0_xzz_xxyyz, g_x_0_xzz_xxyz, g_x_0_xzz_xxyzz, g_x_0_xzz_xxzz, g_x_0_xzz_xxzzz, g_x_0_xzz_xyyy, g_x_0_xzz_xyyyz, g_x_0_xzz_xyyz, g_x_0_xzz_xyyzz, g_x_0_xzz_xyzz, g_x_0_xzz_xyzzz, g_x_0_xzz_xzzz, g_x_0_xzz_xzzzz, g_x_0_xzz_yyyy, g_x_0_xzz_yyyyz, g_x_0_xzz_yyyz, g_x_0_xzz_yyyzz, g_x_0_xzz_yyzz, g_x_0_xzz_yyzzz, g_x_0_xzz_yzzz, g_x_0_xzz_yzzzz, g_x_0_xzz_zzzz, g_x_0_xzz_zzzzz, g_x_0_xzzz_xxxx, g_x_0_xzzz_xxxy, g_x_0_xzzz_xxxz, g_x_0_xzzz_xxyy, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxzz, g_x_0_xzzz_xyyy, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xzzz, g_x_0_xzzz_yyyy, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzz_xxxx[k] = -g_x_0_xzz_xxxx[k] * ab_z + g_x_0_xzz_xxxxz[k];

                g_x_0_xzzz_xxxy[k] = -g_x_0_xzz_xxxy[k] * ab_z + g_x_0_xzz_xxxyz[k];

                g_x_0_xzzz_xxxz[k] = -g_x_0_xzz_xxxz[k] * ab_z + g_x_0_xzz_xxxzz[k];

                g_x_0_xzzz_xxyy[k] = -g_x_0_xzz_xxyy[k] * ab_z + g_x_0_xzz_xxyyz[k];

                g_x_0_xzzz_xxyz[k] = -g_x_0_xzz_xxyz[k] * ab_z + g_x_0_xzz_xxyzz[k];

                g_x_0_xzzz_xxzz[k] = -g_x_0_xzz_xxzz[k] * ab_z + g_x_0_xzz_xxzzz[k];

                g_x_0_xzzz_xyyy[k] = -g_x_0_xzz_xyyy[k] * ab_z + g_x_0_xzz_xyyyz[k];

                g_x_0_xzzz_xyyz[k] = -g_x_0_xzz_xyyz[k] * ab_z + g_x_0_xzz_xyyzz[k];

                g_x_0_xzzz_xyzz[k] = -g_x_0_xzz_xyzz[k] * ab_z + g_x_0_xzz_xyzzz[k];

                g_x_0_xzzz_xzzz[k] = -g_x_0_xzz_xzzz[k] * ab_z + g_x_0_xzz_xzzzz[k];

                g_x_0_xzzz_yyyy[k] = -g_x_0_xzz_yyyy[k] * ab_z + g_x_0_xzz_yyyyz[k];

                g_x_0_xzzz_yyyz[k] = -g_x_0_xzz_yyyz[k] * ab_z + g_x_0_xzz_yyyzz[k];

                g_x_0_xzzz_yyzz[k] = -g_x_0_xzz_yyzz[k] * ab_z + g_x_0_xzz_yyzzz[k];

                g_x_0_xzzz_yzzz[k] = -g_x_0_xzz_yzzz[k] * ab_z + g_x_0_xzz_yzzzz[k];

                g_x_0_xzzz_zzzz[k] = -g_x_0_xzz_zzzz[k] * ab_z + g_x_0_xzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyy_xxxx, g_x_0_yyy_xxxxy, g_x_0_yyy_xxxy, g_x_0_yyy_xxxyy, g_x_0_yyy_xxxyz, g_x_0_yyy_xxxz, g_x_0_yyy_xxyy, g_x_0_yyy_xxyyy, g_x_0_yyy_xxyyz, g_x_0_yyy_xxyz, g_x_0_yyy_xxyzz, g_x_0_yyy_xxzz, g_x_0_yyy_xyyy, g_x_0_yyy_xyyyy, g_x_0_yyy_xyyyz, g_x_0_yyy_xyyz, g_x_0_yyy_xyyzz, g_x_0_yyy_xyzz, g_x_0_yyy_xyzzz, g_x_0_yyy_xzzz, g_x_0_yyy_yyyy, g_x_0_yyy_yyyyy, g_x_0_yyy_yyyyz, g_x_0_yyy_yyyz, g_x_0_yyy_yyyzz, g_x_0_yyy_yyzz, g_x_0_yyy_yyzzz, g_x_0_yyy_yzzz, g_x_0_yyy_yzzzz, g_x_0_yyy_zzzz, g_x_0_yyyy_xxxx, g_x_0_yyyy_xxxy, g_x_0_yyyy_xxxz, g_x_0_yyyy_xxyy, g_x_0_yyyy_xxyz, g_x_0_yyyy_xxzz, g_x_0_yyyy_xyyy, g_x_0_yyyy_xyyz, g_x_0_yyyy_xyzz, g_x_0_yyyy_xzzz, g_x_0_yyyy_yyyy, g_x_0_yyyy_yyyz, g_x_0_yyyy_yyzz, g_x_0_yyyy_yzzz, g_x_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyy_xxxx[k] = -g_x_0_yyy_xxxx[k] * ab_y + g_x_0_yyy_xxxxy[k];

                g_x_0_yyyy_xxxy[k] = -g_x_0_yyy_xxxy[k] * ab_y + g_x_0_yyy_xxxyy[k];

                g_x_0_yyyy_xxxz[k] = -g_x_0_yyy_xxxz[k] * ab_y + g_x_0_yyy_xxxyz[k];

                g_x_0_yyyy_xxyy[k] = -g_x_0_yyy_xxyy[k] * ab_y + g_x_0_yyy_xxyyy[k];

                g_x_0_yyyy_xxyz[k] = -g_x_0_yyy_xxyz[k] * ab_y + g_x_0_yyy_xxyyz[k];

                g_x_0_yyyy_xxzz[k] = -g_x_0_yyy_xxzz[k] * ab_y + g_x_0_yyy_xxyzz[k];

                g_x_0_yyyy_xyyy[k] = -g_x_0_yyy_xyyy[k] * ab_y + g_x_0_yyy_xyyyy[k];

                g_x_0_yyyy_xyyz[k] = -g_x_0_yyy_xyyz[k] * ab_y + g_x_0_yyy_xyyyz[k];

                g_x_0_yyyy_xyzz[k] = -g_x_0_yyy_xyzz[k] * ab_y + g_x_0_yyy_xyyzz[k];

                g_x_0_yyyy_xzzz[k] = -g_x_0_yyy_xzzz[k] * ab_y + g_x_0_yyy_xyzzz[k];

                g_x_0_yyyy_yyyy[k] = -g_x_0_yyy_yyyy[k] * ab_y + g_x_0_yyy_yyyyy[k];

                g_x_0_yyyy_yyyz[k] = -g_x_0_yyy_yyyz[k] * ab_y + g_x_0_yyy_yyyyz[k];

                g_x_0_yyyy_yyzz[k] = -g_x_0_yyy_yyzz[k] * ab_y + g_x_0_yyy_yyyzz[k];

                g_x_0_yyyy_yzzz[k] = -g_x_0_yyy_yzzz[k] * ab_y + g_x_0_yyy_yyzzz[k];

                g_x_0_yyyy_zzzz[k] = -g_x_0_yyy_zzzz[k] * ab_y + g_x_0_yyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyz_xxxx, g_x_0_yyyz_xxxy, g_x_0_yyyz_xxxz, g_x_0_yyyz_xxyy, g_x_0_yyyz_xxyz, g_x_0_yyyz_xxzz, g_x_0_yyyz_xyyy, g_x_0_yyyz_xyyz, g_x_0_yyyz_xyzz, g_x_0_yyyz_xzzz, g_x_0_yyyz_yyyy, g_x_0_yyyz_yyyz, g_x_0_yyyz_yyzz, g_x_0_yyyz_yzzz, g_x_0_yyyz_zzzz, g_x_0_yyz_xxxx, g_x_0_yyz_xxxxy, g_x_0_yyz_xxxy, g_x_0_yyz_xxxyy, g_x_0_yyz_xxxyz, g_x_0_yyz_xxxz, g_x_0_yyz_xxyy, g_x_0_yyz_xxyyy, g_x_0_yyz_xxyyz, g_x_0_yyz_xxyz, g_x_0_yyz_xxyzz, g_x_0_yyz_xxzz, g_x_0_yyz_xyyy, g_x_0_yyz_xyyyy, g_x_0_yyz_xyyyz, g_x_0_yyz_xyyz, g_x_0_yyz_xyyzz, g_x_0_yyz_xyzz, g_x_0_yyz_xyzzz, g_x_0_yyz_xzzz, g_x_0_yyz_yyyy, g_x_0_yyz_yyyyy, g_x_0_yyz_yyyyz, g_x_0_yyz_yyyz, g_x_0_yyz_yyyzz, g_x_0_yyz_yyzz, g_x_0_yyz_yyzzz, g_x_0_yyz_yzzz, g_x_0_yyz_yzzzz, g_x_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyz_xxxx[k] = -g_x_0_yyz_xxxx[k] * ab_y + g_x_0_yyz_xxxxy[k];

                g_x_0_yyyz_xxxy[k] = -g_x_0_yyz_xxxy[k] * ab_y + g_x_0_yyz_xxxyy[k];

                g_x_0_yyyz_xxxz[k] = -g_x_0_yyz_xxxz[k] * ab_y + g_x_0_yyz_xxxyz[k];

                g_x_0_yyyz_xxyy[k] = -g_x_0_yyz_xxyy[k] * ab_y + g_x_0_yyz_xxyyy[k];

                g_x_0_yyyz_xxyz[k] = -g_x_0_yyz_xxyz[k] * ab_y + g_x_0_yyz_xxyyz[k];

                g_x_0_yyyz_xxzz[k] = -g_x_0_yyz_xxzz[k] * ab_y + g_x_0_yyz_xxyzz[k];

                g_x_0_yyyz_xyyy[k] = -g_x_0_yyz_xyyy[k] * ab_y + g_x_0_yyz_xyyyy[k];

                g_x_0_yyyz_xyyz[k] = -g_x_0_yyz_xyyz[k] * ab_y + g_x_0_yyz_xyyyz[k];

                g_x_0_yyyz_xyzz[k] = -g_x_0_yyz_xyzz[k] * ab_y + g_x_0_yyz_xyyzz[k];

                g_x_0_yyyz_xzzz[k] = -g_x_0_yyz_xzzz[k] * ab_y + g_x_0_yyz_xyzzz[k];

                g_x_0_yyyz_yyyy[k] = -g_x_0_yyz_yyyy[k] * ab_y + g_x_0_yyz_yyyyy[k];

                g_x_0_yyyz_yyyz[k] = -g_x_0_yyz_yyyz[k] * ab_y + g_x_0_yyz_yyyyz[k];

                g_x_0_yyyz_yyzz[k] = -g_x_0_yyz_yyzz[k] * ab_y + g_x_0_yyz_yyyzz[k];

                g_x_0_yyyz_yzzz[k] = -g_x_0_yyz_yzzz[k] * ab_y + g_x_0_yyz_yyzzz[k];

                g_x_0_yyyz_zzzz[k] = -g_x_0_yyz_zzzz[k] * ab_y + g_x_0_yyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzz_xxxx, g_x_0_yyzz_xxxy, g_x_0_yyzz_xxxz, g_x_0_yyzz_xxyy, g_x_0_yyzz_xxyz, g_x_0_yyzz_xxzz, g_x_0_yyzz_xyyy, g_x_0_yyzz_xyyz, g_x_0_yyzz_xyzz, g_x_0_yyzz_xzzz, g_x_0_yyzz_yyyy, g_x_0_yyzz_yyyz, g_x_0_yyzz_yyzz, g_x_0_yyzz_yzzz, g_x_0_yyzz_zzzz, g_x_0_yzz_xxxx, g_x_0_yzz_xxxxy, g_x_0_yzz_xxxy, g_x_0_yzz_xxxyy, g_x_0_yzz_xxxyz, g_x_0_yzz_xxxz, g_x_0_yzz_xxyy, g_x_0_yzz_xxyyy, g_x_0_yzz_xxyyz, g_x_0_yzz_xxyz, g_x_0_yzz_xxyzz, g_x_0_yzz_xxzz, g_x_0_yzz_xyyy, g_x_0_yzz_xyyyy, g_x_0_yzz_xyyyz, g_x_0_yzz_xyyz, g_x_0_yzz_xyyzz, g_x_0_yzz_xyzz, g_x_0_yzz_xyzzz, g_x_0_yzz_xzzz, g_x_0_yzz_yyyy, g_x_0_yzz_yyyyy, g_x_0_yzz_yyyyz, g_x_0_yzz_yyyz, g_x_0_yzz_yyyzz, g_x_0_yzz_yyzz, g_x_0_yzz_yyzzz, g_x_0_yzz_yzzz, g_x_0_yzz_yzzzz, g_x_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzz_xxxx[k] = -g_x_0_yzz_xxxx[k] * ab_y + g_x_0_yzz_xxxxy[k];

                g_x_0_yyzz_xxxy[k] = -g_x_0_yzz_xxxy[k] * ab_y + g_x_0_yzz_xxxyy[k];

                g_x_0_yyzz_xxxz[k] = -g_x_0_yzz_xxxz[k] * ab_y + g_x_0_yzz_xxxyz[k];

                g_x_0_yyzz_xxyy[k] = -g_x_0_yzz_xxyy[k] * ab_y + g_x_0_yzz_xxyyy[k];

                g_x_0_yyzz_xxyz[k] = -g_x_0_yzz_xxyz[k] * ab_y + g_x_0_yzz_xxyyz[k];

                g_x_0_yyzz_xxzz[k] = -g_x_0_yzz_xxzz[k] * ab_y + g_x_0_yzz_xxyzz[k];

                g_x_0_yyzz_xyyy[k] = -g_x_0_yzz_xyyy[k] * ab_y + g_x_0_yzz_xyyyy[k];

                g_x_0_yyzz_xyyz[k] = -g_x_0_yzz_xyyz[k] * ab_y + g_x_0_yzz_xyyyz[k];

                g_x_0_yyzz_xyzz[k] = -g_x_0_yzz_xyzz[k] * ab_y + g_x_0_yzz_xyyzz[k];

                g_x_0_yyzz_xzzz[k] = -g_x_0_yzz_xzzz[k] * ab_y + g_x_0_yzz_xyzzz[k];

                g_x_0_yyzz_yyyy[k] = -g_x_0_yzz_yyyy[k] * ab_y + g_x_0_yzz_yyyyy[k];

                g_x_0_yyzz_yyyz[k] = -g_x_0_yzz_yyyz[k] * ab_y + g_x_0_yzz_yyyyz[k];

                g_x_0_yyzz_yyzz[k] = -g_x_0_yzz_yyzz[k] * ab_y + g_x_0_yzz_yyyzz[k];

                g_x_0_yyzz_yzzz[k] = -g_x_0_yzz_yzzz[k] * ab_y + g_x_0_yzz_yyzzz[k];

                g_x_0_yyzz_zzzz[k] = -g_x_0_yzz_zzzz[k] * ab_y + g_x_0_yzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzz_xxxx, g_x_0_yzzz_xxxy, g_x_0_yzzz_xxxz, g_x_0_yzzz_xxyy, g_x_0_yzzz_xxyz, g_x_0_yzzz_xxzz, g_x_0_yzzz_xyyy, g_x_0_yzzz_xyyz, g_x_0_yzzz_xyzz, g_x_0_yzzz_xzzz, g_x_0_yzzz_yyyy, g_x_0_yzzz_yyyz, g_x_0_yzzz_yyzz, g_x_0_yzzz_yzzz, g_x_0_yzzz_zzzz, g_x_0_zzz_xxxx, g_x_0_zzz_xxxxy, g_x_0_zzz_xxxy, g_x_0_zzz_xxxyy, g_x_0_zzz_xxxyz, g_x_0_zzz_xxxz, g_x_0_zzz_xxyy, g_x_0_zzz_xxyyy, g_x_0_zzz_xxyyz, g_x_0_zzz_xxyz, g_x_0_zzz_xxyzz, g_x_0_zzz_xxzz, g_x_0_zzz_xyyy, g_x_0_zzz_xyyyy, g_x_0_zzz_xyyyz, g_x_0_zzz_xyyz, g_x_0_zzz_xyyzz, g_x_0_zzz_xyzz, g_x_0_zzz_xyzzz, g_x_0_zzz_xzzz, g_x_0_zzz_yyyy, g_x_0_zzz_yyyyy, g_x_0_zzz_yyyyz, g_x_0_zzz_yyyz, g_x_0_zzz_yyyzz, g_x_0_zzz_yyzz, g_x_0_zzz_yyzzz, g_x_0_zzz_yzzz, g_x_0_zzz_yzzzz, g_x_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzz_xxxx[k] = -g_x_0_zzz_xxxx[k] * ab_y + g_x_0_zzz_xxxxy[k];

                g_x_0_yzzz_xxxy[k] = -g_x_0_zzz_xxxy[k] * ab_y + g_x_0_zzz_xxxyy[k];

                g_x_0_yzzz_xxxz[k] = -g_x_0_zzz_xxxz[k] * ab_y + g_x_0_zzz_xxxyz[k];

                g_x_0_yzzz_xxyy[k] = -g_x_0_zzz_xxyy[k] * ab_y + g_x_0_zzz_xxyyy[k];

                g_x_0_yzzz_xxyz[k] = -g_x_0_zzz_xxyz[k] * ab_y + g_x_0_zzz_xxyyz[k];

                g_x_0_yzzz_xxzz[k] = -g_x_0_zzz_xxzz[k] * ab_y + g_x_0_zzz_xxyzz[k];

                g_x_0_yzzz_xyyy[k] = -g_x_0_zzz_xyyy[k] * ab_y + g_x_0_zzz_xyyyy[k];

                g_x_0_yzzz_xyyz[k] = -g_x_0_zzz_xyyz[k] * ab_y + g_x_0_zzz_xyyyz[k];

                g_x_0_yzzz_xyzz[k] = -g_x_0_zzz_xyzz[k] * ab_y + g_x_0_zzz_xyyzz[k];

                g_x_0_yzzz_xzzz[k] = -g_x_0_zzz_xzzz[k] * ab_y + g_x_0_zzz_xyzzz[k];

                g_x_0_yzzz_yyyy[k] = -g_x_0_zzz_yyyy[k] * ab_y + g_x_0_zzz_yyyyy[k];

                g_x_0_yzzz_yyyz[k] = -g_x_0_zzz_yyyz[k] * ab_y + g_x_0_zzz_yyyyz[k];

                g_x_0_yzzz_yyzz[k] = -g_x_0_zzz_yyzz[k] * ab_y + g_x_0_zzz_yyyzz[k];

                g_x_0_yzzz_yzzz[k] = -g_x_0_zzz_yzzz[k] * ab_y + g_x_0_zzz_yyzzz[k];

                g_x_0_yzzz_zzzz[k] = -g_x_0_zzz_zzzz[k] * ab_y + g_x_0_zzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzz_xxxx, g_x_0_zzz_xxxxz, g_x_0_zzz_xxxy, g_x_0_zzz_xxxyz, g_x_0_zzz_xxxz, g_x_0_zzz_xxxzz, g_x_0_zzz_xxyy, g_x_0_zzz_xxyyz, g_x_0_zzz_xxyz, g_x_0_zzz_xxyzz, g_x_0_zzz_xxzz, g_x_0_zzz_xxzzz, g_x_0_zzz_xyyy, g_x_0_zzz_xyyyz, g_x_0_zzz_xyyz, g_x_0_zzz_xyyzz, g_x_0_zzz_xyzz, g_x_0_zzz_xyzzz, g_x_0_zzz_xzzz, g_x_0_zzz_xzzzz, g_x_0_zzz_yyyy, g_x_0_zzz_yyyyz, g_x_0_zzz_yyyz, g_x_0_zzz_yyyzz, g_x_0_zzz_yyzz, g_x_0_zzz_yyzzz, g_x_0_zzz_yzzz, g_x_0_zzz_yzzzz, g_x_0_zzz_zzzz, g_x_0_zzz_zzzzz, g_x_0_zzzz_xxxx, g_x_0_zzzz_xxxy, g_x_0_zzzz_xxxz, g_x_0_zzzz_xxyy, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxzz, g_x_0_zzzz_xyyy, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xzzz, g_x_0_zzzz_yyyy, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzz_xxxx[k] = -g_x_0_zzz_xxxx[k] * ab_z + g_x_0_zzz_xxxxz[k];

                g_x_0_zzzz_xxxy[k] = -g_x_0_zzz_xxxy[k] * ab_z + g_x_0_zzz_xxxyz[k];

                g_x_0_zzzz_xxxz[k] = -g_x_0_zzz_xxxz[k] * ab_z + g_x_0_zzz_xxxzz[k];

                g_x_0_zzzz_xxyy[k] = -g_x_0_zzz_xxyy[k] * ab_z + g_x_0_zzz_xxyyz[k];

                g_x_0_zzzz_xxyz[k] = -g_x_0_zzz_xxyz[k] * ab_z + g_x_0_zzz_xxyzz[k];

                g_x_0_zzzz_xxzz[k] = -g_x_0_zzz_xxzz[k] * ab_z + g_x_0_zzz_xxzzz[k];

                g_x_0_zzzz_xyyy[k] = -g_x_0_zzz_xyyy[k] * ab_z + g_x_0_zzz_xyyyz[k];

                g_x_0_zzzz_xyyz[k] = -g_x_0_zzz_xyyz[k] * ab_z + g_x_0_zzz_xyyzz[k];

                g_x_0_zzzz_xyzz[k] = -g_x_0_zzz_xyzz[k] * ab_z + g_x_0_zzz_xyzzz[k];

                g_x_0_zzzz_xzzz[k] = -g_x_0_zzz_xzzz[k] * ab_z + g_x_0_zzz_xzzzz[k];

                g_x_0_zzzz_yyyy[k] = -g_x_0_zzz_yyyy[k] * ab_z + g_x_0_zzz_yyyyz[k];

                g_x_0_zzzz_yyyz[k] = -g_x_0_zzz_yyyz[k] * ab_z + g_x_0_zzz_yyyzz[k];

                g_x_0_zzzz_yyzz[k] = -g_x_0_zzz_yyzz[k] * ab_z + g_x_0_zzz_yyzzz[k];

                g_x_0_zzzz_yzzz[k] = -g_x_0_zzz_yzzz[k] * ab_z + g_x_0_zzz_yzzzz[k];

                g_x_0_zzzz_zzzz[k] = -g_x_0_zzz_zzzz[k] * ab_z + g_x_0_zzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxx_xxxx, g_y_0_xxx_xxxxx, g_y_0_xxx_xxxxy, g_y_0_xxx_xxxxz, g_y_0_xxx_xxxy, g_y_0_xxx_xxxyy, g_y_0_xxx_xxxyz, g_y_0_xxx_xxxz, g_y_0_xxx_xxxzz, g_y_0_xxx_xxyy, g_y_0_xxx_xxyyy, g_y_0_xxx_xxyyz, g_y_0_xxx_xxyz, g_y_0_xxx_xxyzz, g_y_0_xxx_xxzz, g_y_0_xxx_xxzzz, g_y_0_xxx_xyyy, g_y_0_xxx_xyyyy, g_y_0_xxx_xyyyz, g_y_0_xxx_xyyz, g_y_0_xxx_xyyzz, g_y_0_xxx_xyzz, g_y_0_xxx_xyzzz, g_y_0_xxx_xzzz, g_y_0_xxx_xzzzz, g_y_0_xxx_yyyy, g_y_0_xxx_yyyz, g_y_0_xxx_yyzz, g_y_0_xxx_yzzz, g_y_0_xxx_zzzz, g_y_0_xxxx_xxxx, g_y_0_xxxx_xxxy, g_y_0_xxxx_xxxz, g_y_0_xxxx_xxyy, g_y_0_xxxx_xxyz, g_y_0_xxxx_xxzz, g_y_0_xxxx_xyyy, g_y_0_xxxx_xyyz, g_y_0_xxxx_xyzz, g_y_0_xxxx_xzzz, g_y_0_xxxx_yyyy, g_y_0_xxxx_yyyz, g_y_0_xxxx_yyzz, g_y_0_xxxx_yzzz, g_y_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxx_xxxx[k] = -g_y_0_xxx_xxxx[k] * ab_x + g_y_0_xxx_xxxxx[k];

                g_y_0_xxxx_xxxy[k] = -g_y_0_xxx_xxxy[k] * ab_x + g_y_0_xxx_xxxxy[k];

                g_y_0_xxxx_xxxz[k] = -g_y_0_xxx_xxxz[k] * ab_x + g_y_0_xxx_xxxxz[k];

                g_y_0_xxxx_xxyy[k] = -g_y_0_xxx_xxyy[k] * ab_x + g_y_0_xxx_xxxyy[k];

                g_y_0_xxxx_xxyz[k] = -g_y_0_xxx_xxyz[k] * ab_x + g_y_0_xxx_xxxyz[k];

                g_y_0_xxxx_xxzz[k] = -g_y_0_xxx_xxzz[k] * ab_x + g_y_0_xxx_xxxzz[k];

                g_y_0_xxxx_xyyy[k] = -g_y_0_xxx_xyyy[k] * ab_x + g_y_0_xxx_xxyyy[k];

                g_y_0_xxxx_xyyz[k] = -g_y_0_xxx_xyyz[k] * ab_x + g_y_0_xxx_xxyyz[k];

                g_y_0_xxxx_xyzz[k] = -g_y_0_xxx_xyzz[k] * ab_x + g_y_0_xxx_xxyzz[k];

                g_y_0_xxxx_xzzz[k] = -g_y_0_xxx_xzzz[k] * ab_x + g_y_0_xxx_xxzzz[k];

                g_y_0_xxxx_yyyy[k] = -g_y_0_xxx_yyyy[k] * ab_x + g_y_0_xxx_xyyyy[k];

                g_y_0_xxxx_yyyz[k] = -g_y_0_xxx_yyyz[k] * ab_x + g_y_0_xxx_xyyyz[k];

                g_y_0_xxxx_yyzz[k] = -g_y_0_xxx_yyzz[k] * ab_x + g_y_0_xxx_xyyzz[k];

                g_y_0_xxxx_yzzz[k] = -g_y_0_xxx_yzzz[k] * ab_x + g_y_0_xxx_xyzzz[k];

                g_y_0_xxxx_zzzz[k] = -g_y_0_xxx_zzzz[k] * ab_x + g_y_0_xxx_xzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 251 * ccomps * dcomps);

            auto g_y_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxy_xxxx, g_y_0_xxxy_xxxy, g_y_0_xxxy_xxxz, g_y_0_xxxy_xxyy, g_y_0_xxxy_xxyz, g_y_0_xxxy_xxzz, g_y_0_xxxy_xyyy, g_y_0_xxxy_xyyz, g_y_0_xxxy_xyzz, g_y_0_xxxy_xzzz, g_y_0_xxxy_yyyy, g_y_0_xxxy_yyyz, g_y_0_xxxy_yyzz, g_y_0_xxxy_yzzz, g_y_0_xxxy_zzzz, g_y_0_xxy_xxxx, g_y_0_xxy_xxxxx, g_y_0_xxy_xxxxy, g_y_0_xxy_xxxxz, g_y_0_xxy_xxxy, g_y_0_xxy_xxxyy, g_y_0_xxy_xxxyz, g_y_0_xxy_xxxz, g_y_0_xxy_xxxzz, g_y_0_xxy_xxyy, g_y_0_xxy_xxyyy, g_y_0_xxy_xxyyz, g_y_0_xxy_xxyz, g_y_0_xxy_xxyzz, g_y_0_xxy_xxzz, g_y_0_xxy_xxzzz, g_y_0_xxy_xyyy, g_y_0_xxy_xyyyy, g_y_0_xxy_xyyyz, g_y_0_xxy_xyyz, g_y_0_xxy_xyyzz, g_y_0_xxy_xyzz, g_y_0_xxy_xyzzz, g_y_0_xxy_xzzz, g_y_0_xxy_xzzzz, g_y_0_xxy_yyyy, g_y_0_xxy_yyyz, g_y_0_xxy_yyzz, g_y_0_xxy_yzzz, g_y_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxy_xxxx[k] = -g_y_0_xxy_xxxx[k] * ab_x + g_y_0_xxy_xxxxx[k];

                g_y_0_xxxy_xxxy[k] = -g_y_0_xxy_xxxy[k] * ab_x + g_y_0_xxy_xxxxy[k];

                g_y_0_xxxy_xxxz[k] = -g_y_0_xxy_xxxz[k] * ab_x + g_y_0_xxy_xxxxz[k];

                g_y_0_xxxy_xxyy[k] = -g_y_0_xxy_xxyy[k] * ab_x + g_y_0_xxy_xxxyy[k];

                g_y_0_xxxy_xxyz[k] = -g_y_0_xxy_xxyz[k] * ab_x + g_y_0_xxy_xxxyz[k];

                g_y_0_xxxy_xxzz[k] = -g_y_0_xxy_xxzz[k] * ab_x + g_y_0_xxy_xxxzz[k];

                g_y_0_xxxy_xyyy[k] = -g_y_0_xxy_xyyy[k] * ab_x + g_y_0_xxy_xxyyy[k];

                g_y_0_xxxy_xyyz[k] = -g_y_0_xxy_xyyz[k] * ab_x + g_y_0_xxy_xxyyz[k];

                g_y_0_xxxy_xyzz[k] = -g_y_0_xxy_xyzz[k] * ab_x + g_y_0_xxy_xxyzz[k];

                g_y_0_xxxy_xzzz[k] = -g_y_0_xxy_xzzz[k] * ab_x + g_y_0_xxy_xxzzz[k];

                g_y_0_xxxy_yyyy[k] = -g_y_0_xxy_yyyy[k] * ab_x + g_y_0_xxy_xyyyy[k];

                g_y_0_xxxy_yyyz[k] = -g_y_0_xxy_yyyz[k] * ab_x + g_y_0_xxy_xyyyz[k];

                g_y_0_xxxy_yyzz[k] = -g_y_0_xxy_yyzz[k] * ab_x + g_y_0_xxy_xyyzz[k];

                g_y_0_xxxy_yzzz[k] = -g_y_0_xxy_yzzz[k] * ab_x + g_y_0_xxy_xyzzz[k];

                g_y_0_xxxy_zzzz[k] = -g_y_0_xxy_zzzz[k] * ab_x + g_y_0_xxy_xzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxz_xxxx, g_y_0_xxxz_xxxy, g_y_0_xxxz_xxxz, g_y_0_xxxz_xxyy, g_y_0_xxxz_xxyz, g_y_0_xxxz_xxzz, g_y_0_xxxz_xyyy, g_y_0_xxxz_xyyz, g_y_0_xxxz_xyzz, g_y_0_xxxz_xzzz, g_y_0_xxxz_yyyy, g_y_0_xxxz_yyyz, g_y_0_xxxz_yyzz, g_y_0_xxxz_yzzz, g_y_0_xxxz_zzzz, g_y_0_xxz_xxxx, g_y_0_xxz_xxxxx, g_y_0_xxz_xxxxy, g_y_0_xxz_xxxxz, g_y_0_xxz_xxxy, g_y_0_xxz_xxxyy, g_y_0_xxz_xxxyz, g_y_0_xxz_xxxz, g_y_0_xxz_xxxzz, g_y_0_xxz_xxyy, g_y_0_xxz_xxyyy, g_y_0_xxz_xxyyz, g_y_0_xxz_xxyz, g_y_0_xxz_xxyzz, g_y_0_xxz_xxzz, g_y_0_xxz_xxzzz, g_y_0_xxz_xyyy, g_y_0_xxz_xyyyy, g_y_0_xxz_xyyyz, g_y_0_xxz_xyyz, g_y_0_xxz_xyyzz, g_y_0_xxz_xyzz, g_y_0_xxz_xyzzz, g_y_0_xxz_xzzz, g_y_0_xxz_xzzzz, g_y_0_xxz_yyyy, g_y_0_xxz_yyyz, g_y_0_xxz_yyzz, g_y_0_xxz_yzzz, g_y_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxz_xxxx[k] = -g_y_0_xxz_xxxx[k] * ab_x + g_y_0_xxz_xxxxx[k];

                g_y_0_xxxz_xxxy[k] = -g_y_0_xxz_xxxy[k] * ab_x + g_y_0_xxz_xxxxy[k];

                g_y_0_xxxz_xxxz[k] = -g_y_0_xxz_xxxz[k] * ab_x + g_y_0_xxz_xxxxz[k];

                g_y_0_xxxz_xxyy[k] = -g_y_0_xxz_xxyy[k] * ab_x + g_y_0_xxz_xxxyy[k];

                g_y_0_xxxz_xxyz[k] = -g_y_0_xxz_xxyz[k] * ab_x + g_y_0_xxz_xxxyz[k];

                g_y_0_xxxz_xxzz[k] = -g_y_0_xxz_xxzz[k] * ab_x + g_y_0_xxz_xxxzz[k];

                g_y_0_xxxz_xyyy[k] = -g_y_0_xxz_xyyy[k] * ab_x + g_y_0_xxz_xxyyy[k];

                g_y_0_xxxz_xyyz[k] = -g_y_0_xxz_xyyz[k] * ab_x + g_y_0_xxz_xxyyz[k];

                g_y_0_xxxz_xyzz[k] = -g_y_0_xxz_xyzz[k] * ab_x + g_y_0_xxz_xxyzz[k];

                g_y_0_xxxz_xzzz[k] = -g_y_0_xxz_xzzz[k] * ab_x + g_y_0_xxz_xxzzz[k];

                g_y_0_xxxz_yyyy[k] = -g_y_0_xxz_yyyy[k] * ab_x + g_y_0_xxz_xyyyy[k];

                g_y_0_xxxz_yyyz[k] = -g_y_0_xxz_yyyz[k] * ab_x + g_y_0_xxz_xyyyz[k];

                g_y_0_xxxz_yyzz[k] = -g_y_0_xxz_yyzz[k] * ab_x + g_y_0_xxz_xyyzz[k];

                g_y_0_xxxz_yzzz[k] = -g_y_0_xxz_yzzz[k] * ab_x + g_y_0_xxz_xyzzz[k];

                g_y_0_xxxz_zzzz[k] = -g_y_0_xxz_zzzz[k] * ab_x + g_y_0_xxz_xzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyy_xxxx, g_y_0_xxyy_xxxy, g_y_0_xxyy_xxxz, g_y_0_xxyy_xxyy, g_y_0_xxyy_xxyz, g_y_0_xxyy_xxzz, g_y_0_xxyy_xyyy, g_y_0_xxyy_xyyz, g_y_0_xxyy_xyzz, g_y_0_xxyy_xzzz, g_y_0_xxyy_yyyy, g_y_0_xxyy_yyyz, g_y_0_xxyy_yyzz, g_y_0_xxyy_yzzz, g_y_0_xxyy_zzzz, g_y_0_xyy_xxxx, g_y_0_xyy_xxxxx, g_y_0_xyy_xxxxy, g_y_0_xyy_xxxxz, g_y_0_xyy_xxxy, g_y_0_xyy_xxxyy, g_y_0_xyy_xxxyz, g_y_0_xyy_xxxz, g_y_0_xyy_xxxzz, g_y_0_xyy_xxyy, g_y_0_xyy_xxyyy, g_y_0_xyy_xxyyz, g_y_0_xyy_xxyz, g_y_0_xyy_xxyzz, g_y_0_xyy_xxzz, g_y_0_xyy_xxzzz, g_y_0_xyy_xyyy, g_y_0_xyy_xyyyy, g_y_0_xyy_xyyyz, g_y_0_xyy_xyyz, g_y_0_xyy_xyyzz, g_y_0_xyy_xyzz, g_y_0_xyy_xyzzz, g_y_0_xyy_xzzz, g_y_0_xyy_xzzzz, g_y_0_xyy_yyyy, g_y_0_xyy_yyyz, g_y_0_xyy_yyzz, g_y_0_xyy_yzzz, g_y_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyy_xxxx[k] = -g_y_0_xyy_xxxx[k] * ab_x + g_y_0_xyy_xxxxx[k];

                g_y_0_xxyy_xxxy[k] = -g_y_0_xyy_xxxy[k] * ab_x + g_y_0_xyy_xxxxy[k];

                g_y_0_xxyy_xxxz[k] = -g_y_0_xyy_xxxz[k] * ab_x + g_y_0_xyy_xxxxz[k];

                g_y_0_xxyy_xxyy[k] = -g_y_0_xyy_xxyy[k] * ab_x + g_y_0_xyy_xxxyy[k];

                g_y_0_xxyy_xxyz[k] = -g_y_0_xyy_xxyz[k] * ab_x + g_y_0_xyy_xxxyz[k];

                g_y_0_xxyy_xxzz[k] = -g_y_0_xyy_xxzz[k] * ab_x + g_y_0_xyy_xxxzz[k];

                g_y_0_xxyy_xyyy[k] = -g_y_0_xyy_xyyy[k] * ab_x + g_y_0_xyy_xxyyy[k];

                g_y_0_xxyy_xyyz[k] = -g_y_0_xyy_xyyz[k] * ab_x + g_y_0_xyy_xxyyz[k];

                g_y_0_xxyy_xyzz[k] = -g_y_0_xyy_xyzz[k] * ab_x + g_y_0_xyy_xxyzz[k];

                g_y_0_xxyy_xzzz[k] = -g_y_0_xyy_xzzz[k] * ab_x + g_y_0_xyy_xxzzz[k];

                g_y_0_xxyy_yyyy[k] = -g_y_0_xyy_yyyy[k] * ab_x + g_y_0_xyy_xyyyy[k];

                g_y_0_xxyy_yyyz[k] = -g_y_0_xyy_yyyz[k] * ab_x + g_y_0_xyy_xyyyz[k];

                g_y_0_xxyy_yyzz[k] = -g_y_0_xyy_yyzz[k] * ab_x + g_y_0_xyy_xyyzz[k];

                g_y_0_xxyy_yzzz[k] = -g_y_0_xyy_yzzz[k] * ab_x + g_y_0_xyy_xyzzz[k];

                g_y_0_xxyy_zzzz[k] = -g_y_0_xyy_zzzz[k] * ab_x + g_y_0_xyy_xzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyz_xxxx, g_y_0_xxyz_xxxy, g_y_0_xxyz_xxxz, g_y_0_xxyz_xxyy, g_y_0_xxyz_xxyz, g_y_0_xxyz_xxzz, g_y_0_xxyz_xyyy, g_y_0_xxyz_xyyz, g_y_0_xxyz_xyzz, g_y_0_xxyz_xzzz, g_y_0_xxyz_yyyy, g_y_0_xxyz_yyyz, g_y_0_xxyz_yyzz, g_y_0_xxyz_yzzz, g_y_0_xxyz_zzzz, g_y_0_xyz_xxxx, g_y_0_xyz_xxxxx, g_y_0_xyz_xxxxy, g_y_0_xyz_xxxxz, g_y_0_xyz_xxxy, g_y_0_xyz_xxxyy, g_y_0_xyz_xxxyz, g_y_0_xyz_xxxz, g_y_0_xyz_xxxzz, g_y_0_xyz_xxyy, g_y_0_xyz_xxyyy, g_y_0_xyz_xxyyz, g_y_0_xyz_xxyz, g_y_0_xyz_xxyzz, g_y_0_xyz_xxzz, g_y_0_xyz_xxzzz, g_y_0_xyz_xyyy, g_y_0_xyz_xyyyy, g_y_0_xyz_xyyyz, g_y_0_xyz_xyyz, g_y_0_xyz_xyyzz, g_y_0_xyz_xyzz, g_y_0_xyz_xyzzz, g_y_0_xyz_xzzz, g_y_0_xyz_xzzzz, g_y_0_xyz_yyyy, g_y_0_xyz_yyyz, g_y_0_xyz_yyzz, g_y_0_xyz_yzzz, g_y_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyz_xxxx[k] = -g_y_0_xyz_xxxx[k] * ab_x + g_y_0_xyz_xxxxx[k];

                g_y_0_xxyz_xxxy[k] = -g_y_0_xyz_xxxy[k] * ab_x + g_y_0_xyz_xxxxy[k];

                g_y_0_xxyz_xxxz[k] = -g_y_0_xyz_xxxz[k] * ab_x + g_y_0_xyz_xxxxz[k];

                g_y_0_xxyz_xxyy[k] = -g_y_0_xyz_xxyy[k] * ab_x + g_y_0_xyz_xxxyy[k];

                g_y_0_xxyz_xxyz[k] = -g_y_0_xyz_xxyz[k] * ab_x + g_y_0_xyz_xxxyz[k];

                g_y_0_xxyz_xxzz[k] = -g_y_0_xyz_xxzz[k] * ab_x + g_y_0_xyz_xxxzz[k];

                g_y_0_xxyz_xyyy[k] = -g_y_0_xyz_xyyy[k] * ab_x + g_y_0_xyz_xxyyy[k];

                g_y_0_xxyz_xyyz[k] = -g_y_0_xyz_xyyz[k] * ab_x + g_y_0_xyz_xxyyz[k];

                g_y_0_xxyz_xyzz[k] = -g_y_0_xyz_xyzz[k] * ab_x + g_y_0_xyz_xxyzz[k];

                g_y_0_xxyz_xzzz[k] = -g_y_0_xyz_xzzz[k] * ab_x + g_y_0_xyz_xxzzz[k];

                g_y_0_xxyz_yyyy[k] = -g_y_0_xyz_yyyy[k] * ab_x + g_y_0_xyz_xyyyy[k];

                g_y_0_xxyz_yyyz[k] = -g_y_0_xyz_yyyz[k] * ab_x + g_y_0_xyz_xyyyz[k];

                g_y_0_xxyz_yyzz[k] = -g_y_0_xyz_yyzz[k] * ab_x + g_y_0_xyz_xyyzz[k];

                g_y_0_xxyz_yzzz[k] = -g_y_0_xyz_yzzz[k] * ab_x + g_y_0_xyz_xyzzz[k];

                g_y_0_xxyz_zzzz[k] = -g_y_0_xyz_zzzz[k] * ab_x + g_y_0_xyz_xzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzz_xxxx, g_y_0_xxzz_xxxy, g_y_0_xxzz_xxxz, g_y_0_xxzz_xxyy, g_y_0_xxzz_xxyz, g_y_0_xxzz_xxzz, g_y_0_xxzz_xyyy, g_y_0_xxzz_xyyz, g_y_0_xxzz_xyzz, g_y_0_xxzz_xzzz, g_y_0_xxzz_yyyy, g_y_0_xxzz_yyyz, g_y_0_xxzz_yyzz, g_y_0_xxzz_yzzz, g_y_0_xxzz_zzzz, g_y_0_xzz_xxxx, g_y_0_xzz_xxxxx, g_y_0_xzz_xxxxy, g_y_0_xzz_xxxxz, g_y_0_xzz_xxxy, g_y_0_xzz_xxxyy, g_y_0_xzz_xxxyz, g_y_0_xzz_xxxz, g_y_0_xzz_xxxzz, g_y_0_xzz_xxyy, g_y_0_xzz_xxyyy, g_y_0_xzz_xxyyz, g_y_0_xzz_xxyz, g_y_0_xzz_xxyzz, g_y_0_xzz_xxzz, g_y_0_xzz_xxzzz, g_y_0_xzz_xyyy, g_y_0_xzz_xyyyy, g_y_0_xzz_xyyyz, g_y_0_xzz_xyyz, g_y_0_xzz_xyyzz, g_y_0_xzz_xyzz, g_y_0_xzz_xyzzz, g_y_0_xzz_xzzz, g_y_0_xzz_xzzzz, g_y_0_xzz_yyyy, g_y_0_xzz_yyyz, g_y_0_xzz_yyzz, g_y_0_xzz_yzzz, g_y_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzz_xxxx[k] = -g_y_0_xzz_xxxx[k] * ab_x + g_y_0_xzz_xxxxx[k];

                g_y_0_xxzz_xxxy[k] = -g_y_0_xzz_xxxy[k] * ab_x + g_y_0_xzz_xxxxy[k];

                g_y_0_xxzz_xxxz[k] = -g_y_0_xzz_xxxz[k] * ab_x + g_y_0_xzz_xxxxz[k];

                g_y_0_xxzz_xxyy[k] = -g_y_0_xzz_xxyy[k] * ab_x + g_y_0_xzz_xxxyy[k];

                g_y_0_xxzz_xxyz[k] = -g_y_0_xzz_xxyz[k] * ab_x + g_y_0_xzz_xxxyz[k];

                g_y_0_xxzz_xxzz[k] = -g_y_0_xzz_xxzz[k] * ab_x + g_y_0_xzz_xxxzz[k];

                g_y_0_xxzz_xyyy[k] = -g_y_0_xzz_xyyy[k] * ab_x + g_y_0_xzz_xxyyy[k];

                g_y_0_xxzz_xyyz[k] = -g_y_0_xzz_xyyz[k] * ab_x + g_y_0_xzz_xxyyz[k];

                g_y_0_xxzz_xyzz[k] = -g_y_0_xzz_xyzz[k] * ab_x + g_y_0_xzz_xxyzz[k];

                g_y_0_xxzz_xzzz[k] = -g_y_0_xzz_xzzz[k] * ab_x + g_y_0_xzz_xxzzz[k];

                g_y_0_xxzz_yyyy[k] = -g_y_0_xzz_yyyy[k] * ab_x + g_y_0_xzz_xyyyy[k];

                g_y_0_xxzz_yyyz[k] = -g_y_0_xzz_yyyz[k] * ab_x + g_y_0_xzz_xyyyz[k];

                g_y_0_xxzz_yyzz[k] = -g_y_0_xzz_yyzz[k] * ab_x + g_y_0_xzz_xyyzz[k];

                g_y_0_xxzz_yzzz[k] = -g_y_0_xzz_yzzz[k] * ab_x + g_y_0_xzz_xyzzz[k];

                g_y_0_xxzz_zzzz[k] = -g_y_0_xzz_zzzz[k] * ab_x + g_y_0_xzz_xzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyy_xxxx, g_y_0_xyyy_xxxy, g_y_0_xyyy_xxxz, g_y_0_xyyy_xxyy, g_y_0_xyyy_xxyz, g_y_0_xyyy_xxzz, g_y_0_xyyy_xyyy, g_y_0_xyyy_xyyz, g_y_0_xyyy_xyzz, g_y_0_xyyy_xzzz, g_y_0_xyyy_yyyy, g_y_0_xyyy_yyyz, g_y_0_xyyy_yyzz, g_y_0_xyyy_yzzz, g_y_0_xyyy_zzzz, g_y_0_yyy_xxxx, g_y_0_yyy_xxxxx, g_y_0_yyy_xxxxy, g_y_0_yyy_xxxxz, g_y_0_yyy_xxxy, g_y_0_yyy_xxxyy, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxz, g_y_0_yyy_xxxzz, g_y_0_yyy_xxyy, g_y_0_yyy_xxyyy, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxzz, g_y_0_yyy_xxzzz, g_y_0_yyy_xyyy, g_y_0_yyy_xyyyy, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xzzz, g_y_0_yyy_xzzzz, g_y_0_yyy_yyyy, g_y_0_yyy_yyyz, g_y_0_yyy_yyzz, g_y_0_yyy_yzzz, g_y_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyy_xxxx[k] = -g_y_0_yyy_xxxx[k] * ab_x + g_y_0_yyy_xxxxx[k];

                g_y_0_xyyy_xxxy[k] = -g_y_0_yyy_xxxy[k] * ab_x + g_y_0_yyy_xxxxy[k];

                g_y_0_xyyy_xxxz[k] = -g_y_0_yyy_xxxz[k] * ab_x + g_y_0_yyy_xxxxz[k];

                g_y_0_xyyy_xxyy[k] = -g_y_0_yyy_xxyy[k] * ab_x + g_y_0_yyy_xxxyy[k];

                g_y_0_xyyy_xxyz[k] = -g_y_0_yyy_xxyz[k] * ab_x + g_y_0_yyy_xxxyz[k];

                g_y_0_xyyy_xxzz[k] = -g_y_0_yyy_xxzz[k] * ab_x + g_y_0_yyy_xxxzz[k];

                g_y_0_xyyy_xyyy[k] = -g_y_0_yyy_xyyy[k] * ab_x + g_y_0_yyy_xxyyy[k];

                g_y_0_xyyy_xyyz[k] = -g_y_0_yyy_xyyz[k] * ab_x + g_y_0_yyy_xxyyz[k];

                g_y_0_xyyy_xyzz[k] = -g_y_0_yyy_xyzz[k] * ab_x + g_y_0_yyy_xxyzz[k];

                g_y_0_xyyy_xzzz[k] = -g_y_0_yyy_xzzz[k] * ab_x + g_y_0_yyy_xxzzz[k];

                g_y_0_xyyy_yyyy[k] = -g_y_0_yyy_yyyy[k] * ab_x + g_y_0_yyy_xyyyy[k];

                g_y_0_xyyy_yyyz[k] = -g_y_0_yyy_yyyz[k] * ab_x + g_y_0_yyy_xyyyz[k];

                g_y_0_xyyy_yyzz[k] = -g_y_0_yyy_yyzz[k] * ab_x + g_y_0_yyy_xyyzz[k];

                g_y_0_xyyy_yzzz[k] = -g_y_0_yyy_yzzz[k] * ab_x + g_y_0_yyy_xyzzz[k];

                g_y_0_xyyy_zzzz[k] = -g_y_0_yyy_zzzz[k] * ab_x + g_y_0_yyy_xzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyz_xxxx, g_y_0_xyyz_xxxy, g_y_0_xyyz_xxxz, g_y_0_xyyz_xxyy, g_y_0_xyyz_xxyz, g_y_0_xyyz_xxzz, g_y_0_xyyz_xyyy, g_y_0_xyyz_xyyz, g_y_0_xyyz_xyzz, g_y_0_xyyz_xzzz, g_y_0_xyyz_yyyy, g_y_0_xyyz_yyyz, g_y_0_xyyz_yyzz, g_y_0_xyyz_yzzz, g_y_0_xyyz_zzzz, g_y_0_yyz_xxxx, g_y_0_yyz_xxxxx, g_y_0_yyz_xxxxy, g_y_0_yyz_xxxxz, g_y_0_yyz_xxxy, g_y_0_yyz_xxxyy, g_y_0_yyz_xxxyz, g_y_0_yyz_xxxz, g_y_0_yyz_xxxzz, g_y_0_yyz_xxyy, g_y_0_yyz_xxyyy, g_y_0_yyz_xxyyz, g_y_0_yyz_xxyz, g_y_0_yyz_xxyzz, g_y_0_yyz_xxzz, g_y_0_yyz_xxzzz, g_y_0_yyz_xyyy, g_y_0_yyz_xyyyy, g_y_0_yyz_xyyyz, g_y_0_yyz_xyyz, g_y_0_yyz_xyyzz, g_y_0_yyz_xyzz, g_y_0_yyz_xyzzz, g_y_0_yyz_xzzz, g_y_0_yyz_xzzzz, g_y_0_yyz_yyyy, g_y_0_yyz_yyyz, g_y_0_yyz_yyzz, g_y_0_yyz_yzzz, g_y_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyz_xxxx[k] = -g_y_0_yyz_xxxx[k] * ab_x + g_y_0_yyz_xxxxx[k];

                g_y_0_xyyz_xxxy[k] = -g_y_0_yyz_xxxy[k] * ab_x + g_y_0_yyz_xxxxy[k];

                g_y_0_xyyz_xxxz[k] = -g_y_0_yyz_xxxz[k] * ab_x + g_y_0_yyz_xxxxz[k];

                g_y_0_xyyz_xxyy[k] = -g_y_0_yyz_xxyy[k] * ab_x + g_y_0_yyz_xxxyy[k];

                g_y_0_xyyz_xxyz[k] = -g_y_0_yyz_xxyz[k] * ab_x + g_y_0_yyz_xxxyz[k];

                g_y_0_xyyz_xxzz[k] = -g_y_0_yyz_xxzz[k] * ab_x + g_y_0_yyz_xxxzz[k];

                g_y_0_xyyz_xyyy[k] = -g_y_0_yyz_xyyy[k] * ab_x + g_y_0_yyz_xxyyy[k];

                g_y_0_xyyz_xyyz[k] = -g_y_0_yyz_xyyz[k] * ab_x + g_y_0_yyz_xxyyz[k];

                g_y_0_xyyz_xyzz[k] = -g_y_0_yyz_xyzz[k] * ab_x + g_y_0_yyz_xxyzz[k];

                g_y_0_xyyz_xzzz[k] = -g_y_0_yyz_xzzz[k] * ab_x + g_y_0_yyz_xxzzz[k];

                g_y_0_xyyz_yyyy[k] = -g_y_0_yyz_yyyy[k] * ab_x + g_y_0_yyz_xyyyy[k];

                g_y_0_xyyz_yyyz[k] = -g_y_0_yyz_yyyz[k] * ab_x + g_y_0_yyz_xyyyz[k];

                g_y_0_xyyz_yyzz[k] = -g_y_0_yyz_yyzz[k] * ab_x + g_y_0_yyz_xyyzz[k];

                g_y_0_xyyz_yzzz[k] = -g_y_0_yyz_yzzz[k] * ab_x + g_y_0_yyz_xyzzz[k];

                g_y_0_xyyz_zzzz[k] = -g_y_0_yyz_zzzz[k] * ab_x + g_y_0_yyz_xzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzz_xxxx, g_y_0_xyzz_xxxy, g_y_0_xyzz_xxxz, g_y_0_xyzz_xxyy, g_y_0_xyzz_xxyz, g_y_0_xyzz_xxzz, g_y_0_xyzz_xyyy, g_y_0_xyzz_xyyz, g_y_0_xyzz_xyzz, g_y_0_xyzz_xzzz, g_y_0_xyzz_yyyy, g_y_0_xyzz_yyyz, g_y_0_xyzz_yyzz, g_y_0_xyzz_yzzz, g_y_0_xyzz_zzzz, g_y_0_yzz_xxxx, g_y_0_yzz_xxxxx, g_y_0_yzz_xxxxy, g_y_0_yzz_xxxxz, g_y_0_yzz_xxxy, g_y_0_yzz_xxxyy, g_y_0_yzz_xxxyz, g_y_0_yzz_xxxz, g_y_0_yzz_xxxzz, g_y_0_yzz_xxyy, g_y_0_yzz_xxyyy, g_y_0_yzz_xxyyz, g_y_0_yzz_xxyz, g_y_0_yzz_xxyzz, g_y_0_yzz_xxzz, g_y_0_yzz_xxzzz, g_y_0_yzz_xyyy, g_y_0_yzz_xyyyy, g_y_0_yzz_xyyyz, g_y_0_yzz_xyyz, g_y_0_yzz_xyyzz, g_y_0_yzz_xyzz, g_y_0_yzz_xyzzz, g_y_0_yzz_xzzz, g_y_0_yzz_xzzzz, g_y_0_yzz_yyyy, g_y_0_yzz_yyyz, g_y_0_yzz_yyzz, g_y_0_yzz_yzzz, g_y_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzz_xxxx[k] = -g_y_0_yzz_xxxx[k] * ab_x + g_y_0_yzz_xxxxx[k];

                g_y_0_xyzz_xxxy[k] = -g_y_0_yzz_xxxy[k] * ab_x + g_y_0_yzz_xxxxy[k];

                g_y_0_xyzz_xxxz[k] = -g_y_0_yzz_xxxz[k] * ab_x + g_y_0_yzz_xxxxz[k];

                g_y_0_xyzz_xxyy[k] = -g_y_0_yzz_xxyy[k] * ab_x + g_y_0_yzz_xxxyy[k];

                g_y_0_xyzz_xxyz[k] = -g_y_0_yzz_xxyz[k] * ab_x + g_y_0_yzz_xxxyz[k];

                g_y_0_xyzz_xxzz[k] = -g_y_0_yzz_xxzz[k] * ab_x + g_y_0_yzz_xxxzz[k];

                g_y_0_xyzz_xyyy[k] = -g_y_0_yzz_xyyy[k] * ab_x + g_y_0_yzz_xxyyy[k];

                g_y_0_xyzz_xyyz[k] = -g_y_0_yzz_xyyz[k] * ab_x + g_y_0_yzz_xxyyz[k];

                g_y_0_xyzz_xyzz[k] = -g_y_0_yzz_xyzz[k] * ab_x + g_y_0_yzz_xxyzz[k];

                g_y_0_xyzz_xzzz[k] = -g_y_0_yzz_xzzz[k] * ab_x + g_y_0_yzz_xxzzz[k];

                g_y_0_xyzz_yyyy[k] = -g_y_0_yzz_yyyy[k] * ab_x + g_y_0_yzz_xyyyy[k];

                g_y_0_xyzz_yyyz[k] = -g_y_0_yzz_yyyz[k] * ab_x + g_y_0_yzz_xyyyz[k];

                g_y_0_xyzz_yyzz[k] = -g_y_0_yzz_yyzz[k] * ab_x + g_y_0_yzz_xyyzz[k];

                g_y_0_xyzz_yzzz[k] = -g_y_0_yzz_yzzz[k] * ab_x + g_y_0_yzz_xyzzz[k];

                g_y_0_xyzz_zzzz[k] = -g_y_0_yzz_zzzz[k] * ab_x + g_y_0_yzz_xzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzz_xxxx, g_y_0_xzzz_xxxy, g_y_0_xzzz_xxxz, g_y_0_xzzz_xxyy, g_y_0_xzzz_xxyz, g_y_0_xzzz_xxzz, g_y_0_xzzz_xyyy, g_y_0_xzzz_xyyz, g_y_0_xzzz_xyzz, g_y_0_xzzz_xzzz, g_y_0_xzzz_yyyy, g_y_0_xzzz_yyyz, g_y_0_xzzz_yyzz, g_y_0_xzzz_yzzz, g_y_0_xzzz_zzzz, g_y_0_zzz_xxxx, g_y_0_zzz_xxxxx, g_y_0_zzz_xxxxy, g_y_0_zzz_xxxxz, g_y_0_zzz_xxxy, g_y_0_zzz_xxxyy, g_y_0_zzz_xxxyz, g_y_0_zzz_xxxz, g_y_0_zzz_xxxzz, g_y_0_zzz_xxyy, g_y_0_zzz_xxyyy, g_y_0_zzz_xxyyz, g_y_0_zzz_xxyz, g_y_0_zzz_xxyzz, g_y_0_zzz_xxzz, g_y_0_zzz_xxzzz, g_y_0_zzz_xyyy, g_y_0_zzz_xyyyy, g_y_0_zzz_xyyyz, g_y_0_zzz_xyyz, g_y_0_zzz_xyyzz, g_y_0_zzz_xyzz, g_y_0_zzz_xyzzz, g_y_0_zzz_xzzz, g_y_0_zzz_xzzzz, g_y_0_zzz_yyyy, g_y_0_zzz_yyyz, g_y_0_zzz_yyzz, g_y_0_zzz_yzzz, g_y_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzz_xxxx[k] = -g_y_0_zzz_xxxx[k] * ab_x + g_y_0_zzz_xxxxx[k];

                g_y_0_xzzz_xxxy[k] = -g_y_0_zzz_xxxy[k] * ab_x + g_y_0_zzz_xxxxy[k];

                g_y_0_xzzz_xxxz[k] = -g_y_0_zzz_xxxz[k] * ab_x + g_y_0_zzz_xxxxz[k];

                g_y_0_xzzz_xxyy[k] = -g_y_0_zzz_xxyy[k] * ab_x + g_y_0_zzz_xxxyy[k];

                g_y_0_xzzz_xxyz[k] = -g_y_0_zzz_xxyz[k] * ab_x + g_y_0_zzz_xxxyz[k];

                g_y_0_xzzz_xxzz[k] = -g_y_0_zzz_xxzz[k] * ab_x + g_y_0_zzz_xxxzz[k];

                g_y_0_xzzz_xyyy[k] = -g_y_0_zzz_xyyy[k] * ab_x + g_y_0_zzz_xxyyy[k];

                g_y_0_xzzz_xyyz[k] = -g_y_0_zzz_xyyz[k] * ab_x + g_y_0_zzz_xxyyz[k];

                g_y_0_xzzz_xyzz[k] = -g_y_0_zzz_xyzz[k] * ab_x + g_y_0_zzz_xxyzz[k];

                g_y_0_xzzz_xzzz[k] = -g_y_0_zzz_xzzz[k] * ab_x + g_y_0_zzz_xxzzz[k];

                g_y_0_xzzz_yyyy[k] = -g_y_0_zzz_yyyy[k] * ab_x + g_y_0_zzz_xyyyy[k];

                g_y_0_xzzz_yyyz[k] = -g_y_0_zzz_yyyz[k] * ab_x + g_y_0_zzz_xyyyz[k];

                g_y_0_xzzz_yyzz[k] = -g_y_0_zzz_yyzz[k] * ab_x + g_y_0_zzz_xyyzz[k];

                g_y_0_xzzz_yzzz[k] = -g_y_0_zzz_yzzz[k] * ab_x + g_y_0_zzz_xyzzz[k];

                g_y_0_xzzz_zzzz[k] = -g_y_0_zzz_zzzz[k] * ab_x + g_y_0_zzz_xzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_xxxx, g_y_0_yyy_xxxxy, g_y_0_yyy_xxxy, g_y_0_yyy_xxxyy, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxz, g_y_0_yyy_xxyy, g_y_0_yyy_xxyyy, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxzz, g_y_0_yyy_xyyy, g_y_0_yyy_xyyyy, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xzzz, g_y_0_yyy_yyyy, g_y_0_yyy_yyyyy, g_y_0_yyy_yyyyz, g_y_0_yyy_yyyz, g_y_0_yyy_yyyzz, g_y_0_yyy_yyzz, g_y_0_yyy_yyzzz, g_y_0_yyy_yzzz, g_y_0_yyy_yzzzz, g_y_0_yyy_zzzz, g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzzz, g_yyy_xxxx, g_yyy_xxxy, g_yyy_xxxz, g_yyy_xxyy, g_yyy_xxyz, g_yyy_xxzz, g_yyy_xyyy, g_yyy_xyyz, g_yyy_xyzz, g_yyy_xzzz, g_yyy_yyyy, g_yyy_yyyz, g_yyy_yyzz, g_yyy_yzzz, g_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyy_xxxx[k] = -g_yyy_xxxx[k] - g_y_0_yyy_xxxx[k] * ab_y + g_y_0_yyy_xxxxy[k];

                g_y_0_yyyy_xxxy[k] = -g_yyy_xxxy[k] - g_y_0_yyy_xxxy[k] * ab_y + g_y_0_yyy_xxxyy[k];

                g_y_0_yyyy_xxxz[k] = -g_yyy_xxxz[k] - g_y_0_yyy_xxxz[k] * ab_y + g_y_0_yyy_xxxyz[k];

                g_y_0_yyyy_xxyy[k] = -g_yyy_xxyy[k] - g_y_0_yyy_xxyy[k] * ab_y + g_y_0_yyy_xxyyy[k];

                g_y_0_yyyy_xxyz[k] = -g_yyy_xxyz[k] - g_y_0_yyy_xxyz[k] * ab_y + g_y_0_yyy_xxyyz[k];

                g_y_0_yyyy_xxzz[k] = -g_yyy_xxzz[k] - g_y_0_yyy_xxzz[k] * ab_y + g_y_0_yyy_xxyzz[k];

                g_y_0_yyyy_xyyy[k] = -g_yyy_xyyy[k] - g_y_0_yyy_xyyy[k] * ab_y + g_y_0_yyy_xyyyy[k];

                g_y_0_yyyy_xyyz[k] = -g_yyy_xyyz[k] - g_y_0_yyy_xyyz[k] * ab_y + g_y_0_yyy_xyyyz[k];

                g_y_0_yyyy_xyzz[k] = -g_yyy_xyzz[k] - g_y_0_yyy_xyzz[k] * ab_y + g_y_0_yyy_xyyzz[k];

                g_y_0_yyyy_xzzz[k] = -g_yyy_xzzz[k] - g_y_0_yyy_xzzz[k] * ab_y + g_y_0_yyy_xyzzz[k];

                g_y_0_yyyy_yyyy[k] = -g_yyy_yyyy[k] - g_y_0_yyy_yyyy[k] * ab_y + g_y_0_yyy_yyyyy[k];

                g_y_0_yyyy_yyyz[k] = -g_yyy_yyyz[k] - g_y_0_yyy_yyyz[k] * ab_y + g_y_0_yyy_yyyyz[k];

                g_y_0_yyyy_yyzz[k] = -g_yyy_yyzz[k] - g_y_0_yyy_yyzz[k] * ab_y + g_y_0_yyy_yyyzz[k];

                g_y_0_yyyy_yzzz[k] = -g_yyy_yzzz[k] - g_y_0_yyy_yzzz[k] * ab_y + g_y_0_yyy_yyzzz[k];

                g_y_0_yyyy_zzzz[k] = -g_yyy_zzzz[k] - g_y_0_yyy_zzzz[k] * ab_y + g_y_0_yyy_yzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_xxxx, g_y_0_yyy_xxxxz, g_y_0_yyy_xxxy, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxz, g_y_0_yyy_xxxzz, g_y_0_yyy_xxyy, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxzz, g_y_0_yyy_xxzzz, g_y_0_yyy_xyyy, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xzzz, g_y_0_yyy_xzzzz, g_y_0_yyy_yyyy, g_y_0_yyy_yyyyz, g_y_0_yyy_yyyz, g_y_0_yyy_yyyzz, g_y_0_yyy_yyzz, g_y_0_yyy_yyzzz, g_y_0_yyy_yzzz, g_y_0_yyy_yzzzz, g_y_0_yyy_zzzz, g_y_0_yyy_zzzzz, g_y_0_yyyz_xxxx, g_y_0_yyyz_xxxy, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxyy, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xyyy, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_yyyy, g_y_0_yyyz_yyyz, g_y_0_yyyz_yyzz, g_y_0_yyyz_yzzz, g_y_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyz_xxxx[k] = -g_y_0_yyy_xxxx[k] * ab_z + g_y_0_yyy_xxxxz[k];

                g_y_0_yyyz_xxxy[k] = -g_y_0_yyy_xxxy[k] * ab_z + g_y_0_yyy_xxxyz[k];

                g_y_0_yyyz_xxxz[k] = -g_y_0_yyy_xxxz[k] * ab_z + g_y_0_yyy_xxxzz[k];

                g_y_0_yyyz_xxyy[k] = -g_y_0_yyy_xxyy[k] * ab_z + g_y_0_yyy_xxyyz[k];

                g_y_0_yyyz_xxyz[k] = -g_y_0_yyy_xxyz[k] * ab_z + g_y_0_yyy_xxyzz[k];

                g_y_0_yyyz_xxzz[k] = -g_y_0_yyy_xxzz[k] * ab_z + g_y_0_yyy_xxzzz[k];

                g_y_0_yyyz_xyyy[k] = -g_y_0_yyy_xyyy[k] * ab_z + g_y_0_yyy_xyyyz[k];

                g_y_0_yyyz_xyyz[k] = -g_y_0_yyy_xyyz[k] * ab_z + g_y_0_yyy_xyyzz[k];

                g_y_0_yyyz_xyzz[k] = -g_y_0_yyy_xyzz[k] * ab_z + g_y_0_yyy_xyzzz[k];

                g_y_0_yyyz_xzzz[k] = -g_y_0_yyy_xzzz[k] * ab_z + g_y_0_yyy_xzzzz[k];

                g_y_0_yyyz_yyyy[k] = -g_y_0_yyy_yyyy[k] * ab_z + g_y_0_yyy_yyyyz[k];

                g_y_0_yyyz_yyyz[k] = -g_y_0_yyy_yyyz[k] * ab_z + g_y_0_yyy_yyyzz[k];

                g_y_0_yyyz_yyzz[k] = -g_y_0_yyy_yyzz[k] * ab_z + g_y_0_yyy_yyzzz[k];

                g_y_0_yyyz_yzzz[k] = -g_y_0_yyy_yzzz[k] * ab_z + g_y_0_yyy_yzzzz[k];

                g_y_0_yyyz_zzzz[k] = -g_y_0_yyy_zzzz[k] * ab_z + g_y_0_yyy_zzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyz_xxxx, g_y_0_yyz_xxxxz, g_y_0_yyz_xxxy, g_y_0_yyz_xxxyz, g_y_0_yyz_xxxz, g_y_0_yyz_xxxzz, g_y_0_yyz_xxyy, g_y_0_yyz_xxyyz, g_y_0_yyz_xxyz, g_y_0_yyz_xxyzz, g_y_0_yyz_xxzz, g_y_0_yyz_xxzzz, g_y_0_yyz_xyyy, g_y_0_yyz_xyyyz, g_y_0_yyz_xyyz, g_y_0_yyz_xyyzz, g_y_0_yyz_xyzz, g_y_0_yyz_xyzzz, g_y_0_yyz_xzzz, g_y_0_yyz_xzzzz, g_y_0_yyz_yyyy, g_y_0_yyz_yyyyz, g_y_0_yyz_yyyz, g_y_0_yyz_yyyzz, g_y_0_yyz_yyzz, g_y_0_yyz_yyzzz, g_y_0_yyz_yzzz, g_y_0_yyz_yzzzz, g_y_0_yyz_zzzz, g_y_0_yyz_zzzzz, g_y_0_yyzz_xxxx, g_y_0_yyzz_xxxy, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxyy, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xyyy, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_yyyy, g_y_0_yyzz_yyyz, g_y_0_yyzz_yyzz, g_y_0_yyzz_yzzz, g_y_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzz_xxxx[k] = -g_y_0_yyz_xxxx[k] * ab_z + g_y_0_yyz_xxxxz[k];

                g_y_0_yyzz_xxxy[k] = -g_y_0_yyz_xxxy[k] * ab_z + g_y_0_yyz_xxxyz[k];

                g_y_0_yyzz_xxxz[k] = -g_y_0_yyz_xxxz[k] * ab_z + g_y_0_yyz_xxxzz[k];

                g_y_0_yyzz_xxyy[k] = -g_y_0_yyz_xxyy[k] * ab_z + g_y_0_yyz_xxyyz[k];

                g_y_0_yyzz_xxyz[k] = -g_y_0_yyz_xxyz[k] * ab_z + g_y_0_yyz_xxyzz[k];

                g_y_0_yyzz_xxzz[k] = -g_y_0_yyz_xxzz[k] * ab_z + g_y_0_yyz_xxzzz[k];

                g_y_0_yyzz_xyyy[k] = -g_y_0_yyz_xyyy[k] * ab_z + g_y_0_yyz_xyyyz[k];

                g_y_0_yyzz_xyyz[k] = -g_y_0_yyz_xyyz[k] * ab_z + g_y_0_yyz_xyyzz[k];

                g_y_0_yyzz_xyzz[k] = -g_y_0_yyz_xyzz[k] * ab_z + g_y_0_yyz_xyzzz[k];

                g_y_0_yyzz_xzzz[k] = -g_y_0_yyz_xzzz[k] * ab_z + g_y_0_yyz_xzzzz[k];

                g_y_0_yyzz_yyyy[k] = -g_y_0_yyz_yyyy[k] * ab_z + g_y_0_yyz_yyyyz[k];

                g_y_0_yyzz_yyyz[k] = -g_y_0_yyz_yyyz[k] * ab_z + g_y_0_yyz_yyyzz[k];

                g_y_0_yyzz_yyzz[k] = -g_y_0_yyz_yyzz[k] * ab_z + g_y_0_yyz_yyzzz[k];

                g_y_0_yyzz_yzzz[k] = -g_y_0_yyz_yzzz[k] * ab_z + g_y_0_yyz_yzzzz[k];

                g_y_0_yyzz_zzzz[k] = -g_y_0_yyz_zzzz[k] * ab_z + g_y_0_yyz_zzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzz_xxxx, g_y_0_yzz_xxxxz, g_y_0_yzz_xxxy, g_y_0_yzz_xxxyz, g_y_0_yzz_xxxz, g_y_0_yzz_xxxzz, g_y_0_yzz_xxyy, g_y_0_yzz_xxyyz, g_y_0_yzz_xxyz, g_y_0_yzz_xxyzz, g_y_0_yzz_xxzz, g_y_0_yzz_xxzzz, g_y_0_yzz_xyyy, g_y_0_yzz_xyyyz, g_y_0_yzz_xyyz, g_y_0_yzz_xyyzz, g_y_0_yzz_xyzz, g_y_0_yzz_xyzzz, g_y_0_yzz_xzzz, g_y_0_yzz_xzzzz, g_y_0_yzz_yyyy, g_y_0_yzz_yyyyz, g_y_0_yzz_yyyz, g_y_0_yzz_yyyzz, g_y_0_yzz_yyzz, g_y_0_yzz_yyzzz, g_y_0_yzz_yzzz, g_y_0_yzz_yzzzz, g_y_0_yzz_zzzz, g_y_0_yzz_zzzzz, g_y_0_yzzz_xxxx, g_y_0_yzzz_xxxy, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxyy, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xyyy, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_yyyy, g_y_0_yzzz_yyyz, g_y_0_yzzz_yyzz, g_y_0_yzzz_yzzz, g_y_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzz_xxxx[k] = -g_y_0_yzz_xxxx[k] * ab_z + g_y_0_yzz_xxxxz[k];

                g_y_0_yzzz_xxxy[k] = -g_y_0_yzz_xxxy[k] * ab_z + g_y_0_yzz_xxxyz[k];

                g_y_0_yzzz_xxxz[k] = -g_y_0_yzz_xxxz[k] * ab_z + g_y_0_yzz_xxxzz[k];

                g_y_0_yzzz_xxyy[k] = -g_y_0_yzz_xxyy[k] * ab_z + g_y_0_yzz_xxyyz[k];

                g_y_0_yzzz_xxyz[k] = -g_y_0_yzz_xxyz[k] * ab_z + g_y_0_yzz_xxyzz[k];

                g_y_0_yzzz_xxzz[k] = -g_y_0_yzz_xxzz[k] * ab_z + g_y_0_yzz_xxzzz[k];

                g_y_0_yzzz_xyyy[k] = -g_y_0_yzz_xyyy[k] * ab_z + g_y_0_yzz_xyyyz[k];

                g_y_0_yzzz_xyyz[k] = -g_y_0_yzz_xyyz[k] * ab_z + g_y_0_yzz_xyyzz[k];

                g_y_0_yzzz_xyzz[k] = -g_y_0_yzz_xyzz[k] * ab_z + g_y_0_yzz_xyzzz[k];

                g_y_0_yzzz_xzzz[k] = -g_y_0_yzz_xzzz[k] * ab_z + g_y_0_yzz_xzzzz[k];

                g_y_0_yzzz_yyyy[k] = -g_y_0_yzz_yyyy[k] * ab_z + g_y_0_yzz_yyyyz[k];

                g_y_0_yzzz_yyyz[k] = -g_y_0_yzz_yyyz[k] * ab_z + g_y_0_yzz_yyyzz[k];

                g_y_0_yzzz_yyzz[k] = -g_y_0_yzz_yyzz[k] * ab_z + g_y_0_yzz_yyzzz[k];

                g_y_0_yzzz_yzzz[k] = -g_y_0_yzz_yzzz[k] * ab_z + g_y_0_yzz_yzzzz[k];

                g_y_0_yzzz_zzzz[k] = -g_y_0_yzz_zzzz[k] * ab_z + g_y_0_yzz_zzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzz_xxxx, g_y_0_zzz_xxxxz, g_y_0_zzz_xxxy, g_y_0_zzz_xxxyz, g_y_0_zzz_xxxz, g_y_0_zzz_xxxzz, g_y_0_zzz_xxyy, g_y_0_zzz_xxyyz, g_y_0_zzz_xxyz, g_y_0_zzz_xxyzz, g_y_0_zzz_xxzz, g_y_0_zzz_xxzzz, g_y_0_zzz_xyyy, g_y_0_zzz_xyyyz, g_y_0_zzz_xyyz, g_y_0_zzz_xyyzz, g_y_0_zzz_xyzz, g_y_0_zzz_xyzzz, g_y_0_zzz_xzzz, g_y_0_zzz_xzzzz, g_y_0_zzz_yyyy, g_y_0_zzz_yyyyz, g_y_0_zzz_yyyz, g_y_0_zzz_yyyzz, g_y_0_zzz_yyzz, g_y_0_zzz_yyzzz, g_y_0_zzz_yzzz, g_y_0_zzz_yzzzz, g_y_0_zzz_zzzz, g_y_0_zzz_zzzzz, g_y_0_zzzz_xxxx, g_y_0_zzzz_xxxy, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxyy, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xyyy, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_yyyy, g_y_0_zzzz_yyyz, g_y_0_zzzz_yyzz, g_y_0_zzzz_yzzz, g_y_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzz_xxxx[k] = -g_y_0_zzz_xxxx[k] * ab_z + g_y_0_zzz_xxxxz[k];

                g_y_0_zzzz_xxxy[k] = -g_y_0_zzz_xxxy[k] * ab_z + g_y_0_zzz_xxxyz[k];

                g_y_0_zzzz_xxxz[k] = -g_y_0_zzz_xxxz[k] * ab_z + g_y_0_zzz_xxxzz[k];

                g_y_0_zzzz_xxyy[k] = -g_y_0_zzz_xxyy[k] * ab_z + g_y_0_zzz_xxyyz[k];

                g_y_0_zzzz_xxyz[k] = -g_y_0_zzz_xxyz[k] * ab_z + g_y_0_zzz_xxyzz[k];

                g_y_0_zzzz_xxzz[k] = -g_y_0_zzz_xxzz[k] * ab_z + g_y_0_zzz_xxzzz[k];

                g_y_0_zzzz_xyyy[k] = -g_y_0_zzz_xyyy[k] * ab_z + g_y_0_zzz_xyyyz[k];

                g_y_0_zzzz_xyyz[k] = -g_y_0_zzz_xyyz[k] * ab_z + g_y_0_zzz_xyyzz[k];

                g_y_0_zzzz_xyzz[k] = -g_y_0_zzz_xyzz[k] * ab_z + g_y_0_zzz_xyzzz[k];

                g_y_0_zzzz_xzzz[k] = -g_y_0_zzz_xzzz[k] * ab_z + g_y_0_zzz_xzzzz[k];

                g_y_0_zzzz_yyyy[k] = -g_y_0_zzz_yyyy[k] * ab_z + g_y_0_zzz_yyyyz[k];

                g_y_0_zzzz_yyyz[k] = -g_y_0_zzz_yyyz[k] * ab_z + g_y_0_zzz_yyyzz[k];

                g_y_0_zzzz_yyzz[k] = -g_y_0_zzz_yyzz[k] * ab_z + g_y_0_zzz_yyzzz[k];

                g_y_0_zzzz_yzzz[k] = -g_y_0_zzz_yzzz[k] * ab_z + g_y_0_zzz_yzzzz[k];

                g_y_0_zzzz_zzzz[k] = -g_y_0_zzz_zzzz[k] * ab_z + g_y_0_zzz_zzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 450 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 451 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 452 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 453 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 454 * ccomps * dcomps);

            auto g_z_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 455 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 456 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 457 * ccomps * dcomps);

            auto g_z_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 458 * ccomps * dcomps);

            auto g_z_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 459 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 460 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 461 * ccomps * dcomps);

            auto g_z_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 462 * ccomps * dcomps);

            auto g_z_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 463 * ccomps * dcomps);

            auto g_z_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxx_xxxx, g_z_0_xxx_xxxxx, g_z_0_xxx_xxxxy, g_z_0_xxx_xxxxz, g_z_0_xxx_xxxy, g_z_0_xxx_xxxyy, g_z_0_xxx_xxxyz, g_z_0_xxx_xxxz, g_z_0_xxx_xxxzz, g_z_0_xxx_xxyy, g_z_0_xxx_xxyyy, g_z_0_xxx_xxyyz, g_z_0_xxx_xxyz, g_z_0_xxx_xxyzz, g_z_0_xxx_xxzz, g_z_0_xxx_xxzzz, g_z_0_xxx_xyyy, g_z_0_xxx_xyyyy, g_z_0_xxx_xyyyz, g_z_0_xxx_xyyz, g_z_0_xxx_xyyzz, g_z_0_xxx_xyzz, g_z_0_xxx_xyzzz, g_z_0_xxx_xzzz, g_z_0_xxx_xzzzz, g_z_0_xxx_yyyy, g_z_0_xxx_yyyz, g_z_0_xxx_yyzz, g_z_0_xxx_yzzz, g_z_0_xxx_zzzz, g_z_0_xxxx_xxxx, g_z_0_xxxx_xxxy, g_z_0_xxxx_xxxz, g_z_0_xxxx_xxyy, g_z_0_xxxx_xxyz, g_z_0_xxxx_xxzz, g_z_0_xxxx_xyyy, g_z_0_xxxx_xyyz, g_z_0_xxxx_xyzz, g_z_0_xxxx_xzzz, g_z_0_xxxx_yyyy, g_z_0_xxxx_yyyz, g_z_0_xxxx_yyzz, g_z_0_xxxx_yzzz, g_z_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxx_xxxx[k] = -g_z_0_xxx_xxxx[k] * ab_x + g_z_0_xxx_xxxxx[k];

                g_z_0_xxxx_xxxy[k] = -g_z_0_xxx_xxxy[k] * ab_x + g_z_0_xxx_xxxxy[k];

                g_z_0_xxxx_xxxz[k] = -g_z_0_xxx_xxxz[k] * ab_x + g_z_0_xxx_xxxxz[k];

                g_z_0_xxxx_xxyy[k] = -g_z_0_xxx_xxyy[k] * ab_x + g_z_0_xxx_xxxyy[k];

                g_z_0_xxxx_xxyz[k] = -g_z_0_xxx_xxyz[k] * ab_x + g_z_0_xxx_xxxyz[k];

                g_z_0_xxxx_xxzz[k] = -g_z_0_xxx_xxzz[k] * ab_x + g_z_0_xxx_xxxzz[k];

                g_z_0_xxxx_xyyy[k] = -g_z_0_xxx_xyyy[k] * ab_x + g_z_0_xxx_xxyyy[k];

                g_z_0_xxxx_xyyz[k] = -g_z_0_xxx_xyyz[k] * ab_x + g_z_0_xxx_xxyyz[k];

                g_z_0_xxxx_xyzz[k] = -g_z_0_xxx_xyzz[k] * ab_x + g_z_0_xxx_xxyzz[k];

                g_z_0_xxxx_xzzz[k] = -g_z_0_xxx_xzzz[k] * ab_x + g_z_0_xxx_xxzzz[k];

                g_z_0_xxxx_yyyy[k] = -g_z_0_xxx_yyyy[k] * ab_x + g_z_0_xxx_xyyyy[k];

                g_z_0_xxxx_yyyz[k] = -g_z_0_xxx_yyyz[k] * ab_x + g_z_0_xxx_xyyyz[k];

                g_z_0_xxxx_yyzz[k] = -g_z_0_xxx_yyzz[k] * ab_x + g_z_0_xxx_xyyzz[k];

                g_z_0_xxxx_yzzz[k] = -g_z_0_xxx_yzzz[k] * ab_x + g_z_0_xxx_xyzzz[k];

                g_z_0_xxxx_zzzz[k] = -g_z_0_xxx_zzzz[k] * ab_x + g_z_0_xxx_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 465 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 466 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 467 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 468 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 469 * ccomps * dcomps);

            auto g_z_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 470 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 471 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 472 * ccomps * dcomps);

            auto g_z_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 473 * ccomps * dcomps);

            auto g_z_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 474 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 475 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 476 * ccomps * dcomps);

            auto g_z_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 477 * ccomps * dcomps);

            auto g_z_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 478 * ccomps * dcomps);

            auto g_z_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxy_xxxx, g_z_0_xxxy_xxxy, g_z_0_xxxy_xxxz, g_z_0_xxxy_xxyy, g_z_0_xxxy_xxyz, g_z_0_xxxy_xxzz, g_z_0_xxxy_xyyy, g_z_0_xxxy_xyyz, g_z_0_xxxy_xyzz, g_z_0_xxxy_xzzz, g_z_0_xxxy_yyyy, g_z_0_xxxy_yyyz, g_z_0_xxxy_yyzz, g_z_0_xxxy_yzzz, g_z_0_xxxy_zzzz, g_z_0_xxy_xxxx, g_z_0_xxy_xxxxx, g_z_0_xxy_xxxxy, g_z_0_xxy_xxxxz, g_z_0_xxy_xxxy, g_z_0_xxy_xxxyy, g_z_0_xxy_xxxyz, g_z_0_xxy_xxxz, g_z_0_xxy_xxxzz, g_z_0_xxy_xxyy, g_z_0_xxy_xxyyy, g_z_0_xxy_xxyyz, g_z_0_xxy_xxyz, g_z_0_xxy_xxyzz, g_z_0_xxy_xxzz, g_z_0_xxy_xxzzz, g_z_0_xxy_xyyy, g_z_0_xxy_xyyyy, g_z_0_xxy_xyyyz, g_z_0_xxy_xyyz, g_z_0_xxy_xyyzz, g_z_0_xxy_xyzz, g_z_0_xxy_xyzzz, g_z_0_xxy_xzzz, g_z_0_xxy_xzzzz, g_z_0_xxy_yyyy, g_z_0_xxy_yyyz, g_z_0_xxy_yyzz, g_z_0_xxy_yzzz, g_z_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxy_xxxx[k] = -g_z_0_xxy_xxxx[k] * ab_x + g_z_0_xxy_xxxxx[k];

                g_z_0_xxxy_xxxy[k] = -g_z_0_xxy_xxxy[k] * ab_x + g_z_0_xxy_xxxxy[k];

                g_z_0_xxxy_xxxz[k] = -g_z_0_xxy_xxxz[k] * ab_x + g_z_0_xxy_xxxxz[k];

                g_z_0_xxxy_xxyy[k] = -g_z_0_xxy_xxyy[k] * ab_x + g_z_0_xxy_xxxyy[k];

                g_z_0_xxxy_xxyz[k] = -g_z_0_xxy_xxyz[k] * ab_x + g_z_0_xxy_xxxyz[k];

                g_z_0_xxxy_xxzz[k] = -g_z_0_xxy_xxzz[k] * ab_x + g_z_0_xxy_xxxzz[k];

                g_z_0_xxxy_xyyy[k] = -g_z_0_xxy_xyyy[k] * ab_x + g_z_0_xxy_xxyyy[k];

                g_z_0_xxxy_xyyz[k] = -g_z_0_xxy_xyyz[k] * ab_x + g_z_0_xxy_xxyyz[k];

                g_z_0_xxxy_xyzz[k] = -g_z_0_xxy_xyzz[k] * ab_x + g_z_0_xxy_xxyzz[k];

                g_z_0_xxxy_xzzz[k] = -g_z_0_xxy_xzzz[k] * ab_x + g_z_0_xxy_xxzzz[k];

                g_z_0_xxxy_yyyy[k] = -g_z_0_xxy_yyyy[k] * ab_x + g_z_0_xxy_xyyyy[k];

                g_z_0_xxxy_yyyz[k] = -g_z_0_xxy_yyyz[k] * ab_x + g_z_0_xxy_xyyyz[k];

                g_z_0_xxxy_yyzz[k] = -g_z_0_xxy_yyzz[k] * ab_x + g_z_0_xxy_xyyzz[k];

                g_z_0_xxxy_yzzz[k] = -g_z_0_xxy_yzzz[k] * ab_x + g_z_0_xxy_xyzzz[k];

                g_z_0_xxxy_zzzz[k] = -g_z_0_xxy_zzzz[k] * ab_x + g_z_0_xxy_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 480 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 481 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 482 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 483 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 484 * ccomps * dcomps);

            auto g_z_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 485 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 486 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 487 * ccomps * dcomps);

            auto g_z_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 488 * ccomps * dcomps);

            auto g_z_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 489 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 490 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 491 * ccomps * dcomps);

            auto g_z_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 492 * ccomps * dcomps);

            auto g_z_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 493 * ccomps * dcomps);

            auto g_z_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxz_xxxx, g_z_0_xxxz_xxxy, g_z_0_xxxz_xxxz, g_z_0_xxxz_xxyy, g_z_0_xxxz_xxyz, g_z_0_xxxz_xxzz, g_z_0_xxxz_xyyy, g_z_0_xxxz_xyyz, g_z_0_xxxz_xyzz, g_z_0_xxxz_xzzz, g_z_0_xxxz_yyyy, g_z_0_xxxz_yyyz, g_z_0_xxxz_yyzz, g_z_0_xxxz_yzzz, g_z_0_xxxz_zzzz, g_z_0_xxz_xxxx, g_z_0_xxz_xxxxx, g_z_0_xxz_xxxxy, g_z_0_xxz_xxxxz, g_z_0_xxz_xxxy, g_z_0_xxz_xxxyy, g_z_0_xxz_xxxyz, g_z_0_xxz_xxxz, g_z_0_xxz_xxxzz, g_z_0_xxz_xxyy, g_z_0_xxz_xxyyy, g_z_0_xxz_xxyyz, g_z_0_xxz_xxyz, g_z_0_xxz_xxyzz, g_z_0_xxz_xxzz, g_z_0_xxz_xxzzz, g_z_0_xxz_xyyy, g_z_0_xxz_xyyyy, g_z_0_xxz_xyyyz, g_z_0_xxz_xyyz, g_z_0_xxz_xyyzz, g_z_0_xxz_xyzz, g_z_0_xxz_xyzzz, g_z_0_xxz_xzzz, g_z_0_xxz_xzzzz, g_z_0_xxz_yyyy, g_z_0_xxz_yyyz, g_z_0_xxz_yyzz, g_z_0_xxz_yzzz, g_z_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxz_xxxx[k] = -g_z_0_xxz_xxxx[k] * ab_x + g_z_0_xxz_xxxxx[k];

                g_z_0_xxxz_xxxy[k] = -g_z_0_xxz_xxxy[k] * ab_x + g_z_0_xxz_xxxxy[k];

                g_z_0_xxxz_xxxz[k] = -g_z_0_xxz_xxxz[k] * ab_x + g_z_0_xxz_xxxxz[k];

                g_z_0_xxxz_xxyy[k] = -g_z_0_xxz_xxyy[k] * ab_x + g_z_0_xxz_xxxyy[k];

                g_z_0_xxxz_xxyz[k] = -g_z_0_xxz_xxyz[k] * ab_x + g_z_0_xxz_xxxyz[k];

                g_z_0_xxxz_xxzz[k] = -g_z_0_xxz_xxzz[k] * ab_x + g_z_0_xxz_xxxzz[k];

                g_z_0_xxxz_xyyy[k] = -g_z_0_xxz_xyyy[k] * ab_x + g_z_0_xxz_xxyyy[k];

                g_z_0_xxxz_xyyz[k] = -g_z_0_xxz_xyyz[k] * ab_x + g_z_0_xxz_xxyyz[k];

                g_z_0_xxxz_xyzz[k] = -g_z_0_xxz_xyzz[k] * ab_x + g_z_0_xxz_xxyzz[k];

                g_z_0_xxxz_xzzz[k] = -g_z_0_xxz_xzzz[k] * ab_x + g_z_0_xxz_xxzzz[k];

                g_z_0_xxxz_yyyy[k] = -g_z_0_xxz_yyyy[k] * ab_x + g_z_0_xxz_xyyyy[k];

                g_z_0_xxxz_yyyz[k] = -g_z_0_xxz_yyyz[k] * ab_x + g_z_0_xxz_xyyyz[k];

                g_z_0_xxxz_yyzz[k] = -g_z_0_xxz_yyzz[k] * ab_x + g_z_0_xxz_xyyzz[k];

                g_z_0_xxxz_yzzz[k] = -g_z_0_xxz_yzzz[k] * ab_x + g_z_0_xxz_xyzzz[k];

                g_z_0_xxxz_zzzz[k] = -g_z_0_xxz_zzzz[k] * ab_x + g_z_0_xxz_xzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 495 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 496 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 497 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 498 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 499 * ccomps * dcomps);

            auto g_z_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 500 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 501 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 502 * ccomps * dcomps);

            auto g_z_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 503 * ccomps * dcomps);

            auto g_z_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 504 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 505 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 506 * ccomps * dcomps);

            auto g_z_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 507 * ccomps * dcomps);

            auto g_z_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 508 * ccomps * dcomps);

            auto g_z_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyy_xxxx, g_z_0_xxyy_xxxy, g_z_0_xxyy_xxxz, g_z_0_xxyy_xxyy, g_z_0_xxyy_xxyz, g_z_0_xxyy_xxzz, g_z_0_xxyy_xyyy, g_z_0_xxyy_xyyz, g_z_0_xxyy_xyzz, g_z_0_xxyy_xzzz, g_z_0_xxyy_yyyy, g_z_0_xxyy_yyyz, g_z_0_xxyy_yyzz, g_z_0_xxyy_yzzz, g_z_0_xxyy_zzzz, g_z_0_xyy_xxxx, g_z_0_xyy_xxxxx, g_z_0_xyy_xxxxy, g_z_0_xyy_xxxxz, g_z_0_xyy_xxxy, g_z_0_xyy_xxxyy, g_z_0_xyy_xxxyz, g_z_0_xyy_xxxz, g_z_0_xyy_xxxzz, g_z_0_xyy_xxyy, g_z_0_xyy_xxyyy, g_z_0_xyy_xxyyz, g_z_0_xyy_xxyz, g_z_0_xyy_xxyzz, g_z_0_xyy_xxzz, g_z_0_xyy_xxzzz, g_z_0_xyy_xyyy, g_z_0_xyy_xyyyy, g_z_0_xyy_xyyyz, g_z_0_xyy_xyyz, g_z_0_xyy_xyyzz, g_z_0_xyy_xyzz, g_z_0_xyy_xyzzz, g_z_0_xyy_xzzz, g_z_0_xyy_xzzzz, g_z_0_xyy_yyyy, g_z_0_xyy_yyyz, g_z_0_xyy_yyzz, g_z_0_xyy_yzzz, g_z_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyy_xxxx[k] = -g_z_0_xyy_xxxx[k] * ab_x + g_z_0_xyy_xxxxx[k];

                g_z_0_xxyy_xxxy[k] = -g_z_0_xyy_xxxy[k] * ab_x + g_z_0_xyy_xxxxy[k];

                g_z_0_xxyy_xxxz[k] = -g_z_0_xyy_xxxz[k] * ab_x + g_z_0_xyy_xxxxz[k];

                g_z_0_xxyy_xxyy[k] = -g_z_0_xyy_xxyy[k] * ab_x + g_z_0_xyy_xxxyy[k];

                g_z_0_xxyy_xxyz[k] = -g_z_0_xyy_xxyz[k] * ab_x + g_z_0_xyy_xxxyz[k];

                g_z_0_xxyy_xxzz[k] = -g_z_0_xyy_xxzz[k] * ab_x + g_z_0_xyy_xxxzz[k];

                g_z_0_xxyy_xyyy[k] = -g_z_0_xyy_xyyy[k] * ab_x + g_z_0_xyy_xxyyy[k];

                g_z_0_xxyy_xyyz[k] = -g_z_0_xyy_xyyz[k] * ab_x + g_z_0_xyy_xxyyz[k];

                g_z_0_xxyy_xyzz[k] = -g_z_0_xyy_xyzz[k] * ab_x + g_z_0_xyy_xxyzz[k];

                g_z_0_xxyy_xzzz[k] = -g_z_0_xyy_xzzz[k] * ab_x + g_z_0_xyy_xxzzz[k];

                g_z_0_xxyy_yyyy[k] = -g_z_0_xyy_yyyy[k] * ab_x + g_z_0_xyy_xyyyy[k];

                g_z_0_xxyy_yyyz[k] = -g_z_0_xyy_yyyz[k] * ab_x + g_z_0_xyy_xyyyz[k];

                g_z_0_xxyy_yyzz[k] = -g_z_0_xyy_yyzz[k] * ab_x + g_z_0_xyy_xyyzz[k];

                g_z_0_xxyy_yzzz[k] = -g_z_0_xyy_yzzz[k] * ab_x + g_z_0_xyy_xyzzz[k];

                g_z_0_xxyy_zzzz[k] = -g_z_0_xyy_zzzz[k] * ab_x + g_z_0_xyy_xzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 510 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 511 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 512 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 513 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 514 * ccomps * dcomps);

            auto g_z_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 515 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 516 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 517 * ccomps * dcomps);

            auto g_z_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 518 * ccomps * dcomps);

            auto g_z_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 519 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 520 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 521 * ccomps * dcomps);

            auto g_z_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 522 * ccomps * dcomps);

            auto g_z_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 523 * ccomps * dcomps);

            auto g_z_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyz_xxxx, g_z_0_xxyz_xxxy, g_z_0_xxyz_xxxz, g_z_0_xxyz_xxyy, g_z_0_xxyz_xxyz, g_z_0_xxyz_xxzz, g_z_0_xxyz_xyyy, g_z_0_xxyz_xyyz, g_z_0_xxyz_xyzz, g_z_0_xxyz_xzzz, g_z_0_xxyz_yyyy, g_z_0_xxyz_yyyz, g_z_0_xxyz_yyzz, g_z_0_xxyz_yzzz, g_z_0_xxyz_zzzz, g_z_0_xyz_xxxx, g_z_0_xyz_xxxxx, g_z_0_xyz_xxxxy, g_z_0_xyz_xxxxz, g_z_0_xyz_xxxy, g_z_0_xyz_xxxyy, g_z_0_xyz_xxxyz, g_z_0_xyz_xxxz, g_z_0_xyz_xxxzz, g_z_0_xyz_xxyy, g_z_0_xyz_xxyyy, g_z_0_xyz_xxyyz, g_z_0_xyz_xxyz, g_z_0_xyz_xxyzz, g_z_0_xyz_xxzz, g_z_0_xyz_xxzzz, g_z_0_xyz_xyyy, g_z_0_xyz_xyyyy, g_z_0_xyz_xyyyz, g_z_0_xyz_xyyz, g_z_0_xyz_xyyzz, g_z_0_xyz_xyzz, g_z_0_xyz_xyzzz, g_z_0_xyz_xzzz, g_z_0_xyz_xzzzz, g_z_0_xyz_yyyy, g_z_0_xyz_yyyz, g_z_0_xyz_yyzz, g_z_0_xyz_yzzz, g_z_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyz_xxxx[k] = -g_z_0_xyz_xxxx[k] * ab_x + g_z_0_xyz_xxxxx[k];

                g_z_0_xxyz_xxxy[k] = -g_z_0_xyz_xxxy[k] * ab_x + g_z_0_xyz_xxxxy[k];

                g_z_0_xxyz_xxxz[k] = -g_z_0_xyz_xxxz[k] * ab_x + g_z_0_xyz_xxxxz[k];

                g_z_0_xxyz_xxyy[k] = -g_z_0_xyz_xxyy[k] * ab_x + g_z_0_xyz_xxxyy[k];

                g_z_0_xxyz_xxyz[k] = -g_z_0_xyz_xxyz[k] * ab_x + g_z_0_xyz_xxxyz[k];

                g_z_0_xxyz_xxzz[k] = -g_z_0_xyz_xxzz[k] * ab_x + g_z_0_xyz_xxxzz[k];

                g_z_0_xxyz_xyyy[k] = -g_z_0_xyz_xyyy[k] * ab_x + g_z_0_xyz_xxyyy[k];

                g_z_0_xxyz_xyyz[k] = -g_z_0_xyz_xyyz[k] * ab_x + g_z_0_xyz_xxyyz[k];

                g_z_0_xxyz_xyzz[k] = -g_z_0_xyz_xyzz[k] * ab_x + g_z_0_xyz_xxyzz[k];

                g_z_0_xxyz_xzzz[k] = -g_z_0_xyz_xzzz[k] * ab_x + g_z_0_xyz_xxzzz[k];

                g_z_0_xxyz_yyyy[k] = -g_z_0_xyz_yyyy[k] * ab_x + g_z_0_xyz_xyyyy[k];

                g_z_0_xxyz_yyyz[k] = -g_z_0_xyz_yyyz[k] * ab_x + g_z_0_xyz_xyyyz[k];

                g_z_0_xxyz_yyzz[k] = -g_z_0_xyz_yyzz[k] * ab_x + g_z_0_xyz_xyyzz[k];

                g_z_0_xxyz_yzzz[k] = -g_z_0_xyz_yzzz[k] * ab_x + g_z_0_xyz_xyzzz[k];

                g_z_0_xxyz_zzzz[k] = -g_z_0_xyz_zzzz[k] * ab_x + g_z_0_xyz_xzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 525 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 526 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 527 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 528 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 529 * ccomps * dcomps);

            auto g_z_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 530 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 531 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 532 * ccomps * dcomps);

            auto g_z_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 533 * ccomps * dcomps);

            auto g_z_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 534 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 535 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 536 * ccomps * dcomps);

            auto g_z_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 537 * ccomps * dcomps);

            auto g_z_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 538 * ccomps * dcomps);

            auto g_z_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzz_xxxx, g_z_0_xxzz_xxxy, g_z_0_xxzz_xxxz, g_z_0_xxzz_xxyy, g_z_0_xxzz_xxyz, g_z_0_xxzz_xxzz, g_z_0_xxzz_xyyy, g_z_0_xxzz_xyyz, g_z_0_xxzz_xyzz, g_z_0_xxzz_xzzz, g_z_0_xxzz_yyyy, g_z_0_xxzz_yyyz, g_z_0_xxzz_yyzz, g_z_0_xxzz_yzzz, g_z_0_xxzz_zzzz, g_z_0_xzz_xxxx, g_z_0_xzz_xxxxx, g_z_0_xzz_xxxxy, g_z_0_xzz_xxxxz, g_z_0_xzz_xxxy, g_z_0_xzz_xxxyy, g_z_0_xzz_xxxyz, g_z_0_xzz_xxxz, g_z_0_xzz_xxxzz, g_z_0_xzz_xxyy, g_z_0_xzz_xxyyy, g_z_0_xzz_xxyyz, g_z_0_xzz_xxyz, g_z_0_xzz_xxyzz, g_z_0_xzz_xxzz, g_z_0_xzz_xxzzz, g_z_0_xzz_xyyy, g_z_0_xzz_xyyyy, g_z_0_xzz_xyyyz, g_z_0_xzz_xyyz, g_z_0_xzz_xyyzz, g_z_0_xzz_xyzz, g_z_0_xzz_xyzzz, g_z_0_xzz_xzzz, g_z_0_xzz_xzzzz, g_z_0_xzz_yyyy, g_z_0_xzz_yyyz, g_z_0_xzz_yyzz, g_z_0_xzz_yzzz, g_z_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzz_xxxx[k] = -g_z_0_xzz_xxxx[k] * ab_x + g_z_0_xzz_xxxxx[k];

                g_z_0_xxzz_xxxy[k] = -g_z_0_xzz_xxxy[k] * ab_x + g_z_0_xzz_xxxxy[k];

                g_z_0_xxzz_xxxz[k] = -g_z_0_xzz_xxxz[k] * ab_x + g_z_0_xzz_xxxxz[k];

                g_z_0_xxzz_xxyy[k] = -g_z_0_xzz_xxyy[k] * ab_x + g_z_0_xzz_xxxyy[k];

                g_z_0_xxzz_xxyz[k] = -g_z_0_xzz_xxyz[k] * ab_x + g_z_0_xzz_xxxyz[k];

                g_z_0_xxzz_xxzz[k] = -g_z_0_xzz_xxzz[k] * ab_x + g_z_0_xzz_xxxzz[k];

                g_z_0_xxzz_xyyy[k] = -g_z_0_xzz_xyyy[k] * ab_x + g_z_0_xzz_xxyyy[k];

                g_z_0_xxzz_xyyz[k] = -g_z_0_xzz_xyyz[k] * ab_x + g_z_0_xzz_xxyyz[k];

                g_z_0_xxzz_xyzz[k] = -g_z_0_xzz_xyzz[k] * ab_x + g_z_0_xzz_xxyzz[k];

                g_z_0_xxzz_xzzz[k] = -g_z_0_xzz_xzzz[k] * ab_x + g_z_0_xzz_xxzzz[k];

                g_z_0_xxzz_yyyy[k] = -g_z_0_xzz_yyyy[k] * ab_x + g_z_0_xzz_xyyyy[k];

                g_z_0_xxzz_yyyz[k] = -g_z_0_xzz_yyyz[k] * ab_x + g_z_0_xzz_xyyyz[k];

                g_z_0_xxzz_yyzz[k] = -g_z_0_xzz_yyzz[k] * ab_x + g_z_0_xzz_xyyzz[k];

                g_z_0_xxzz_yzzz[k] = -g_z_0_xzz_yzzz[k] * ab_x + g_z_0_xzz_xyzzz[k];

                g_z_0_xxzz_zzzz[k] = -g_z_0_xzz_zzzz[k] * ab_x + g_z_0_xzz_xzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 540 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 541 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 542 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 543 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 544 * ccomps * dcomps);

            auto g_z_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 545 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 546 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 547 * ccomps * dcomps);

            auto g_z_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 548 * ccomps * dcomps);

            auto g_z_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 549 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 550 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 551 * ccomps * dcomps);

            auto g_z_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 552 * ccomps * dcomps);

            auto g_z_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 553 * ccomps * dcomps);

            auto g_z_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyy_xxxx, g_z_0_xyyy_xxxy, g_z_0_xyyy_xxxz, g_z_0_xyyy_xxyy, g_z_0_xyyy_xxyz, g_z_0_xyyy_xxzz, g_z_0_xyyy_xyyy, g_z_0_xyyy_xyyz, g_z_0_xyyy_xyzz, g_z_0_xyyy_xzzz, g_z_0_xyyy_yyyy, g_z_0_xyyy_yyyz, g_z_0_xyyy_yyzz, g_z_0_xyyy_yzzz, g_z_0_xyyy_zzzz, g_z_0_yyy_xxxx, g_z_0_yyy_xxxxx, g_z_0_yyy_xxxxy, g_z_0_yyy_xxxxz, g_z_0_yyy_xxxy, g_z_0_yyy_xxxyy, g_z_0_yyy_xxxyz, g_z_0_yyy_xxxz, g_z_0_yyy_xxxzz, g_z_0_yyy_xxyy, g_z_0_yyy_xxyyy, g_z_0_yyy_xxyyz, g_z_0_yyy_xxyz, g_z_0_yyy_xxyzz, g_z_0_yyy_xxzz, g_z_0_yyy_xxzzz, g_z_0_yyy_xyyy, g_z_0_yyy_xyyyy, g_z_0_yyy_xyyyz, g_z_0_yyy_xyyz, g_z_0_yyy_xyyzz, g_z_0_yyy_xyzz, g_z_0_yyy_xyzzz, g_z_0_yyy_xzzz, g_z_0_yyy_xzzzz, g_z_0_yyy_yyyy, g_z_0_yyy_yyyz, g_z_0_yyy_yyzz, g_z_0_yyy_yzzz, g_z_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyy_xxxx[k] = -g_z_0_yyy_xxxx[k] * ab_x + g_z_0_yyy_xxxxx[k];

                g_z_0_xyyy_xxxy[k] = -g_z_0_yyy_xxxy[k] * ab_x + g_z_0_yyy_xxxxy[k];

                g_z_0_xyyy_xxxz[k] = -g_z_0_yyy_xxxz[k] * ab_x + g_z_0_yyy_xxxxz[k];

                g_z_0_xyyy_xxyy[k] = -g_z_0_yyy_xxyy[k] * ab_x + g_z_0_yyy_xxxyy[k];

                g_z_0_xyyy_xxyz[k] = -g_z_0_yyy_xxyz[k] * ab_x + g_z_0_yyy_xxxyz[k];

                g_z_0_xyyy_xxzz[k] = -g_z_0_yyy_xxzz[k] * ab_x + g_z_0_yyy_xxxzz[k];

                g_z_0_xyyy_xyyy[k] = -g_z_0_yyy_xyyy[k] * ab_x + g_z_0_yyy_xxyyy[k];

                g_z_0_xyyy_xyyz[k] = -g_z_0_yyy_xyyz[k] * ab_x + g_z_0_yyy_xxyyz[k];

                g_z_0_xyyy_xyzz[k] = -g_z_0_yyy_xyzz[k] * ab_x + g_z_0_yyy_xxyzz[k];

                g_z_0_xyyy_xzzz[k] = -g_z_0_yyy_xzzz[k] * ab_x + g_z_0_yyy_xxzzz[k];

                g_z_0_xyyy_yyyy[k] = -g_z_0_yyy_yyyy[k] * ab_x + g_z_0_yyy_xyyyy[k];

                g_z_0_xyyy_yyyz[k] = -g_z_0_yyy_yyyz[k] * ab_x + g_z_0_yyy_xyyyz[k];

                g_z_0_xyyy_yyzz[k] = -g_z_0_yyy_yyzz[k] * ab_x + g_z_0_yyy_xyyzz[k];

                g_z_0_xyyy_yzzz[k] = -g_z_0_yyy_yzzz[k] * ab_x + g_z_0_yyy_xyzzz[k];

                g_z_0_xyyy_zzzz[k] = -g_z_0_yyy_zzzz[k] * ab_x + g_z_0_yyy_xzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 555 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 556 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 557 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 558 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 559 * ccomps * dcomps);

            auto g_z_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 566 * ccomps * dcomps);

            auto g_z_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyz_xxxx, g_z_0_xyyz_xxxy, g_z_0_xyyz_xxxz, g_z_0_xyyz_xxyy, g_z_0_xyyz_xxyz, g_z_0_xyyz_xxzz, g_z_0_xyyz_xyyy, g_z_0_xyyz_xyyz, g_z_0_xyyz_xyzz, g_z_0_xyyz_xzzz, g_z_0_xyyz_yyyy, g_z_0_xyyz_yyyz, g_z_0_xyyz_yyzz, g_z_0_xyyz_yzzz, g_z_0_xyyz_zzzz, g_z_0_yyz_xxxx, g_z_0_yyz_xxxxx, g_z_0_yyz_xxxxy, g_z_0_yyz_xxxxz, g_z_0_yyz_xxxy, g_z_0_yyz_xxxyy, g_z_0_yyz_xxxyz, g_z_0_yyz_xxxz, g_z_0_yyz_xxxzz, g_z_0_yyz_xxyy, g_z_0_yyz_xxyyy, g_z_0_yyz_xxyyz, g_z_0_yyz_xxyz, g_z_0_yyz_xxyzz, g_z_0_yyz_xxzz, g_z_0_yyz_xxzzz, g_z_0_yyz_xyyy, g_z_0_yyz_xyyyy, g_z_0_yyz_xyyyz, g_z_0_yyz_xyyz, g_z_0_yyz_xyyzz, g_z_0_yyz_xyzz, g_z_0_yyz_xyzzz, g_z_0_yyz_xzzz, g_z_0_yyz_xzzzz, g_z_0_yyz_yyyy, g_z_0_yyz_yyyz, g_z_0_yyz_yyzz, g_z_0_yyz_yzzz, g_z_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyz_xxxx[k] = -g_z_0_yyz_xxxx[k] * ab_x + g_z_0_yyz_xxxxx[k];

                g_z_0_xyyz_xxxy[k] = -g_z_0_yyz_xxxy[k] * ab_x + g_z_0_yyz_xxxxy[k];

                g_z_0_xyyz_xxxz[k] = -g_z_0_yyz_xxxz[k] * ab_x + g_z_0_yyz_xxxxz[k];

                g_z_0_xyyz_xxyy[k] = -g_z_0_yyz_xxyy[k] * ab_x + g_z_0_yyz_xxxyy[k];

                g_z_0_xyyz_xxyz[k] = -g_z_0_yyz_xxyz[k] * ab_x + g_z_0_yyz_xxxyz[k];

                g_z_0_xyyz_xxzz[k] = -g_z_0_yyz_xxzz[k] * ab_x + g_z_0_yyz_xxxzz[k];

                g_z_0_xyyz_xyyy[k] = -g_z_0_yyz_xyyy[k] * ab_x + g_z_0_yyz_xxyyy[k];

                g_z_0_xyyz_xyyz[k] = -g_z_0_yyz_xyyz[k] * ab_x + g_z_0_yyz_xxyyz[k];

                g_z_0_xyyz_xyzz[k] = -g_z_0_yyz_xyzz[k] * ab_x + g_z_0_yyz_xxyzz[k];

                g_z_0_xyyz_xzzz[k] = -g_z_0_yyz_xzzz[k] * ab_x + g_z_0_yyz_xxzzz[k];

                g_z_0_xyyz_yyyy[k] = -g_z_0_yyz_yyyy[k] * ab_x + g_z_0_yyz_xyyyy[k];

                g_z_0_xyyz_yyyz[k] = -g_z_0_yyz_yyyz[k] * ab_x + g_z_0_yyz_xyyyz[k];

                g_z_0_xyyz_yyzz[k] = -g_z_0_yyz_yyzz[k] * ab_x + g_z_0_yyz_xyyzz[k];

                g_z_0_xyyz_yzzz[k] = -g_z_0_yyz_yzzz[k] * ab_x + g_z_0_yyz_xyzzz[k];

                g_z_0_xyyz_zzzz[k] = -g_z_0_yyz_zzzz[k] * ab_x + g_z_0_yyz_xzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 575 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 579 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzz_xxxx, g_z_0_xyzz_xxxy, g_z_0_xyzz_xxxz, g_z_0_xyzz_xxyy, g_z_0_xyzz_xxyz, g_z_0_xyzz_xxzz, g_z_0_xyzz_xyyy, g_z_0_xyzz_xyyz, g_z_0_xyzz_xyzz, g_z_0_xyzz_xzzz, g_z_0_xyzz_yyyy, g_z_0_xyzz_yyyz, g_z_0_xyzz_yyzz, g_z_0_xyzz_yzzz, g_z_0_xyzz_zzzz, g_z_0_yzz_xxxx, g_z_0_yzz_xxxxx, g_z_0_yzz_xxxxy, g_z_0_yzz_xxxxz, g_z_0_yzz_xxxy, g_z_0_yzz_xxxyy, g_z_0_yzz_xxxyz, g_z_0_yzz_xxxz, g_z_0_yzz_xxxzz, g_z_0_yzz_xxyy, g_z_0_yzz_xxyyy, g_z_0_yzz_xxyyz, g_z_0_yzz_xxyz, g_z_0_yzz_xxyzz, g_z_0_yzz_xxzz, g_z_0_yzz_xxzzz, g_z_0_yzz_xyyy, g_z_0_yzz_xyyyy, g_z_0_yzz_xyyyz, g_z_0_yzz_xyyz, g_z_0_yzz_xyyzz, g_z_0_yzz_xyzz, g_z_0_yzz_xyzzz, g_z_0_yzz_xzzz, g_z_0_yzz_xzzzz, g_z_0_yzz_yyyy, g_z_0_yzz_yyyz, g_z_0_yzz_yyzz, g_z_0_yzz_yzzz, g_z_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzz_xxxx[k] = -g_z_0_yzz_xxxx[k] * ab_x + g_z_0_yzz_xxxxx[k];

                g_z_0_xyzz_xxxy[k] = -g_z_0_yzz_xxxy[k] * ab_x + g_z_0_yzz_xxxxy[k];

                g_z_0_xyzz_xxxz[k] = -g_z_0_yzz_xxxz[k] * ab_x + g_z_0_yzz_xxxxz[k];

                g_z_0_xyzz_xxyy[k] = -g_z_0_yzz_xxyy[k] * ab_x + g_z_0_yzz_xxxyy[k];

                g_z_0_xyzz_xxyz[k] = -g_z_0_yzz_xxyz[k] * ab_x + g_z_0_yzz_xxxyz[k];

                g_z_0_xyzz_xxzz[k] = -g_z_0_yzz_xxzz[k] * ab_x + g_z_0_yzz_xxxzz[k];

                g_z_0_xyzz_xyyy[k] = -g_z_0_yzz_xyyy[k] * ab_x + g_z_0_yzz_xxyyy[k];

                g_z_0_xyzz_xyyz[k] = -g_z_0_yzz_xyyz[k] * ab_x + g_z_0_yzz_xxyyz[k];

                g_z_0_xyzz_xyzz[k] = -g_z_0_yzz_xyzz[k] * ab_x + g_z_0_yzz_xxyzz[k];

                g_z_0_xyzz_xzzz[k] = -g_z_0_yzz_xzzz[k] * ab_x + g_z_0_yzz_xxzzz[k];

                g_z_0_xyzz_yyyy[k] = -g_z_0_yzz_yyyy[k] * ab_x + g_z_0_yzz_xyyyy[k];

                g_z_0_xyzz_yyyz[k] = -g_z_0_yzz_yyyz[k] * ab_x + g_z_0_yzz_xyyyz[k];

                g_z_0_xyzz_yyzz[k] = -g_z_0_yzz_yyzz[k] * ab_x + g_z_0_yzz_xyyzz[k];

                g_z_0_xyzz_yzzz[k] = -g_z_0_yzz_yzzz[k] * ab_x + g_z_0_yzz_xyzzz[k];

                g_z_0_xyzz_zzzz[k] = -g_z_0_yzz_zzzz[k] * ab_x + g_z_0_yzz_xzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 587 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 589 * ccomps * dcomps);

            auto g_z_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzz_xxxx, g_z_0_xzzz_xxxy, g_z_0_xzzz_xxxz, g_z_0_xzzz_xxyy, g_z_0_xzzz_xxyz, g_z_0_xzzz_xxzz, g_z_0_xzzz_xyyy, g_z_0_xzzz_xyyz, g_z_0_xzzz_xyzz, g_z_0_xzzz_xzzz, g_z_0_xzzz_yyyy, g_z_0_xzzz_yyyz, g_z_0_xzzz_yyzz, g_z_0_xzzz_yzzz, g_z_0_xzzz_zzzz, g_z_0_zzz_xxxx, g_z_0_zzz_xxxxx, g_z_0_zzz_xxxxy, g_z_0_zzz_xxxxz, g_z_0_zzz_xxxy, g_z_0_zzz_xxxyy, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxz, g_z_0_zzz_xxxzz, g_z_0_zzz_xxyy, g_z_0_zzz_xxyyy, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxzz, g_z_0_zzz_xxzzz, g_z_0_zzz_xyyy, g_z_0_zzz_xyyyy, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xzzz, g_z_0_zzz_xzzzz, g_z_0_zzz_yyyy, g_z_0_zzz_yyyz, g_z_0_zzz_yyzz, g_z_0_zzz_yzzz, g_z_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzz_xxxx[k] = -g_z_0_zzz_xxxx[k] * ab_x + g_z_0_zzz_xxxxx[k];

                g_z_0_xzzz_xxxy[k] = -g_z_0_zzz_xxxy[k] * ab_x + g_z_0_zzz_xxxxy[k];

                g_z_0_xzzz_xxxz[k] = -g_z_0_zzz_xxxz[k] * ab_x + g_z_0_zzz_xxxxz[k];

                g_z_0_xzzz_xxyy[k] = -g_z_0_zzz_xxyy[k] * ab_x + g_z_0_zzz_xxxyy[k];

                g_z_0_xzzz_xxyz[k] = -g_z_0_zzz_xxyz[k] * ab_x + g_z_0_zzz_xxxyz[k];

                g_z_0_xzzz_xxzz[k] = -g_z_0_zzz_xxzz[k] * ab_x + g_z_0_zzz_xxxzz[k];

                g_z_0_xzzz_xyyy[k] = -g_z_0_zzz_xyyy[k] * ab_x + g_z_0_zzz_xxyyy[k];

                g_z_0_xzzz_xyyz[k] = -g_z_0_zzz_xyyz[k] * ab_x + g_z_0_zzz_xxyyz[k];

                g_z_0_xzzz_xyzz[k] = -g_z_0_zzz_xyzz[k] * ab_x + g_z_0_zzz_xxyzz[k];

                g_z_0_xzzz_xzzz[k] = -g_z_0_zzz_xzzz[k] * ab_x + g_z_0_zzz_xxzzz[k];

                g_z_0_xzzz_yyyy[k] = -g_z_0_zzz_yyyy[k] * ab_x + g_z_0_zzz_xyyyy[k];

                g_z_0_xzzz_yyyz[k] = -g_z_0_zzz_yyyz[k] * ab_x + g_z_0_zzz_xyyyz[k];

                g_z_0_xzzz_yyzz[k] = -g_z_0_zzz_yyzz[k] * ab_x + g_z_0_zzz_xyyzz[k];

                g_z_0_xzzz_yzzz[k] = -g_z_0_zzz_yzzz[k] * ab_x + g_z_0_zzz_xyzzz[k];

                g_z_0_xzzz_zzzz[k] = -g_z_0_zzz_zzzz[k] * ab_x + g_z_0_zzz_xzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 608 * ccomps * dcomps);

            auto g_z_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 609 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 611 * ccomps * dcomps);

            auto g_z_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyy_xxxx, g_z_0_yyy_xxxxy, g_z_0_yyy_xxxy, g_z_0_yyy_xxxyy, g_z_0_yyy_xxxyz, g_z_0_yyy_xxxz, g_z_0_yyy_xxyy, g_z_0_yyy_xxyyy, g_z_0_yyy_xxyyz, g_z_0_yyy_xxyz, g_z_0_yyy_xxyzz, g_z_0_yyy_xxzz, g_z_0_yyy_xyyy, g_z_0_yyy_xyyyy, g_z_0_yyy_xyyyz, g_z_0_yyy_xyyz, g_z_0_yyy_xyyzz, g_z_0_yyy_xyzz, g_z_0_yyy_xyzzz, g_z_0_yyy_xzzz, g_z_0_yyy_yyyy, g_z_0_yyy_yyyyy, g_z_0_yyy_yyyyz, g_z_0_yyy_yyyz, g_z_0_yyy_yyyzz, g_z_0_yyy_yyzz, g_z_0_yyy_yyzzz, g_z_0_yyy_yzzz, g_z_0_yyy_yzzzz, g_z_0_yyy_zzzz, g_z_0_yyyy_xxxx, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxxz, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxzz, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xzzz, g_z_0_yyyy_yyyy, g_z_0_yyyy_yyyz, g_z_0_yyyy_yyzz, g_z_0_yyyy_yzzz, g_z_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyy_xxxx[k] = -g_z_0_yyy_xxxx[k] * ab_y + g_z_0_yyy_xxxxy[k];

                g_z_0_yyyy_xxxy[k] = -g_z_0_yyy_xxxy[k] * ab_y + g_z_0_yyy_xxxyy[k];

                g_z_0_yyyy_xxxz[k] = -g_z_0_yyy_xxxz[k] * ab_y + g_z_0_yyy_xxxyz[k];

                g_z_0_yyyy_xxyy[k] = -g_z_0_yyy_xxyy[k] * ab_y + g_z_0_yyy_xxyyy[k];

                g_z_0_yyyy_xxyz[k] = -g_z_0_yyy_xxyz[k] * ab_y + g_z_0_yyy_xxyyz[k];

                g_z_0_yyyy_xxzz[k] = -g_z_0_yyy_xxzz[k] * ab_y + g_z_0_yyy_xxyzz[k];

                g_z_0_yyyy_xyyy[k] = -g_z_0_yyy_xyyy[k] * ab_y + g_z_0_yyy_xyyyy[k];

                g_z_0_yyyy_xyyz[k] = -g_z_0_yyy_xyyz[k] * ab_y + g_z_0_yyy_xyyyz[k];

                g_z_0_yyyy_xyzz[k] = -g_z_0_yyy_xyzz[k] * ab_y + g_z_0_yyy_xyyzz[k];

                g_z_0_yyyy_xzzz[k] = -g_z_0_yyy_xzzz[k] * ab_y + g_z_0_yyy_xyzzz[k];

                g_z_0_yyyy_yyyy[k] = -g_z_0_yyy_yyyy[k] * ab_y + g_z_0_yyy_yyyyy[k];

                g_z_0_yyyy_yyyz[k] = -g_z_0_yyy_yyyz[k] * ab_y + g_z_0_yyy_yyyyz[k];

                g_z_0_yyyy_yyzz[k] = -g_z_0_yyy_yyzz[k] * ab_y + g_z_0_yyy_yyyzz[k];

                g_z_0_yyyy_yzzz[k] = -g_z_0_yyy_yzzz[k] * ab_y + g_z_0_yyy_yyzzz[k];

                g_z_0_yyyy_zzzz[k] = -g_z_0_yyy_zzzz[k] * ab_y + g_z_0_yyy_yzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 615 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 619 * ccomps * dcomps);

            auto g_z_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyz_xxxx, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxxz, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxzz, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xzzz, g_z_0_yyyz_yyyy, g_z_0_yyyz_yyyz, g_z_0_yyyz_yyzz, g_z_0_yyyz_yzzz, g_z_0_yyyz_zzzz, g_z_0_yyz_xxxx, g_z_0_yyz_xxxxy, g_z_0_yyz_xxxy, g_z_0_yyz_xxxyy, g_z_0_yyz_xxxyz, g_z_0_yyz_xxxz, g_z_0_yyz_xxyy, g_z_0_yyz_xxyyy, g_z_0_yyz_xxyyz, g_z_0_yyz_xxyz, g_z_0_yyz_xxyzz, g_z_0_yyz_xxzz, g_z_0_yyz_xyyy, g_z_0_yyz_xyyyy, g_z_0_yyz_xyyyz, g_z_0_yyz_xyyz, g_z_0_yyz_xyyzz, g_z_0_yyz_xyzz, g_z_0_yyz_xyzzz, g_z_0_yyz_xzzz, g_z_0_yyz_yyyy, g_z_0_yyz_yyyyy, g_z_0_yyz_yyyyz, g_z_0_yyz_yyyz, g_z_0_yyz_yyyzz, g_z_0_yyz_yyzz, g_z_0_yyz_yyzzz, g_z_0_yyz_yzzz, g_z_0_yyz_yzzzz, g_z_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyz_xxxx[k] = -g_z_0_yyz_xxxx[k] * ab_y + g_z_0_yyz_xxxxy[k];

                g_z_0_yyyz_xxxy[k] = -g_z_0_yyz_xxxy[k] * ab_y + g_z_0_yyz_xxxyy[k];

                g_z_0_yyyz_xxxz[k] = -g_z_0_yyz_xxxz[k] * ab_y + g_z_0_yyz_xxxyz[k];

                g_z_0_yyyz_xxyy[k] = -g_z_0_yyz_xxyy[k] * ab_y + g_z_0_yyz_xxyyy[k];

                g_z_0_yyyz_xxyz[k] = -g_z_0_yyz_xxyz[k] * ab_y + g_z_0_yyz_xxyyz[k];

                g_z_0_yyyz_xxzz[k] = -g_z_0_yyz_xxzz[k] * ab_y + g_z_0_yyz_xxyzz[k];

                g_z_0_yyyz_xyyy[k] = -g_z_0_yyz_xyyy[k] * ab_y + g_z_0_yyz_xyyyy[k];

                g_z_0_yyyz_xyyz[k] = -g_z_0_yyz_xyyz[k] * ab_y + g_z_0_yyz_xyyyz[k];

                g_z_0_yyyz_xyzz[k] = -g_z_0_yyz_xyzz[k] * ab_y + g_z_0_yyz_xyyzz[k];

                g_z_0_yyyz_xzzz[k] = -g_z_0_yyz_xzzz[k] * ab_y + g_z_0_yyz_xyzzz[k];

                g_z_0_yyyz_yyyy[k] = -g_z_0_yyz_yyyy[k] * ab_y + g_z_0_yyz_yyyyy[k];

                g_z_0_yyyz_yyyz[k] = -g_z_0_yyz_yyyz[k] * ab_y + g_z_0_yyz_yyyyz[k];

                g_z_0_yyyz_yyzz[k] = -g_z_0_yyz_yyzz[k] * ab_y + g_z_0_yyz_yyyzz[k];

                g_z_0_yyyz_yzzz[k] = -g_z_0_yyz_yzzz[k] * ab_y + g_z_0_yyz_yyzzz[k];

                g_z_0_yyyz_zzzz[k] = -g_z_0_yyz_zzzz[k] * ab_y + g_z_0_yyz_yzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 630 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 631 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 632 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 633 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 634 * ccomps * dcomps);

            auto g_z_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 635 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 636 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 637 * ccomps * dcomps);

            auto g_z_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 638 * ccomps * dcomps);

            auto g_z_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 639 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 640 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 641 * ccomps * dcomps);

            auto g_z_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 642 * ccomps * dcomps);

            auto g_z_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 643 * ccomps * dcomps);

            auto g_z_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzz_xxxx, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxxz, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxzz, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xzzz, g_z_0_yyzz_yyyy, g_z_0_yyzz_yyyz, g_z_0_yyzz_yyzz, g_z_0_yyzz_yzzz, g_z_0_yyzz_zzzz, g_z_0_yzz_xxxx, g_z_0_yzz_xxxxy, g_z_0_yzz_xxxy, g_z_0_yzz_xxxyy, g_z_0_yzz_xxxyz, g_z_0_yzz_xxxz, g_z_0_yzz_xxyy, g_z_0_yzz_xxyyy, g_z_0_yzz_xxyyz, g_z_0_yzz_xxyz, g_z_0_yzz_xxyzz, g_z_0_yzz_xxzz, g_z_0_yzz_xyyy, g_z_0_yzz_xyyyy, g_z_0_yzz_xyyyz, g_z_0_yzz_xyyz, g_z_0_yzz_xyyzz, g_z_0_yzz_xyzz, g_z_0_yzz_xyzzz, g_z_0_yzz_xzzz, g_z_0_yzz_yyyy, g_z_0_yzz_yyyyy, g_z_0_yzz_yyyyz, g_z_0_yzz_yyyz, g_z_0_yzz_yyyzz, g_z_0_yzz_yyzz, g_z_0_yzz_yyzzz, g_z_0_yzz_yzzz, g_z_0_yzz_yzzzz, g_z_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzz_xxxx[k] = -g_z_0_yzz_xxxx[k] * ab_y + g_z_0_yzz_xxxxy[k];

                g_z_0_yyzz_xxxy[k] = -g_z_0_yzz_xxxy[k] * ab_y + g_z_0_yzz_xxxyy[k];

                g_z_0_yyzz_xxxz[k] = -g_z_0_yzz_xxxz[k] * ab_y + g_z_0_yzz_xxxyz[k];

                g_z_0_yyzz_xxyy[k] = -g_z_0_yzz_xxyy[k] * ab_y + g_z_0_yzz_xxyyy[k];

                g_z_0_yyzz_xxyz[k] = -g_z_0_yzz_xxyz[k] * ab_y + g_z_0_yzz_xxyyz[k];

                g_z_0_yyzz_xxzz[k] = -g_z_0_yzz_xxzz[k] * ab_y + g_z_0_yzz_xxyzz[k];

                g_z_0_yyzz_xyyy[k] = -g_z_0_yzz_xyyy[k] * ab_y + g_z_0_yzz_xyyyy[k];

                g_z_0_yyzz_xyyz[k] = -g_z_0_yzz_xyyz[k] * ab_y + g_z_0_yzz_xyyyz[k];

                g_z_0_yyzz_xyzz[k] = -g_z_0_yzz_xyzz[k] * ab_y + g_z_0_yzz_xyyzz[k];

                g_z_0_yyzz_xzzz[k] = -g_z_0_yzz_xzzz[k] * ab_y + g_z_0_yzz_xyzzz[k];

                g_z_0_yyzz_yyyy[k] = -g_z_0_yzz_yyyy[k] * ab_y + g_z_0_yzz_yyyyy[k];

                g_z_0_yyzz_yyyz[k] = -g_z_0_yzz_yyyz[k] * ab_y + g_z_0_yzz_yyyyz[k];

                g_z_0_yyzz_yyzz[k] = -g_z_0_yzz_yyzz[k] * ab_y + g_z_0_yzz_yyyzz[k];

                g_z_0_yyzz_yzzz[k] = -g_z_0_yzz_yzzz[k] * ab_y + g_z_0_yzz_yyzzz[k];

                g_z_0_yyzz_zzzz[k] = -g_z_0_yzz_zzzz[k] * ab_y + g_z_0_yzz_yzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 645 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 646 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 647 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 648 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 649 * ccomps * dcomps);

            auto g_z_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 650 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 651 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 652 * ccomps * dcomps);

            auto g_z_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 653 * ccomps * dcomps);

            auto g_z_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 654 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 655 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 656 * ccomps * dcomps);

            auto g_z_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 657 * ccomps * dcomps);

            auto g_z_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 658 * ccomps * dcomps);

            auto g_z_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzz_xxxx, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxxz, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxzz, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xzzz, g_z_0_yzzz_yyyy, g_z_0_yzzz_yyyz, g_z_0_yzzz_yyzz, g_z_0_yzzz_yzzz, g_z_0_yzzz_zzzz, g_z_0_zzz_xxxx, g_z_0_zzz_xxxxy, g_z_0_zzz_xxxy, g_z_0_zzz_xxxyy, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxz, g_z_0_zzz_xxyy, g_z_0_zzz_xxyyy, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxzz, g_z_0_zzz_xyyy, g_z_0_zzz_xyyyy, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xzzz, g_z_0_zzz_yyyy, g_z_0_zzz_yyyyy, g_z_0_zzz_yyyyz, g_z_0_zzz_yyyz, g_z_0_zzz_yyyzz, g_z_0_zzz_yyzz, g_z_0_zzz_yyzzz, g_z_0_zzz_yzzz, g_z_0_zzz_yzzzz, g_z_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzz_xxxx[k] = -g_z_0_zzz_xxxx[k] * ab_y + g_z_0_zzz_xxxxy[k];

                g_z_0_yzzz_xxxy[k] = -g_z_0_zzz_xxxy[k] * ab_y + g_z_0_zzz_xxxyy[k];

                g_z_0_yzzz_xxxz[k] = -g_z_0_zzz_xxxz[k] * ab_y + g_z_0_zzz_xxxyz[k];

                g_z_0_yzzz_xxyy[k] = -g_z_0_zzz_xxyy[k] * ab_y + g_z_0_zzz_xxyyy[k];

                g_z_0_yzzz_xxyz[k] = -g_z_0_zzz_xxyz[k] * ab_y + g_z_0_zzz_xxyyz[k];

                g_z_0_yzzz_xxzz[k] = -g_z_0_zzz_xxzz[k] * ab_y + g_z_0_zzz_xxyzz[k];

                g_z_0_yzzz_xyyy[k] = -g_z_0_zzz_xyyy[k] * ab_y + g_z_0_zzz_xyyyy[k];

                g_z_0_yzzz_xyyz[k] = -g_z_0_zzz_xyyz[k] * ab_y + g_z_0_zzz_xyyyz[k];

                g_z_0_yzzz_xyzz[k] = -g_z_0_zzz_xyzz[k] * ab_y + g_z_0_zzz_xyyzz[k];

                g_z_0_yzzz_xzzz[k] = -g_z_0_zzz_xzzz[k] * ab_y + g_z_0_zzz_xyzzz[k];

                g_z_0_yzzz_yyyy[k] = -g_z_0_zzz_yyyy[k] * ab_y + g_z_0_zzz_yyyyy[k];

                g_z_0_yzzz_yyyz[k] = -g_z_0_zzz_yyyz[k] * ab_y + g_z_0_zzz_yyyyz[k];

                g_z_0_yzzz_yyzz[k] = -g_z_0_zzz_yyzz[k] * ab_y + g_z_0_zzz_yyyzz[k];

                g_z_0_yzzz_yzzz[k] = -g_z_0_zzz_yzzz[k] * ab_y + g_z_0_zzz_yyzzz[k];

                g_z_0_yzzz_zzzz[k] = -g_z_0_zzz_zzzz[k] * ab_y + g_z_0_zzz_yzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 660 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 661 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 662 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 663 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 664 * ccomps * dcomps);

            auto g_z_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 665 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 666 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 667 * ccomps * dcomps);

            auto g_z_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 668 * ccomps * dcomps);

            auto g_z_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 669 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 670 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 671 * ccomps * dcomps);

            auto g_z_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 672 * ccomps * dcomps);

            auto g_z_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 673 * ccomps * dcomps);

            auto g_z_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_xxxx, g_z_0_zzz_xxxxz, g_z_0_zzz_xxxy, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxz, g_z_0_zzz_xxxzz, g_z_0_zzz_xxyy, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxzz, g_z_0_zzz_xxzzz, g_z_0_zzz_xyyy, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xzzz, g_z_0_zzz_xzzzz, g_z_0_zzz_yyyy, g_z_0_zzz_yyyyz, g_z_0_zzz_yyyz, g_z_0_zzz_yyyzz, g_z_0_zzz_yyzz, g_z_0_zzz_yyzzz, g_z_0_zzz_yzzz, g_z_0_zzz_yzzzz, g_z_0_zzz_zzzz, g_z_0_zzz_zzzzz, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzzz, g_zzz_xxxx, g_zzz_xxxy, g_zzz_xxxz, g_zzz_xxyy, g_zzz_xxyz, g_zzz_xxzz, g_zzz_xyyy, g_zzz_xyyz, g_zzz_xyzz, g_zzz_xzzz, g_zzz_yyyy, g_zzz_yyyz, g_zzz_yyzz, g_zzz_yzzz, g_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzz_xxxx[k] = -g_zzz_xxxx[k] - g_z_0_zzz_xxxx[k] * ab_z + g_z_0_zzz_xxxxz[k];

                g_z_0_zzzz_xxxy[k] = -g_zzz_xxxy[k] - g_z_0_zzz_xxxy[k] * ab_z + g_z_0_zzz_xxxyz[k];

                g_z_0_zzzz_xxxz[k] = -g_zzz_xxxz[k] - g_z_0_zzz_xxxz[k] * ab_z + g_z_0_zzz_xxxzz[k];

                g_z_0_zzzz_xxyy[k] = -g_zzz_xxyy[k] - g_z_0_zzz_xxyy[k] * ab_z + g_z_0_zzz_xxyyz[k];

                g_z_0_zzzz_xxyz[k] = -g_zzz_xxyz[k] - g_z_0_zzz_xxyz[k] * ab_z + g_z_0_zzz_xxyzz[k];

                g_z_0_zzzz_xxzz[k] = -g_zzz_xxzz[k] - g_z_0_zzz_xxzz[k] * ab_z + g_z_0_zzz_xxzzz[k];

                g_z_0_zzzz_xyyy[k] = -g_zzz_xyyy[k] - g_z_0_zzz_xyyy[k] * ab_z + g_z_0_zzz_xyyyz[k];

                g_z_0_zzzz_xyyz[k] = -g_zzz_xyyz[k] - g_z_0_zzz_xyyz[k] * ab_z + g_z_0_zzz_xyyzz[k];

                g_z_0_zzzz_xyzz[k] = -g_zzz_xyzz[k] - g_z_0_zzz_xyzz[k] * ab_z + g_z_0_zzz_xyzzz[k];

                g_z_0_zzzz_xzzz[k] = -g_zzz_xzzz[k] - g_z_0_zzz_xzzz[k] * ab_z + g_z_0_zzz_xzzzz[k];

                g_z_0_zzzz_yyyy[k] = -g_zzz_yyyy[k] - g_z_0_zzz_yyyy[k] * ab_z + g_z_0_zzz_yyyyz[k];

                g_z_0_zzzz_yyyz[k] = -g_zzz_yyyz[k] - g_z_0_zzz_yyyz[k] * ab_z + g_z_0_zzz_yyyzz[k];

                g_z_0_zzzz_yyzz[k] = -g_zzz_yyzz[k] - g_z_0_zzz_yyzz[k] * ab_z + g_z_0_zzz_yyzzz[k];

                g_z_0_zzzz_yzzz[k] = -g_zzz_yzzz[k] - g_z_0_zzz_yzzz[k] * ab_z + g_z_0_zzz_yzzzz[k];

                g_z_0_zzzz_zzzz[k] = -g_zzz_zzzz[k] - g_z_0_zzz_zzzz[k] * ab_z + g_z_0_zzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

