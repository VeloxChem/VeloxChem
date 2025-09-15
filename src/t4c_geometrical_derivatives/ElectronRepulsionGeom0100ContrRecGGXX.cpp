#include "ElectronRepulsionGeom0100ContrRecGGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_ggxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_ggxx,
                                            const size_t idx_fgxx,
                                            const size_t idx_geom_01_fgxx,
                                            const size_t idx_geom_01_fhxx,
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

            const auto fg_geom_01_off = idx_geom_01_fgxx + i * dcomps + j;

            auto g_0_x_xxx_xxxx = cbuffer.data(fg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxxy = cbuffer.data(fg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxxz = cbuffer.data(fg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xxyy = cbuffer.data(fg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xxyz = cbuffer.data(fg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xxzz = cbuffer.data(fg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_xyyy = cbuffer.data(fg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_xyyz = cbuffer.data(fg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_xyzz = cbuffer.data(fg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_xzzz = cbuffer.data(fg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxx_yyyy = cbuffer.data(fg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxx_yyyz = cbuffer.data(fg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxx_yyzz = cbuffer.data(fg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxx_yzzz = cbuffer.data(fg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxx_zzzz = cbuffer.data(fg_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxy_xxxx = cbuffer.data(fg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxy_xxxy = cbuffer.data(fg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxy_xxxz = cbuffer.data(fg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxy_xxyy = cbuffer.data(fg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxy_xxyz = cbuffer.data(fg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxy_xxzz = cbuffer.data(fg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxy_xyyy = cbuffer.data(fg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxy_xyyz = cbuffer.data(fg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxy_xyzz = cbuffer.data(fg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxy_xzzz = cbuffer.data(fg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxy_yyyy = cbuffer.data(fg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxy_yyyz = cbuffer.data(fg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxy_yyzz = cbuffer.data(fg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxy_yzzz = cbuffer.data(fg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxy_zzzz = cbuffer.data(fg_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxz_xxxx = cbuffer.data(fg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxz_xxxy = cbuffer.data(fg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxz_xxxz = cbuffer.data(fg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxz_xxyy = cbuffer.data(fg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxz_xxyz = cbuffer.data(fg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxz_xxzz = cbuffer.data(fg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxz_xyyy = cbuffer.data(fg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxz_xyyz = cbuffer.data(fg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxz_xyzz = cbuffer.data(fg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxz_xzzz = cbuffer.data(fg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxz_yyyy = cbuffer.data(fg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxz_yyyz = cbuffer.data(fg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxz_yyzz = cbuffer.data(fg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxz_yzzz = cbuffer.data(fg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxz_zzzz = cbuffer.data(fg_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xyy_xxxx = cbuffer.data(fg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyy_xxxy = cbuffer.data(fg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyy_xxxz = cbuffer.data(fg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xyy_xxyy = cbuffer.data(fg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyy_xxyz = cbuffer.data(fg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xyy_xxzz = cbuffer.data(fg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xyy_xyyy = cbuffer.data(fg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xyy_xyyz = cbuffer.data(fg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xyy_xyzz = cbuffer.data(fg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xyy_xzzz = cbuffer.data(fg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xyy_yyyy = cbuffer.data(fg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xyy_yyyz = cbuffer.data(fg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xyy_yyzz = cbuffer.data(fg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xyy_yzzz = cbuffer.data(fg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xyy_zzzz = cbuffer.data(fg_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xyz_xxxx = cbuffer.data(fg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xyz_xxxy = cbuffer.data(fg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xyz_xxxz = cbuffer.data(fg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyz_xxyy = cbuffer.data(fg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyz_xxyz = cbuffer.data(fg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyz_xxzz = cbuffer.data(fg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xyz_xyyy = cbuffer.data(fg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyz_xyyz = cbuffer.data(fg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyz_xyzz = cbuffer.data(fg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyz_xzzz = cbuffer.data(fg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyz_yyyy = cbuffer.data(fg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyz_yyyz = cbuffer.data(fg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xyz_yyzz = cbuffer.data(fg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyz_yzzz = cbuffer.data(fg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyz_zzzz = cbuffer.data(fg_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xzz_xxxx = cbuffer.data(fg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xzz_xxxy = cbuffer.data(fg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xzz_xxxz = cbuffer.data(fg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xzz_xxyy = cbuffer.data(fg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xzz_xxyz = cbuffer.data(fg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xzz_xxzz = cbuffer.data(fg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xzz_xyyy = cbuffer.data(fg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xzz_xyyz = cbuffer.data(fg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xzz_xyzz = cbuffer.data(fg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xzz_xzzz = cbuffer.data(fg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xzz_yyyy = cbuffer.data(fg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xzz_yyyz = cbuffer.data(fg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xzz_yyzz = cbuffer.data(fg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xzz_yzzz = cbuffer.data(fg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xzz_zzzz = cbuffer.data(fg_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_yyy_xxxx = cbuffer.data(fg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_yyy_xxxy = cbuffer.data(fg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_yyy_xxxz = cbuffer.data(fg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_yyy_xxyy = cbuffer.data(fg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_yyy_xxyz = cbuffer.data(fg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_yyy_xxzz = cbuffer.data(fg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_yyy_xyyy = cbuffer.data(fg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_yyy_xyyz = cbuffer.data(fg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_yyy_xyzz = cbuffer.data(fg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_yyy_xzzz = cbuffer.data(fg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yyy_yyyy = cbuffer.data(fg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yyy_yyyz = cbuffer.data(fg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_yyy_yyzz = cbuffer.data(fg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yyy_yzzz = cbuffer.data(fg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yyy_zzzz = cbuffer.data(fg_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_yyz_xxxx = cbuffer.data(fg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_yyz_xxxy = cbuffer.data(fg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_yyz_xxxz = cbuffer.data(fg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_yyz_xxyy = cbuffer.data(fg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yyz_xxyz = cbuffer.data(fg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_yyz_xxzz = cbuffer.data(fg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yyz_xyyy = cbuffer.data(fg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yyz_xyyz = cbuffer.data(fg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yyz_xyzz = cbuffer.data(fg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_yyz_xzzz = cbuffer.data(fg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yyz_yyyy = cbuffer.data(fg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yyz_yyyz = cbuffer.data(fg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yyz_yyzz = cbuffer.data(fg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yyz_yzzz = cbuffer.data(fg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yyz_zzzz = cbuffer.data(fg_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_yzz_xxxx = cbuffer.data(fg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_yzz_xxxy = cbuffer.data(fg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_yzz_xxxz = cbuffer.data(fg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_yzz_xxyy = cbuffer.data(fg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_yzz_xxyz = cbuffer.data(fg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_yzz_xxzz = cbuffer.data(fg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yzz_xyyy = cbuffer.data(fg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yzz_xyyz = cbuffer.data(fg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yzz_xyzz = cbuffer.data(fg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yzz_xzzz = cbuffer.data(fg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yzz_yyyy = cbuffer.data(fg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yzz_yyyz = cbuffer.data(fg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yzz_yyzz = cbuffer.data(fg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yzz_yzzz = cbuffer.data(fg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yzz_zzzz = cbuffer.data(fg_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_zzz_xxxx = cbuffer.data(fg_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_zzz_xxxy = cbuffer.data(fg_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_zzz_xxxz = cbuffer.data(fg_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_zzz_xxyy = cbuffer.data(fg_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_zzz_xxyz = cbuffer.data(fg_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_zzz_xxzz = cbuffer.data(fg_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_zzz_xyyy = cbuffer.data(fg_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_zzz_xyyz = cbuffer.data(fg_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_zzz_xyzz = cbuffer.data(fg_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_zzz_xzzz = cbuffer.data(fg_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_zzz_yyyy = cbuffer.data(fg_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_zzz_yyyz = cbuffer.data(fg_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_zzz_yyzz = cbuffer.data(fg_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_zzz_yzzz = cbuffer.data(fg_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_zzz_zzzz = cbuffer.data(fg_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_xxx_xxxx = cbuffer.data(fg_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xxx_xxxy = cbuffer.data(fg_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xxx_xxxz = cbuffer.data(fg_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xxx_xxyy = cbuffer.data(fg_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xxx_xxyz = cbuffer.data(fg_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xxx_xxzz = cbuffer.data(fg_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xxx_xyyy = cbuffer.data(fg_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xxx_xyyz = cbuffer.data(fg_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xxx_xyzz = cbuffer.data(fg_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xxx_xzzz = cbuffer.data(fg_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xxx_yyyy = cbuffer.data(fg_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xxx_yyyz = cbuffer.data(fg_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_xxx_yyzz = cbuffer.data(fg_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xxx_yzzz = cbuffer.data(fg_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xxx_zzzz = cbuffer.data(fg_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_xxy_xxxx = cbuffer.data(fg_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xxy_xxxy = cbuffer.data(fg_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xxy_xxxz = cbuffer.data(fg_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xxy_xxyy = cbuffer.data(fg_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxy_xxyz = cbuffer.data(fg_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxy_xxzz = cbuffer.data(fg_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxy_xyyy = cbuffer.data(fg_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxy_xyyz = cbuffer.data(fg_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxy_xyzz = cbuffer.data(fg_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xxy_xzzz = cbuffer.data(fg_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxy_yyyy = cbuffer.data(fg_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxy_yyyz = cbuffer.data(fg_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxy_yyzz = cbuffer.data(fg_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxy_yzzz = cbuffer.data(fg_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxy_zzzz = cbuffer.data(fg_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_xxz_xxxx = cbuffer.data(fg_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxz_xxxy = cbuffer.data(fg_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxz_xxxz = cbuffer.data(fg_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxz_xxyy = cbuffer.data(fg_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxz_xxyz = cbuffer.data(fg_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxz_xxzz = cbuffer.data(fg_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xxz_xyyy = cbuffer.data(fg_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xxz_xyyz = cbuffer.data(fg_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xxz_xyzz = cbuffer.data(fg_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xxz_xzzz = cbuffer.data(fg_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xxz_yyyy = cbuffer.data(fg_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xxz_yyyz = cbuffer.data(fg_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_xxz_yyzz = cbuffer.data(fg_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xxz_yzzz = cbuffer.data(fg_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xxz_zzzz = cbuffer.data(fg_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xyy_xxxx = cbuffer.data(fg_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xyy_xxxy = cbuffer.data(fg_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xyy_xxxz = cbuffer.data(fg_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_xyy_xxyy = cbuffer.data(fg_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xyy_xxyz = cbuffer.data(fg_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xyy_xxzz = cbuffer.data(fg_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xyy_xyyy = cbuffer.data(fg_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xyy_xyyz = cbuffer.data(fg_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xyy_xyzz = cbuffer.data(fg_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_xyy_xzzz = cbuffer.data(fg_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xyy_yyyy = cbuffer.data(fg_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xyy_yyyz = cbuffer.data(fg_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xyy_yyzz = cbuffer.data(fg_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xyy_yzzz = cbuffer.data(fg_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xyy_zzzz = cbuffer.data(fg_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xyz_xxxx = cbuffer.data(fg_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xyz_xxxy = cbuffer.data(fg_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xyz_xxxz = cbuffer.data(fg_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xyz_xxyy = cbuffer.data(fg_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xyz_xxyz = cbuffer.data(fg_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xyz_xxzz = cbuffer.data(fg_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xyz_xyyy = cbuffer.data(fg_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xyz_xyyz = cbuffer.data(fg_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xyz_xyzz = cbuffer.data(fg_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xyz_xzzz = cbuffer.data(fg_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xyz_yyyy = cbuffer.data(fg_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xyz_yyyz = cbuffer.data(fg_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xyz_yyzz = cbuffer.data(fg_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xyz_yzzz = cbuffer.data(fg_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xyz_zzzz = cbuffer.data(fg_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xzz_xxxx = cbuffer.data(fg_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xzz_xxxy = cbuffer.data(fg_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xzz_xxxz = cbuffer.data(fg_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xzz_xxyy = cbuffer.data(fg_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xzz_xxyz = cbuffer.data(fg_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xzz_xxzz = cbuffer.data(fg_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xzz_xyyy = cbuffer.data(fg_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xzz_xyyz = cbuffer.data(fg_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xzz_xyzz = cbuffer.data(fg_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xzz_xzzz = cbuffer.data(fg_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xzz_yyyy = cbuffer.data(fg_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xzz_yyyz = cbuffer.data(fg_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xzz_yyzz = cbuffer.data(fg_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xzz_yzzz = cbuffer.data(fg_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xzz_zzzz = cbuffer.data(fg_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_yyy_xxxx = cbuffer.data(fg_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_yyy_xxxy = cbuffer.data(fg_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_yyy_xxxz = cbuffer.data(fg_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_yyy_xxyy = cbuffer.data(fg_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_yyy_xxyz = cbuffer.data(fg_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_yyy_xxzz = cbuffer.data(fg_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_yyy_xyyy = cbuffer.data(fg_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_yyy_xyyz = cbuffer.data(fg_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_yyy_xyzz = cbuffer.data(fg_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_yyy_xzzz = cbuffer.data(fg_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_yyy_yyyy = cbuffer.data(fg_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_yyy_yyyz = cbuffer.data(fg_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_yyy_yyzz = cbuffer.data(fg_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_yyy_yzzz = cbuffer.data(fg_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_yyy_zzzz = cbuffer.data(fg_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_yyz_xxxx = cbuffer.data(fg_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_yyz_xxxy = cbuffer.data(fg_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_yyz_xxxz = cbuffer.data(fg_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_yyz_xxyy = cbuffer.data(fg_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_yyz_xxyz = cbuffer.data(fg_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_yyz_xxzz = cbuffer.data(fg_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_yyz_xyyy = cbuffer.data(fg_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_yyz_xyyz = cbuffer.data(fg_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_yyz_xyzz = cbuffer.data(fg_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_yyz_xzzz = cbuffer.data(fg_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_yyz_yyyy = cbuffer.data(fg_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_yyz_yyyz = cbuffer.data(fg_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_yyz_yyzz = cbuffer.data(fg_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_yyz_yzzz = cbuffer.data(fg_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_yyz_zzzz = cbuffer.data(fg_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_yzz_xxxx = cbuffer.data(fg_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_yzz_xxxy = cbuffer.data(fg_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_yzz_xxxz = cbuffer.data(fg_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_yzz_xxyy = cbuffer.data(fg_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_yzz_xxyz = cbuffer.data(fg_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_yzz_xxzz = cbuffer.data(fg_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_yzz_xyyy = cbuffer.data(fg_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_yzz_xyyz = cbuffer.data(fg_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_yzz_xyzz = cbuffer.data(fg_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_yzz_xzzz = cbuffer.data(fg_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_yzz_yyyy = cbuffer.data(fg_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_yzz_yyyz = cbuffer.data(fg_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_yzz_yyzz = cbuffer.data(fg_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_yzz_yzzz = cbuffer.data(fg_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_yzz_zzzz = cbuffer.data(fg_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_zzz_xxxx = cbuffer.data(fg_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_zzz_xxxy = cbuffer.data(fg_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_zzz_xxxz = cbuffer.data(fg_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_zzz_xxyy = cbuffer.data(fg_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_zzz_xxyz = cbuffer.data(fg_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_zzz_xxzz = cbuffer.data(fg_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_zzz_xyyy = cbuffer.data(fg_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_zzz_xyyz = cbuffer.data(fg_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_zzz_xyzz = cbuffer.data(fg_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_zzz_xzzz = cbuffer.data(fg_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_zzz_yyyy = cbuffer.data(fg_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_zzz_yyyz = cbuffer.data(fg_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_zzz_yyzz = cbuffer.data(fg_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_zzz_yzzz = cbuffer.data(fg_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_zzz_zzzz = cbuffer.data(fg_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_xxx_xxxx = cbuffer.data(fg_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xxx_xxxy = cbuffer.data(fg_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xxx_xxxz = cbuffer.data(fg_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xxx_xxyy = cbuffer.data(fg_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xxx_xxyz = cbuffer.data(fg_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xxx_xxzz = cbuffer.data(fg_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_xxx_xyyy = cbuffer.data(fg_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xxx_xyyz = cbuffer.data(fg_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xxx_xyzz = cbuffer.data(fg_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xxx_xzzz = cbuffer.data(fg_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_xxx_yyyy = cbuffer.data(fg_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xxx_yyyz = cbuffer.data(fg_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_xxx_yyzz = cbuffer.data(fg_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xxx_yzzz = cbuffer.data(fg_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xxx_zzzz = cbuffer.data(fg_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_xxy_xxxx = cbuffer.data(fg_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_xxy_xxxy = cbuffer.data(fg_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_xxy_xxxz = cbuffer.data(fg_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_xxy_xxyy = cbuffer.data(fg_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_xxy_xxyz = cbuffer.data(fg_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_xxy_xxzz = cbuffer.data(fg_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_xxy_xyyy = cbuffer.data(fg_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_xxy_xyyz = cbuffer.data(fg_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_xxy_xyzz = cbuffer.data(fg_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_xxy_xzzz = cbuffer.data(fg_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_xxy_yyyy = cbuffer.data(fg_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_xxy_yyyz = cbuffer.data(fg_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_xxy_yyzz = cbuffer.data(fg_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_xxy_yzzz = cbuffer.data(fg_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_xxy_zzzz = cbuffer.data(fg_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_z_xxz_xxxx = cbuffer.data(fg_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_xxz_xxxy = cbuffer.data(fg_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_xxz_xxxz = cbuffer.data(fg_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_xxz_xxyy = cbuffer.data(fg_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_xxz_xxyz = cbuffer.data(fg_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_xxz_xxzz = cbuffer.data(fg_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_xxz_xyyy = cbuffer.data(fg_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xxz_xyyz = cbuffer.data(fg_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xxz_xyzz = cbuffer.data(fg_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xxz_xzzz = cbuffer.data(fg_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xxz_yyyy = cbuffer.data(fg_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xxz_yyyz = cbuffer.data(fg_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_xxz_yyzz = cbuffer.data(fg_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_xxz_yzzz = cbuffer.data(fg_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_xxz_zzzz = cbuffer.data(fg_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_xyy_xxxx = cbuffer.data(fg_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_xyy_xxxy = cbuffer.data(fg_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_xyy_xxxz = cbuffer.data(fg_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_xyy_xxyy = cbuffer.data(fg_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_xyy_xxyz = cbuffer.data(fg_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_xyy_xxzz = cbuffer.data(fg_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_xyy_xyyy = cbuffer.data(fg_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_xyy_xyyz = cbuffer.data(fg_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_xyy_xyzz = cbuffer.data(fg_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_xyy_xzzz = cbuffer.data(fg_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_xyy_yyyy = cbuffer.data(fg_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_xyy_yyyz = cbuffer.data(fg_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_xyy_yyzz = cbuffer.data(fg_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_xyy_yzzz = cbuffer.data(fg_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_xyy_zzzz = cbuffer.data(fg_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_xyz_xxxx = cbuffer.data(fg_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_xyz_xxxy = cbuffer.data(fg_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_xyz_xxxz = cbuffer.data(fg_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_xyz_xxyy = cbuffer.data(fg_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_xyz_xxyz = cbuffer.data(fg_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_xyz_xxzz = cbuffer.data(fg_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_xyz_xyyy = cbuffer.data(fg_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_xyz_xyyz = cbuffer.data(fg_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_xyz_xyzz = cbuffer.data(fg_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_xyz_xzzz = cbuffer.data(fg_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_xyz_yyyy = cbuffer.data(fg_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_xyz_yyyz = cbuffer.data(fg_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_xyz_yyzz = cbuffer.data(fg_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_xyz_yzzz = cbuffer.data(fg_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_xyz_zzzz = cbuffer.data(fg_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_xzz_xxxx = cbuffer.data(fg_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_xzz_xxxy = cbuffer.data(fg_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_xzz_xxxz = cbuffer.data(fg_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_z_xzz_xxyy = cbuffer.data(fg_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_z_xzz_xxyz = cbuffer.data(fg_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_z_xzz_xxzz = cbuffer.data(fg_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_z_xzz_xyyy = cbuffer.data(fg_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_z_xzz_xyyz = cbuffer.data(fg_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_z_xzz_xyzz = cbuffer.data(fg_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_z_xzz_xzzz = cbuffer.data(fg_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_z_xzz_yyyy = cbuffer.data(fg_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_z_xzz_yyyz = cbuffer.data(fg_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_z_xzz_yyzz = cbuffer.data(fg_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_z_xzz_yzzz = cbuffer.data(fg_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_z_xzz_zzzz = cbuffer.data(fg_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_z_yyy_xxxx = cbuffer.data(fg_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_z_yyy_xxxy = cbuffer.data(fg_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_z_yyy_xxxz = cbuffer.data(fg_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_z_yyy_xxyy = cbuffer.data(fg_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_z_yyy_xxyz = cbuffer.data(fg_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_z_yyy_xxzz = cbuffer.data(fg_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_z_yyy_xyyy = cbuffer.data(fg_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_z_yyy_xyyz = cbuffer.data(fg_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_z_yyy_xyzz = cbuffer.data(fg_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_z_yyy_xzzz = cbuffer.data(fg_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_z_yyy_yyyy = cbuffer.data(fg_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_z_yyy_yyyz = cbuffer.data(fg_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_z_yyy_yyzz = cbuffer.data(fg_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_z_yyy_yzzz = cbuffer.data(fg_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_z_yyy_zzzz = cbuffer.data(fg_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_z_yyz_xxxx = cbuffer.data(fg_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_z_yyz_xxxy = cbuffer.data(fg_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_z_yyz_xxxz = cbuffer.data(fg_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_z_yyz_xxyy = cbuffer.data(fg_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_z_yyz_xxyz = cbuffer.data(fg_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_z_yyz_xxzz = cbuffer.data(fg_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_z_yyz_xyyy = cbuffer.data(fg_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_z_yyz_xyyz = cbuffer.data(fg_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_z_yyz_xyzz = cbuffer.data(fg_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_z_yyz_xzzz = cbuffer.data(fg_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_z_yyz_yyyy = cbuffer.data(fg_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_z_yyz_yyyz = cbuffer.data(fg_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_z_yyz_yyzz = cbuffer.data(fg_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_z_yyz_yzzz = cbuffer.data(fg_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_z_yyz_zzzz = cbuffer.data(fg_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_z_yzz_xxxx = cbuffer.data(fg_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_yzz_xxxy = cbuffer.data(fg_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_yzz_xxxz = cbuffer.data(fg_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_yzz_xxyy = cbuffer.data(fg_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_yzz_xxyz = cbuffer.data(fg_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_yzz_xxzz = cbuffer.data(fg_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_yzz_xyyy = cbuffer.data(fg_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_yzz_xyyz = cbuffer.data(fg_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_yzz_xyzz = cbuffer.data(fg_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_yzz_xzzz = cbuffer.data(fg_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_yzz_yyyy = cbuffer.data(fg_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_yzz_yyyz = cbuffer.data(fg_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_yzz_yyzz = cbuffer.data(fg_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_yzz_yzzz = cbuffer.data(fg_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_yzz_zzzz = cbuffer.data(fg_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_zzz_xxxx = cbuffer.data(fg_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_zzz_xxxy = cbuffer.data(fg_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_zzz_xxxz = cbuffer.data(fg_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_zzz_xxyy = cbuffer.data(fg_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_zzz_xxyz = cbuffer.data(fg_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_zzz_xxzz = cbuffer.data(fg_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_zzz_xyyy = cbuffer.data(fg_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_zzz_xyyz = cbuffer.data(fg_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_zzz_xyzz = cbuffer.data(fg_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_zzz_xzzz = cbuffer.data(fg_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_zzz_yyyy = cbuffer.data(fg_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_zzz_yyyz = cbuffer.data(fg_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_zzz_yyzz = cbuffer.data(fg_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_zzz_yzzz = cbuffer.data(fg_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_zzz_zzzz = cbuffer.data(fg_geom_01_off + 449 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FHSS

            const auto fh_geom_01_off = idx_geom_01_fhxx + i * dcomps + j;

            auto g_0_x_xxx_xxxxx = cbuffer.data(fh_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxy = cbuffer.data(fh_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxz = cbuffer.data(fh_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyy = cbuffer.data(fh_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyz = cbuffer.data(fh_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xxxzz = cbuffer.data(fh_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyy = cbuffer.data(fh_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyz = cbuffer.data(fh_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_xxyzz = cbuffer.data(fh_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_xxzzz = cbuffer.data(fh_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyy = cbuffer.data(fh_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyz = cbuffer.data(fh_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxx_xyyzz = cbuffer.data(fh_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxx_xyzzz = cbuffer.data(fh_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxx_xzzzz = cbuffer.data(fh_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyy = cbuffer.data(fh_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyz = cbuffer.data(fh_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxx_yyyzz = cbuffer.data(fh_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxx_yyzzz = cbuffer.data(fh_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxx_yzzzz = cbuffer.data(fh_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxx_zzzzz = cbuffer.data(fh_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxx = cbuffer.data(fh_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxy = cbuffer.data(fh_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxz = cbuffer.data(fh_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyy = cbuffer.data(fh_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyz = cbuffer.data(fh_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxy_xxxzz = cbuffer.data(fh_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyy = cbuffer.data(fh_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyz = cbuffer.data(fh_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxy_xxyzz = cbuffer.data(fh_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxy_xxzzz = cbuffer.data(fh_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyy = cbuffer.data(fh_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyz = cbuffer.data(fh_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxy_xyyzz = cbuffer.data(fh_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxy_xyzzz = cbuffer.data(fh_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxy_xzzzz = cbuffer.data(fh_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyy = cbuffer.data(fh_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyz = cbuffer.data(fh_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxy_yyyzz = cbuffer.data(fh_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxy_yyzzz = cbuffer.data(fh_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxy_yzzzz = cbuffer.data(fh_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxy_zzzzz = cbuffer.data(fh_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxx = cbuffer.data(fh_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxy = cbuffer.data(fh_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxz = cbuffer.data(fh_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyy = cbuffer.data(fh_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyz = cbuffer.data(fh_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxz_xxxzz = cbuffer.data(fh_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyy = cbuffer.data(fh_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyz = cbuffer.data(fh_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxz_xxyzz = cbuffer.data(fh_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxz_xxzzz = cbuffer.data(fh_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyy = cbuffer.data(fh_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyz = cbuffer.data(fh_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxz_xyyzz = cbuffer.data(fh_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxz_xyzzz = cbuffer.data(fh_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxz_xzzzz = cbuffer.data(fh_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyy = cbuffer.data(fh_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyz = cbuffer.data(fh_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxz_yyyzz = cbuffer.data(fh_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxz_yyzzz = cbuffer.data(fh_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxz_yzzzz = cbuffer.data(fh_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxz_zzzzz = cbuffer.data(fh_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxx = cbuffer.data(fh_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxy = cbuffer.data(fh_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxz = cbuffer.data(fh_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyy = cbuffer.data(fh_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyz = cbuffer.data(fh_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyy_xxxzz = cbuffer.data(fh_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyy = cbuffer.data(fh_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyz = cbuffer.data(fh_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyy_xxyzz = cbuffer.data(fh_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xyy_xxzzz = cbuffer.data(fh_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyy = cbuffer.data(fh_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyz = cbuffer.data(fh_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xyy_xyyzz = cbuffer.data(fh_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xyy_xyzzz = cbuffer.data(fh_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xyy_xzzzz = cbuffer.data(fh_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyy = cbuffer.data(fh_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyz = cbuffer.data(fh_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xyy_yyyzz = cbuffer.data(fh_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xyy_yyzzz = cbuffer.data(fh_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xyy_yzzzz = cbuffer.data(fh_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xyy_zzzzz = cbuffer.data(fh_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxx = cbuffer.data(fh_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxy = cbuffer.data(fh_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxz = cbuffer.data(fh_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyy = cbuffer.data(fh_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyz = cbuffer.data(fh_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xyz_xxxzz = cbuffer.data(fh_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyy = cbuffer.data(fh_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyz = cbuffer.data(fh_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xyz_xxyzz = cbuffer.data(fh_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xyz_xxzzz = cbuffer.data(fh_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyy = cbuffer.data(fh_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyz = cbuffer.data(fh_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xyz_xyyzz = cbuffer.data(fh_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xyz_xyzzz = cbuffer.data(fh_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xyz_xzzzz = cbuffer.data(fh_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyy = cbuffer.data(fh_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyz = cbuffer.data(fh_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xyz_yyyzz = cbuffer.data(fh_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xyz_yyzzz = cbuffer.data(fh_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xyz_yzzzz = cbuffer.data(fh_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xyz_zzzzz = cbuffer.data(fh_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxx = cbuffer.data(fh_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxy = cbuffer.data(fh_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxz = cbuffer.data(fh_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyy = cbuffer.data(fh_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyz = cbuffer.data(fh_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xzz_xxxzz = cbuffer.data(fh_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyy = cbuffer.data(fh_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyz = cbuffer.data(fh_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xzz_xxyzz = cbuffer.data(fh_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xzz_xxzzz = cbuffer.data(fh_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyy = cbuffer.data(fh_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyz = cbuffer.data(fh_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xzz_xyyzz = cbuffer.data(fh_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xzz_xyzzz = cbuffer.data(fh_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xzz_xzzzz = cbuffer.data(fh_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyy = cbuffer.data(fh_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyz = cbuffer.data(fh_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xzz_yyyzz = cbuffer.data(fh_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xzz_yyzzz = cbuffer.data(fh_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xzz_yzzzz = cbuffer.data(fh_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xzz_zzzzz = cbuffer.data(fh_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxx = cbuffer.data(fh_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxy = cbuffer.data(fh_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxz = cbuffer.data(fh_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyy = cbuffer.data(fh_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyz = cbuffer.data(fh_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yyy_xxxzz = cbuffer.data(fh_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyy = cbuffer.data(fh_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyz = cbuffer.data(fh_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yyy_xxyzz = cbuffer.data(fh_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yyy_xxzzz = cbuffer.data(fh_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyy = cbuffer.data(fh_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyz = cbuffer.data(fh_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_yyy_xyyzz = cbuffer.data(fh_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yyy_xyzzz = cbuffer.data(fh_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_yyy_xzzzz = cbuffer.data(fh_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyy = cbuffer.data(fh_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyz = cbuffer.data(fh_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_yyy_yyyzz = cbuffer.data(fh_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_yyy_yyzzz = cbuffer.data(fh_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_yyy_yzzzz = cbuffer.data(fh_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_yyy_zzzzz = cbuffer.data(fh_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxx = cbuffer.data(fh_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxy = cbuffer.data(fh_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxz = cbuffer.data(fh_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyy = cbuffer.data(fh_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyz = cbuffer.data(fh_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yyz_xxxzz = cbuffer.data(fh_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyy = cbuffer.data(fh_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyz = cbuffer.data(fh_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yyz_xxyzz = cbuffer.data(fh_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_yyz_xxzzz = cbuffer.data(fh_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyy = cbuffer.data(fh_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyz = cbuffer.data(fh_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yyz_xyyzz = cbuffer.data(fh_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yyz_xyzzz = cbuffer.data(fh_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yyz_xzzzz = cbuffer.data(fh_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyy = cbuffer.data(fh_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyz = cbuffer.data(fh_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_yyz_yyyzz = cbuffer.data(fh_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_yyz_yyzzz = cbuffer.data(fh_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_yyz_yzzzz = cbuffer.data(fh_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_yyz_zzzzz = cbuffer.data(fh_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxx = cbuffer.data(fh_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxy = cbuffer.data(fh_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxz = cbuffer.data(fh_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyy = cbuffer.data(fh_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyz = cbuffer.data(fh_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_yzz_xxxzz = cbuffer.data(fh_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyy = cbuffer.data(fh_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyz = cbuffer.data(fh_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_yzz_xxyzz = cbuffer.data(fh_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_yzz_xxzzz = cbuffer.data(fh_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyy = cbuffer.data(fh_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyz = cbuffer.data(fh_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_yzz_xyyzz = cbuffer.data(fh_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_yzz_xyzzz = cbuffer.data(fh_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_yzz_xzzzz = cbuffer.data(fh_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyy = cbuffer.data(fh_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyz = cbuffer.data(fh_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_yzz_yyyzz = cbuffer.data(fh_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_yzz_yyzzz = cbuffer.data(fh_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_yzz_yzzzz = cbuffer.data(fh_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_yzz_zzzzz = cbuffer.data(fh_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxx = cbuffer.data(fh_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxy = cbuffer.data(fh_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxz = cbuffer.data(fh_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyy = cbuffer.data(fh_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyz = cbuffer.data(fh_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_zzz_xxxzz = cbuffer.data(fh_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyy = cbuffer.data(fh_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyz = cbuffer.data(fh_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_zzz_xxyzz = cbuffer.data(fh_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_zzz_xxzzz = cbuffer.data(fh_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyy = cbuffer.data(fh_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyz = cbuffer.data(fh_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_zzz_xyyzz = cbuffer.data(fh_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_zzz_xyzzz = cbuffer.data(fh_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_zzz_xzzzz = cbuffer.data(fh_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyy = cbuffer.data(fh_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyz = cbuffer.data(fh_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_zzz_yyyzz = cbuffer.data(fh_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_zzz_yyzzz = cbuffer.data(fh_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_zzz_yzzzz = cbuffer.data(fh_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_zzz_zzzzz = cbuffer.data(fh_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxx = cbuffer.data(fh_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxy = cbuffer.data(fh_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxz = cbuffer.data(fh_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyy = cbuffer.data(fh_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyz = cbuffer.data(fh_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xxx_xxxzz = cbuffer.data(fh_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyy = cbuffer.data(fh_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyz = cbuffer.data(fh_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xxx_xxyzz = cbuffer.data(fh_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xxx_xxzzz = cbuffer.data(fh_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyy = cbuffer.data(fh_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyz = cbuffer.data(fh_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xxx_xyyzz = cbuffer.data(fh_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xxx_xyzzz = cbuffer.data(fh_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xxx_xzzzz = cbuffer.data(fh_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyy = cbuffer.data(fh_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyz = cbuffer.data(fh_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xxx_yyyzz = cbuffer.data(fh_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xxx_yyzzz = cbuffer.data(fh_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xxx_yzzzz = cbuffer.data(fh_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xxx_zzzzz = cbuffer.data(fh_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxx = cbuffer.data(fh_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxy = cbuffer.data(fh_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxz = cbuffer.data(fh_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyy = cbuffer.data(fh_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyz = cbuffer.data(fh_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xxy_xxxzz = cbuffer.data(fh_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyy = cbuffer.data(fh_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyz = cbuffer.data(fh_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xxy_xxyzz = cbuffer.data(fh_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_xxy_xxzzz = cbuffer.data(fh_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyy = cbuffer.data(fh_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyz = cbuffer.data(fh_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xxy_xyyzz = cbuffer.data(fh_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xxy_xyzzz = cbuffer.data(fh_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xxy_xzzzz = cbuffer.data(fh_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyy = cbuffer.data(fh_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyz = cbuffer.data(fh_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xxy_yyyzz = cbuffer.data(fh_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xxy_yyzzz = cbuffer.data(fh_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xxy_yzzzz = cbuffer.data(fh_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xxy_zzzzz = cbuffer.data(fh_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxx = cbuffer.data(fh_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxy = cbuffer.data(fh_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxz = cbuffer.data(fh_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyy = cbuffer.data(fh_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyz = cbuffer.data(fh_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xxz_xxxzz = cbuffer.data(fh_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyy = cbuffer.data(fh_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyz = cbuffer.data(fh_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xxz_xxyzz = cbuffer.data(fh_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xxz_xxzzz = cbuffer.data(fh_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyy = cbuffer.data(fh_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyz = cbuffer.data(fh_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_xxz_xyyzz = cbuffer.data(fh_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xxz_xyzzz = cbuffer.data(fh_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xxz_xzzzz = cbuffer.data(fh_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyy = cbuffer.data(fh_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyz = cbuffer.data(fh_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xxz_yyyzz = cbuffer.data(fh_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_xxz_yyzzz = cbuffer.data(fh_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xxz_yzzzz = cbuffer.data(fh_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xxz_zzzzz = cbuffer.data(fh_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxx = cbuffer.data(fh_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxy = cbuffer.data(fh_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxz = cbuffer.data(fh_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyy = cbuffer.data(fh_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyz = cbuffer.data(fh_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xyy_xxxzz = cbuffer.data(fh_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyy = cbuffer.data(fh_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyz = cbuffer.data(fh_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xyy_xxyzz = cbuffer.data(fh_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xyy_xxzzz = cbuffer.data(fh_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyy = cbuffer.data(fh_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyz = cbuffer.data(fh_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xyy_xyyzz = cbuffer.data(fh_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xyy_xyzzz = cbuffer.data(fh_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xyy_xzzzz = cbuffer.data(fh_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyy = cbuffer.data(fh_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyz = cbuffer.data(fh_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xyy_yyyzz = cbuffer.data(fh_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xyy_yyzzz = cbuffer.data(fh_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xyy_yzzzz = cbuffer.data(fh_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xyy_zzzzz = cbuffer.data(fh_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxx = cbuffer.data(fh_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxy = cbuffer.data(fh_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxz = cbuffer.data(fh_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyy = cbuffer.data(fh_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyz = cbuffer.data(fh_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xyz_xxxzz = cbuffer.data(fh_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyy = cbuffer.data(fh_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyz = cbuffer.data(fh_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xyz_xxyzz = cbuffer.data(fh_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xyz_xxzzz = cbuffer.data(fh_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyy = cbuffer.data(fh_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyz = cbuffer.data(fh_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xyz_xyyzz = cbuffer.data(fh_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xyz_xyzzz = cbuffer.data(fh_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xyz_xzzzz = cbuffer.data(fh_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyy = cbuffer.data(fh_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyz = cbuffer.data(fh_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xyz_yyyzz = cbuffer.data(fh_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xyz_yyzzz = cbuffer.data(fh_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xyz_yzzzz = cbuffer.data(fh_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xyz_zzzzz = cbuffer.data(fh_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxx = cbuffer.data(fh_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxy = cbuffer.data(fh_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxz = cbuffer.data(fh_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyy = cbuffer.data(fh_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyz = cbuffer.data(fh_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xzz_xxxzz = cbuffer.data(fh_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyy = cbuffer.data(fh_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyz = cbuffer.data(fh_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xzz_xxyzz = cbuffer.data(fh_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xzz_xxzzz = cbuffer.data(fh_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyy = cbuffer.data(fh_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyz = cbuffer.data(fh_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xzz_xyyzz = cbuffer.data(fh_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xzz_xyzzz = cbuffer.data(fh_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xzz_xzzzz = cbuffer.data(fh_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyy = cbuffer.data(fh_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyz = cbuffer.data(fh_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xzz_yyyzz = cbuffer.data(fh_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xzz_yyzzz = cbuffer.data(fh_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xzz_yzzzz = cbuffer.data(fh_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xzz_zzzzz = cbuffer.data(fh_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxx = cbuffer.data(fh_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxy = cbuffer.data(fh_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxz = cbuffer.data(fh_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyy = cbuffer.data(fh_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyz = cbuffer.data(fh_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_yyy_xxxzz = cbuffer.data(fh_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyy = cbuffer.data(fh_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyz = cbuffer.data(fh_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_yyy_xxyzz = cbuffer.data(fh_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_yyy_xxzzz = cbuffer.data(fh_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyy = cbuffer.data(fh_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyz = cbuffer.data(fh_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_yyy_xyyzz = cbuffer.data(fh_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_yyy_xyzzz = cbuffer.data(fh_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_yyy_xzzzz = cbuffer.data(fh_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyy = cbuffer.data(fh_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyz = cbuffer.data(fh_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_yyy_yyyzz = cbuffer.data(fh_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_yyy_yyzzz = cbuffer.data(fh_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_yyy_yzzzz = cbuffer.data(fh_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_yyy_zzzzz = cbuffer.data(fh_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxx = cbuffer.data(fh_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxy = cbuffer.data(fh_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxz = cbuffer.data(fh_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyy = cbuffer.data(fh_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyz = cbuffer.data(fh_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_yyz_xxxzz = cbuffer.data(fh_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyy = cbuffer.data(fh_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyz = cbuffer.data(fh_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_yyz_xxyzz = cbuffer.data(fh_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_yyz_xxzzz = cbuffer.data(fh_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyy = cbuffer.data(fh_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyz = cbuffer.data(fh_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_yyz_xyyzz = cbuffer.data(fh_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_yyz_xyzzz = cbuffer.data(fh_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_yyz_xzzzz = cbuffer.data(fh_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyy = cbuffer.data(fh_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyz = cbuffer.data(fh_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_yyz_yyyzz = cbuffer.data(fh_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_yyz_yyzzz = cbuffer.data(fh_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_yyz_yzzzz = cbuffer.data(fh_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_yyz_zzzzz = cbuffer.data(fh_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxx = cbuffer.data(fh_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxy = cbuffer.data(fh_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxz = cbuffer.data(fh_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyy = cbuffer.data(fh_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyz = cbuffer.data(fh_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_yzz_xxxzz = cbuffer.data(fh_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyy = cbuffer.data(fh_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyz = cbuffer.data(fh_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_yzz_xxyzz = cbuffer.data(fh_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_yzz_xxzzz = cbuffer.data(fh_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyy = cbuffer.data(fh_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyz = cbuffer.data(fh_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_yzz_xyyzz = cbuffer.data(fh_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_yzz_xyzzz = cbuffer.data(fh_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_yzz_xzzzz = cbuffer.data(fh_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyy = cbuffer.data(fh_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyz = cbuffer.data(fh_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_yzz_yyyzz = cbuffer.data(fh_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_yzz_yyzzz = cbuffer.data(fh_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_yzz_yzzzz = cbuffer.data(fh_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_yzz_zzzzz = cbuffer.data(fh_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxx = cbuffer.data(fh_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxy = cbuffer.data(fh_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxz = cbuffer.data(fh_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyy = cbuffer.data(fh_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyz = cbuffer.data(fh_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_zzz_xxxzz = cbuffer.data(fh_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyy = cbuffer.data(fh_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyz = cbuffer.data(fh_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_zzz_xxyzz = cbuffer.data(fh_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_zzz_xxzzz = cbuffer.data(fh_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyy = cbuffer.data(fh_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyz = cbuffer.data(fh_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_zzz_xyyzz = cbuffer.data(fh_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_zzz_xyzzz = cbuffer.data(fh_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_zzz_xzzzz = cbuffer.data(fh_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyy = cbuffer.data(fh_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyz = cbuffer.data(fh_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_zzz_yyyzz = cbuffer.data(fh_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_zzz_yyzzz = cbuffer.data(fh_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_zzz_yzzzz = cbuffer.data(fh_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_zzz_zzzzz = cbuffer.data(fh_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxx = cbuffer.data(fh_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxy = cbuffer.data(fh_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxz = cbuffer.data(fh_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyy = cbuffer.data(fh_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyz = cbuffer.data(fh_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_xxx_xxxzz = cbuffer.data(fh_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyy = cbuffer.data(fh_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyz = cbuffer.data(fh_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_xxx_xxyzz = cbuffer.data(fh_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_xxx_xxzzz = cbuffer.data(fh_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyy = cbuffer.data(fh_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyz = cbuffer.data(fh_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_xxx_xyyzz = cbuffer.data(fh_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_xxx_xyzzz = cbuffer.data(fh_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_xxx_xzzzz = cbuffer.data(fh_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyy = cbuffer.data(fh_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyz = cbuffer.data(fh_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_xxx_yyyzz = cbuffer.data(fh_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_xxx_yyzzz = cbuffer.data(fh_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_xxx_yzzzz = cbuffer.data(fh_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_xxx_zzzzz = cbuffer.data(fh_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxx = cbuffer.data(fh_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxy = cbuffer.data(fh_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxz = cbuffer.data(fh_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyy = cbuffer.data(fh_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyz = cbuffer.data(fh_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_xxy_xxxzz = cbuffer.data(fh_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyy = cbuffer.data(fh_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyz = cbuffer.data(fh_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_xxy_xxyzz = cbuffer.data(fh_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_z_xxy_xxzzz = cbuffer.data(fh_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyy = cbuffer.data(fh_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyz = cbuffer.data(fh_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xxy_xyyzz = cbuffer.data(fh_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xxy_xyzzz = cbuffer.data(fh_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xxy_xzzzz = cbuffer.data(fh_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyy = cbuffer.data(fh_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyz = cbuffer.data(fh_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xxy_yyyzz = cbuffer.data(fh_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xxy_yyzzz = cbuffer.data(fh_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xxy_yzzzz = cbuffer.data(fh_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xxy_zzzzz = cbuffer.data(fh_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxx = cbuffer.data(fh_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxy = cbuffer.data(fh_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxz = cbuffer.data(fh_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyy = cbuffer.data(fh_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyz = cbuffer.data(fh_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_xxz_xxxzz = cbuffer.data(fh_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyy = cbuffer.data(fh_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyz = cbuffer.data(fh_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_xxz_xxyzz = cbuffer.data(fh_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_xxz_xxzzz = cbuffer.data(fh_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyy = cbuffer.data(fh_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyz = cbuffer.data(fh_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_z_xxz_xyyzz = cbuffer.data(fh_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_xxz_xyzzz = cbuffer.data(fh_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_xxz_xzzzz = cbuffer.data(fh_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyy = cbuffer.data(fh_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyz = cbuffer.data(fh_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_xxz_yyyzz = cbuffer.data(fh_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_z_xxz_yyzzz = cbuffer.data(fh_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_xxz_yzzzz = cbuffer.data(fh_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_xxz_zzzzz = cbuffer.data(fh_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxx = cbuffer.data(fh_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxy = cbuffer.data(fh_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxz = cbuffer.data(fh_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyy = cbuffer.data(fh_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyz = cbuffer.data(fh_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_xyy_xxxzz = cbuffer.data(fh_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyy = cbuffer.data(fh_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyz = cbuffer.data(fh_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_xyy_xxyzz = cbuffer.data(fh_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_z_xyy_xxzzz = cbuffer.data(fh_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyy = cbuffer.data(fh_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyz = cbuffer.data(fh_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_z_xyy_xyyzz = cbuffer.data(fh_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_xyy_xyzzz = cbuffer.data(fh_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_xyy_xzzzz = cbuffer.data(fh_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyy = cbuffer.data(fh_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyz = cbuffer.data(fh_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_xyy_yyyzz = cbuffer.data(fh_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_xyy_yyzzz = cbuffer.data(fh_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_xyy_yzzzz = cbuffer.data(fh_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_xyy_zzzzz = cbuffer.data(fh_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxx = cbuffer.data(fh_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxy = cbuffer.data(fh_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxz = cbuffer.data(fh_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyy = cbuffer.data(fh_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyz = cbuffer.data(fh_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_z_xyz_xxxzz = cbuffer.data(fh_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyy = cbuffer.data(fh_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyz = cbuffer.data(fh_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_z_xyz_xxyzz = cbuffer.data(fh_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_z_xyz_xxzzz = cbuffer.data(fh_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyy = cbuffer.data(fh_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyz = cbuffer.data(fh_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_z_xyz_xyyzz = cbuffer.data(fh_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_z_xyz_xyzzz = cbuffer.data(fh_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_z_xyz_xzzzz = cbuffer.data(fh_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyy = cbuffer.data(fh_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyz = cbuffer.data(fh_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_z_xyz_yyyzz = cbuffer.data(fh_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_z_xyz_yyzzz = cbuffer.data(fh_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_z_xyz_yzzzz = cbuffer.data(fh_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_z_xyz_zzzzz = cbuffer.data(fh_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxx = cbuffer.data(fh_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxy = cbuffer.data(fh_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxz = cbuffer.data(fh_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyy = cbuffer.data(fh_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyz = cbuffer.data(fh_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_z_xzz_xxxzz = cbuffer.data(fh_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyy = cbuffer.data(fh_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyz = cbuffer.data(fh_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_z_xzz_xxyzz = cbuffer.data(fh_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_z_xzz_xxzzz = cbuffer.data(fh_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyy = cbuffer.data(fh_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyz = cbuffer.data(fh_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_z_xzz_xyyzz = cbuffer.data(fh_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_z_xzz_xyzzz = cbuffer.data(fh_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_z_xzz_xzzzz = cbuffer.data(fh_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyy = cbuffer.data(fh_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyz = cbuffer.data(fh_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_z_xzz_yyyzz = cbuffer.data(fh_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_z_xzz_yyzzz = cbuffer.data(fh_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_z_xzz_yzzzz = cbuffer.data(fh_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_z_xzz_zzzzz = cbuffer.data(fh_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxx = cbuffer.data(fh_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxy = cbuffer.data(fh_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxz = cbuffer.data(fh_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyy = cbuffer.data(fh_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyz = cbuffer.data(fh_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_z_yyy_xxxzz = cbuffer.data(fh_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyy = cbuffer.data(fh_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyz = cbuffer.data(fh_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_z_yyy_xxyzz = cbuffer.data(fh_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_z_yyy_xxzzz = cbuffer.data(fh_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyy = cbuffer.data(fh_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyz = cbuffer.data(fh_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_z_yyy_xyyzz = cbuffer.data(fh_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_z_yyy_xyzzz = cbuffer.data(fh_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_yyy_xzzzz = cbuffer.data(fh_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyy = cbuffer.data(fh_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyz = cbuffer.data(fh_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_yyy_yyyzz = cbuffer.data(fh_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_yyy_yyzzz = cbuffer.data(fh_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_yyy_yzzzz = cbuffer.data(fh_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_yyy_zzzzz = cbuffer.data(fh_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxx = cbuffer.data(fh_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxy = cbuffer.data(fh_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxz = cbuffer.data(fh_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyy = cbuffer.data(fh_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyz = cbuffer.data(fh_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_yyz_xxxzz = cbuffer.data(fh_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyy = cbuffer.data(fh_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyz = cbuffer.data(fh_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_yyz_xxyzz = cbuffer.data(fh_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_z_yyz_xxzzz = cbuffer.data(fh_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyy = cbuffer.data(fh_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyz = cbuffer.data(fh_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_yyz_xyyzz = cbuffer.data(fh_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_yyz_xyzzz = cbuffer.data(fh_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_yyz_xzzzz = cbuffer.data(fh_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyy = cbuffer.data(fh_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyz = cbuffer.data(fh_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_yyz_yyyzz = cbuffer.data(fh_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_yyz_yyzzz = cbuffer.data(fh_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_yyz_yzzzz = cbuffer.data(fh_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_yyz_zzzzz = cbuffer.data(fh_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxx = cbuffer.data(fh_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxy = cbuffer.data(fh_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxz = cbuffer.data(fh_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyy = cbuffer.data(fh_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyz = cbuffer.data(fh_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_yzz_xxxzz = cbuffer.data(fh_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyy = cbuffer.data(fh_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyz = cbuffer.data(fh_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_yzz_xxyzz = cbuffer.data(fh_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_yzz_xxzzz = cbuffer.data(fh_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyy = cbuffer.data(fh_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyz = cbuffer.data(fh_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_z_yzz_xyyzz = cbuffer.data(fh_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_yzz_xyzzz = cbuffer.data(fh_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_yzz_xzzzz = cbuffer.data(fh_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyy = cbuffer.data(fh_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyz = cbuffer.data(fh_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_yzz_yyyzz = cbuffer.data(fh_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_yzz_yyzzz = cbuffer.data(fh_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_yzz_yzzzz = cbuffer.data(fh_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_yzz_zzzzz = cbuffer.data(fh_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxx = cbuffer.data(fh_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxy = cbuffer.data(fh_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxz = cbuffer.data(fh_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyy = cbuffer.data(fh_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyz = cbuffer.data(fh_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_zzz_xxxzz = cbuffer.data(fh_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyy = cbuffer.data(fh_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyz = cbuffer.data(fh_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_zzz_xxyzz = cbuffer.data(fh_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_zzz_xxzzz = cbuffer.data(fh_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyy = cbuffer.data(fh_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyz = cbuffer.data(fh_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_zzz_xyyzz = cbuffer.data(fh_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_zzz_xyzzz = cbuffer.data(fh_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_zzz_xzzzz = cbuffer.data(fh_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyy = cbuffer.data(fh_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyz = cbuffer.data(fh_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_zzz_yyyzz = cbuffer.data(fh_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_zzz_yyzzz = cbuffer.data(fh_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_zzz_yzzzz = cbuffer.data(fh_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_zzz_zzzzz = cbuffer.data(fh_geom_01_off + 629 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ggxx

            const auto gg_geom_01_off = idx_geom_01_ggxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxx_xxxx = cbuffer.data(gg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxy = cbuffer.data(gg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxz = cbuffer.data(gg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyy = cbuffer.data(gg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyz = cbuffer.data(gg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_xxzz = cbuffer.data(gg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyy = cbuffer.data(gg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyz = cbuffer.data(gg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxx_xyzz = cbuffer.data(gg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxx_xzzz = cbuffer.data(gg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyy = cbuffer.data(gg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyz = cbuffer.data(gg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxx_yyzz = cbuffer.data(gg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxx_yzzz = cbuffer.data(gg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxx_zzzz = cbuffer.data(gg_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xxxx, g_0_x_xxx_xxxxx, g_0_x_xxx_xxxxy, g_0_x_xxx_xxxxz, g_0_x_xxx_xxxy, g_0_x_xxx_xxxyy, g_0_x_xxx_xxxyz, g_0_x_xxx_xxxz, g_0_x_xxx_xxxzz, g_0_x_xxx_xxyy, g_0_x_xxx_xxyyy, g_0_x_xxx_xxyyz, g_0_x_xxx_xxyz, g_0_x_xxx_xxyzz, g_0_x_xxx_xxzz, g_0_x_xxx_xxzzz, g_0_x_xxx_xyyy, g_0_x_xxx_xyyyy, g_0_x_xxx_xyyyz, g_0_x_xxx_xyyz, g_0_x_xxx_xyyzz, g_0_x_xxx_xyzz, g_0_x_xxx_xyzzz, g_0_x_xxx_xzzz, g_0_x_xxx_xzzzz, g_0_x_xxx_yyyy, g_0_x_xxx_yyyz, g_0_x_xxx_yyzz, g_0_x_xxx_yzzz, g_0_x_xxx_zzzz, g_0_x_xxxx_xxxx, g_0_x_xxxx_xxxy, g_0_x_xxxx_xxxz, g_0_x_xxxx_xxyy, g_0_x_xxxx_xxyz, g_0_x_xxxx_xxzz, g_0_x_xxxx_xyyy, g_0_x_xxxx_xyyz, g_0_x_xxxx_xyzz, g_0_x_xxxx_xzzz, g_0_x_xxxx_yyyy, g_0_x_xxxx_yyyz, g_0_x_xxxx_yyzz, g_0_x_xxxx_yzzz, g_0_x_xxxx_zzzz, g_xxx_xxxx, g_xxx_xxxy, g_xxx_xxxz, g_xxx_xxyy, g_xxx_xxyz, g_xxx_xxzz, g_xxx_xyyy, g_xxx_xyyz, g_xxx_xyzz, g_xxx_xzzz, g_xxx_yyyy, g_xxx_yyyz, g_xxx_yyzz, g_xxx_yzzz, g_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxx_xxxx[k] = g_xxx_xxxx[k] - g_0_x_xxx_xxxx[k] * ab_x + g_0_x_xxx_xxxxx[k];

                g_0_x_xxxx_xxxy[k] = g_xxx_xxxy[k] - g_0_x_xxx_xxxy[k] * ab_x + g_0_x_xxx_xxxxy[k];

                g_0_x_xxxx_xxxz[k] = g_xxx_xxxz[k] - g_0_x_xxx_xxxz[k] * ab_x + g_0_x_xxx_xxxxz[k];

                g_0_x_xxxx_xxyy[k] = g_xxx_xxyy[k] - g_0_x_xxx_xxyy[k] * ab_x + g_0_x_xxx_xxxyy[k];

                g_0_x_xxxx_xxyz[k] = g_xxx_xxyz[k] - g_0_x_xxx_xxyz[k] * ab_x + g_0_x_xxx_xxxyz[k];

                g_0_x_xxxx_xxzz[k] = g_xxx_xxzz[k] - g_0_x_xxx_xxzz[k] * ab_x + g_0_x_xxx_xxxzz[k];

                g_0_x_xxxx_xyyy[k] = g_xxx_xyyy[k] - g_0_x_xxx_xyyy[k] * ab_x + g_0_x_xxx_xxyyy[k];

                g_0_x_xxxx_xyyz[k] = g_xxx_xyyz[k] - g_0_x_xxx_xyyz[k] * ab_x + g_0_x_xxx_xxyyz[k];

                g_0_x_xxxx_xyzz[k] = g_xxx_xyzz[k] - g_0_x_xxx_xyzz[k] * ab_x + g_0_x_xxx_xxyzz[k];

                g_0_x_xxxx_xzzz[k] = g_xxx_xzzz[k] - g_0_x_xxx_xzzz[k] * ab_x + g_0_x_xxx_xxzzz[k];

                g_0_x_xxxx_yyyy[k] = g_xxx_yyyy[k] - g_0_x_xxx_yyyy[k] * ab_x + g_0_x_xxx_xyyyy[k];

                g_0_x_xxxx_yyyz[k] = g_xxx_yyyz[k] - g_0_x_xxx_yyyz[k] * ab_x + g_0_x_xxx_xyyyz[k];

                g_0_x_xxxx_yyzz[k] = g_xxx_yyzz[k] - g_0_x_xxx_yyzz[k] * ab_x + g_0_x_xxx_xyyzz[k];

                g_0_x_xxxx_yzzz[k] = g_xxx_yzzz[k] - g_0_x_xxx_yzzz[k] * ab_x + g_0_x_xxx_xyzzz[k];

                g_0_x_xxxx_zzzz[k] = g_xxx_zzzz[k] - g_0_x_xxx_zzzz[k] * ab_x + g_0_x_xxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxy_xxxx = cbuffer.data(gg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxy = cbuffer.data(gg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxz = cbuffer.data(gg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyy = cbuffer.data(gg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyz = cbuffer.data(gg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxy_xxzz = cbuffer.data(gg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyy = cbuffer.data(gg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyz = cbuffer.data(gg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxy_xyzz = cbuffer.data(gg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxy_xzzz = cbuffer.data(gg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyy = cbuffer.data(gg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyz = cbuffer.data(gg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxy_yyzz = cbuffer.data(gg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxy_yzzz = cbuffer.data(gg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxy_zzzz = cbuffer.data(gg_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xxxx, g_0_x_xxx_xxxxy, g_0_x_xxx_xxxy, g_0_x_xxx_xxxyy, g_0_x_xxx_xxxyz, g_0_x_xxx_xxxz, g_0_x_xxx_xxyy, g_0_x_xxx_xxyyy, g_0_x_xxx_xxyyz, g_0_x_xxx_xxyz, g_0_x_xxx_xxyzz, g_0_x_xxx_xxzz, g_0_x_xxx_xyyy, g_0_x_xxx_xyyyy, g_0_x_xxx_xyyyz, g_0_x_xxx_xyyz, g_0_x_xxx_xyyzz, g_0_x_xxx_xyzz, g_0_x_xxx_xyzzz, g_0_x_xxx_xzzz, g_0_x_xxx_yyyy, g_0_x_xxx_yyyyy, g_0_x_xxx_yyyyz, g_0_x_xxx_yyyz, g_0_x_xxx_yyyzz, g_0_x_xxx_yyzz, g_0_x_xxx_yyzzz, g_0_x_xxx_yzzz, g_0_x_xxx_yzzzz, g_0_x_xxx_zzzz, g_0_x_xxxy_xxxx, g_0_x_xxxy_xxxy, g_0_x_xxxy_xxxz, g_0_x_xxxy_xxyy, g_0_x_xxxy_xxyz, g_0_x_xxxy_xxzz, g_0_x_xxxy_xyyy, g_0_x_xxxy_xyyz, g_0_x_xxxy_xyzz, g_0_x_xxxy_xzzz, g_0_x_xxxy_yyyy, g_0_x_xxxy_yyyz, g_0_x_xxxy_yyzz, g_0_x_xxxy_yzzz, g_0_x_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxy_xxxx[k] = -g_0_x_xxx_xxxx[k] * ab_y + g_0_x_xxx_xxxxy[k];

                g_0_x_xxxy_xxxy[k] = -g_0_x_xxx_xxxy[k] * ab_y + g_0_x_xxx_xxxyy[k];

                g_0_x_xxxy_xxxz[k] = -g_0_x_xxx_xxxz[k] * ab_y + g_0_x_xxx_xxxyz[k];

                g_0_x_xxxy_xxyy[k] = -g_0_x_xxx_xxyy[k] * ab_y + g_0_x_xxx_xxyyy[k];

                g_0_x_xxxy_xxyz[k] = -g_0_x_xxx_xxyz[k] * ab_y + g_0_x_xxx_xxyyz[k];

                g_0_x_xxxy_xxzz[k] = -g_0_x_xxx_xxzz[k] * ab_y + g_0_x_xxx_xxyzz[k];

                g_0_x_xxxy_xyyy[k] = -g_0_x_xxx_xyyy[k] * ab_y + g_0_x_xxx_xyyyy[k];

                g_0_x_xxxy_xyyz[k] = -g_0_x_xxx_xyyz[k] * ab_y + g_0_x_xxx_xyyyz[k];

                g_0_x_xxxy_xyzz[k] = -g_0_x_xxx_xyzz[k] * ab_y + g_0_x_xxx_xyyzz[k];

                g_0_x_xxxy_xzzz[k] = -g_0_x_xxx_xzzz[k] * ab_y + g_0_x_xxx_xyzzz[k];

                g_0_x_xxxy_yyyy[k] = -g_0_x_xxx_yyyy[k] * ab_y + g_0_x_xxx_yyyyy[k];

                g_0_x_xxxy_yyyz[k] = -g_0_x_xxx_yyyz[k] * ab_y + g_0_x_xxx_yyyyz[k];

                g_0_x_xxxy_yyzz[k] = -g_0_x_xxx_yyzz[k] * ab_y + g_0_x_xxx_yyyzz[k];

                g_0_x_xxxy_yzzz[k] = -g_0_x_xxx_yzzz[k] * ab_y + g_0_x_xxx_yyzzz[k];

                g_0_x_xxxy_zzzz[k] = -g_0_x_xxx_zzzz[k] * ab_y + g_0_x_xxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxz_xxxx = cbuffer.data(gg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxy = cbuffer.data(gg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxz = cbuffer.data(gg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyy = cbuffer.data(gg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyz = cbuffer.data(gg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxz_xxzz = cbuffer.data(gg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyy = cbuffer.data(gg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyz = cbuffer.data(gg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxz_xyzz = cbuffer.data(gg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxz_xzzz = cbuffer.data(gg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyy = cbuffer.data(gg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyz = cbuffer.data(gg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxz_yyzz = cbuffer.data(gg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxz_yzzz = cbuffer.data(gg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxz_zzzz = cbuffer.data(gg_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xxxx, g_0_x_xxx_xxxxz, g_0_x_xxx_xxxy, g_0_x_xxx_xxxyz, g_0_x_xxx_xxxz, g_0_x_xxx_xxxzz, g_0_x_xxx_xxyy, g_0_x_xxx_xxyyz, g_0_x_xxx_xxyz, g_0_x_xxx_xxyzz, g_0_x_xxx_xxzz, g_0_x_xxx_xxzzz, g_0_x_xxx_xyyy, g_0_x_xxx_xyyyz, g_0_x_xxx_xyyz, g_0_x_xxx_xyyzz, g_0_x_xxx_xyzz, g_0_x_xxx_xyzzz, g_0_x_xxx_xzzz, g_0_x_xxx_xzzzz, g_0_x_xxx_yyyy, g_0_x_xxx_yyyyz, g_0_x_xxx_yyyz, g_0_x_xxx_yyyzz, g_0_x_xxx_yyzz, g_0_x_xxx_yyzzz, g_0_x_xxx_yzzz, g_0_x_xxx_yzzzz, g_0_x_xxx_zzzz, g_0_x_xxx_zzzzz, g_0_x_xxxz_xxxx, g_0_x_xxxz_xxxy, g_0_x_xxxz_xxxz, g_0_x_xxxz_xxyy, g_0_x_xxxz_xxyz, g_0_x_xxxz_xxzz, g_0_x_xxxz_xyyy, g_0_x_xxxz_xyyz, g_0_x_xxxz_xyzz, g_0_x_xxxz_xzzz, g_0_x_xxxz_yyyy, g_0_x_xxxz_yyyz, g_0_x_xxxz_yyzz, g_0_x_xxxz_yzzz, g_0_x_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxz_xxxx[k] = -g_0_x_xxx_xxxx[k] * ab_z + g_0_x_xxx_xxxxz[k];

                g_0_x_xxxz_xxxy[k] = -g_0_x_xxx_xxxy[k] * ab_z + g_0_x_xxx_xxxyz[k];

                g_0_x_xxxz_xxxz[k] = -g_0_x_xxx_xxxz[k] * ab_z + g_0_x_xxx_xxxzz[k];

                g_0_x_xxxz_xxyy[k] = -g_0_x_xxx_xxyy[k] * ab_z + g_0_x_xxx_xxyyz[k];

                g_0_x_xxxz_xxyz[k] = -g_0_x_xxx_xxyz[k] * ab_z + g_0_x_xxx_xxyzz[k];

                g_0_x_xxxz_xxzz[k] = -g_0_x_xxx_xxzz[k] * ab_z + g_0_x_xxx_xxzzz[k];

                g_0_x_xxxz_xyyy[k] = -g_0_x_xxx_xyyy[k] * ab_z + g_0_x_xxx_xyyyz[k];

                g_0_x_xxxz_xyyz[k] = -g_0_x_xxx_xyyz[k] * ab_z + g_0_x_xxx_xyyzz[k];

                g_0_x_xxxz_xyzz[k] = -g_0_x_xxx_xyzz[k] * ab_z + g_0_x_xxx_xyzzz[k];

                g_0_x_xxxz_xzzz[k] = -g_0_x_xxx_xzzz[k] * ab_z + g_0_x_xxx_xzzzz[k];

                g_0_x_xxxz_yyyy[k] = -g_0_x_xxx_yyyy[k] * ab_z + g_0_x_xxx_yyyyz[k];

                g_0_x_xxxz_yyyz[k] = -g_0_x_xxx_yyyz[k] * ab_z + g_0_x_xxx_yyyzz[k];

                g_0_x_xxxz_yyzz[k] = -g_0_x_xxx_yyzz[k] * ab_z + g_0_x_xxx_yyzzz[k];

                g_0_x_xxxz_yzzz[k] = -g_0_x_xxx_yzzz[k] * ab_z + g_0_x_xxx_yzzzz[k];

                g_0_x_xxxz_zzzz[k] = -g_0_x_xxx_zzzz[k] * ab_z + g_0_x_xxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyy_xxxx = cbuffer.data(gg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxy = cbuffer.data(gg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxz = cbuffer.data(gg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyy = cbuffer.data(gg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyz = cbuffer.data(gg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxyy_xxzz = cbuffer.data(gg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyy = cbuffer.data(gg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyz = cbuffer.data(gg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxyy_xyzz = cbuffer.data(gg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxyy_xzzz = cbuffer.data(gg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyy = cbuffer.data(gg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyz = cbuffer.data(gg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxyy_yyzz = cbuffer.data(gg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxyy_yzzz = cbuffer.data(gg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxyy_zzzz = cbuffer.data(gg_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxy_xxxx, g_0_x_xxy_xxxxy, g_0_x_xxy_xxxy, g_0_x_xxy_xxxyy, g_0_x_xxy_xxxyz, g_0_x_xxy_xxxz, g_0_x_xxy_xxyy, g_0_x_xxy_xxyyy, g_0_x_xxy_xxyyz, g_0_x_xxy_xxyz, g_0_x_xxy_xxyzz, g_0_x_xxy_xxzz, g_0_x_xxy_xyyy, g_0_x_xxy_xyyyy, g_0_x_xxy_xyyyz, g_0_x_xxy_xyyz, g_0_x_xxy_xyyzz, g_0_x_xxy_xyzz, g_0_x_xxy_xyzzz, g_0_x_xxy_xzzz, g_0_x_xxy_yyyy, g_0_x_xxy_yyyyy, g_0_x_xxy_yyyyz, g_0_x_xxy_yyyz, g_0_x_xxy_yyyzz, g_0_x_xxy_yyzz, g_0_x_xxy_yyzzz, g_0_x_xxy_yzzz, g_0_x_xxy_yzzzz, g_0_x_xxy_zzzz, g_0_x_xxyy_xxxx, g_0_x_xxyy_xxxy, g_0_x_xxyy_xxxz, g_0_x_xxyy_xxyy, g_0_x_xxyy_xxyz, g_0_x_xxyy_xxzz, g_0_x_xxyy_xyyy, g_0_x_xxyy_xyyz, g_0_x_xxyy_xyzz, g_0_x_xxyy_xzzz, g_0_x_xxyy_yyyy, g_0_x_xxyy_yyyz, g_0_x_xxyy_yyzz, g_0_x_xxyy_yzzz, g_0_x_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyy_xxxx[k] = -g_0_x_xxy_xxxx[k] * ab_y + g_0_x_xxy_xxxxy[k];

                g_0_x_xxyy_xxxy[k] = -g_0_x_xxy_xxxy[k] * ab_y + g_0_x_xxy_xxxyy[k];

                g_0_x_xxyy_xxxz[k] = -g_0_x_xxy_xxxz[k] * ab_y + g_0_x_xxy_xxxyz[k];

                g_0_x_xxyy_xxyy[k] = -g_0_x_xxy_xxyy[k] * ab_y + g_0_x_xxy_xxyyy[k];

                g_0_x_xxyy_xxyz[k] = -g_0_x_xxy_xxyz[k] * ab_y + g_0_x_xxy_xxyyz[k];

                g_0_x_xxyy_xxzz[k] = -g_0_x_xxy_xxzz[k] * ab_y + g_0_x_xxy_xxyzz[k];

                g_0_x_xxyy_xyyy[k] = -g_0_x_xxy_xyyy[k] * ab_y + g_0_x_xxy_xyyyy[k];

                g_0_x_xxyy_xyyz[k] = -g_0_x_xxy_xyyz[k] * ab_y + g_0_x_xxy_xyyyz[k];

                g_0_x_xxyy_xyzz[k] = -g_0_x_xxy_xyzz[k] * ab_y + g_0_x_xxy_xyyzz[k];

                g_0_x_xxyy_xzzz[k] = -g_0_x_xxy_xzzz[k] * ab_y + g_0_x_xxy_xyzzz[k];

                g_0_x_xxyy_yyyy[k] = -g_0_x_xxy_yyyy[k] * ab_y + g_0_x_xxy_yyyyy[k];

                g_0_x_xxyy_yyyz[k] = -g_0_x_xxy_yyyz[k] * ab_y + g_0_x_xxy_yyyyz[k];

                g_0_x_xxyy_yyzz[k] = -g_0_x_xxy_yyzz[k] * ab_y + g_0_x_xxy_yyyzz[k];

                g_0_x_xxyy_yzzz[k] = -g_0_x_xxy_yzzz[k] * ab_y + g_0_x_xxy_yyzzz[k];

                g_0_x_xxyy_zzzz[k] = -g_0_x_xxy_zzzz[k] * ab_y + g_0_x_xxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyz_xxxx = cbuffer.data(gg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxy = cbuffer.data(gg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxz = cbuffer.data(gg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyy = cbuffer.data(gg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyz = cbuffer.data(gg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxyz_xxzz = cbuffer.data(gg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyy = cbuffer.data(gg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyz = cbuffer.data(gg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxyz_xyzz = cbuffer.data(gg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxyz_xzzz = cbuffer.data(gg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyy = cbuffer.data(gg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyz = cbuffer.data(gg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxyz_yyzz = cbuffer.data(gg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxyz_yzzz = cbuffer.data(gg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxyz_zzzz = cbuffer.data(gg_geom_01_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyz_xxxx, g_0_x_xxyz_xxxy, g_0_x_xxyz_xxxz, g_0_x_xxyz_xxyy, g_0_x_xxyz_xxyz, g_0_x_xxyz_xxzz, g_0_x_xxyz_xyyy, g_0_x_xxyz_xyyz, g_0_x_xxyz_xyzz, g_0_x_xxyz_xzzz, g_0_x_xxyz_yyyy, g_0_x_xxyz_yyyz, g_0_x_xxyz_yyzz, g_0_x_xxyz_yzzz, g_0_x_xxyz_zzzz, g_0_x_xxz_xxxx, g_0_x_xxz_xxxxy, g_0_x_xxz_xxxy, g_0_x_xxz_xxxyy, g_0_x_xxz_xxxyz, g_0_x_xxz_xxxz, g_0_x_xxz_xxyy, g_0_x_xxz_xxyyy, g_0_x_xxz_xxyyz, g_0_x_xxz_xxyz, g_0_x_xxz_xxyzz, g_0_x_xxz_xxzz, g_0_x_xxz_xyyy, g_0_x_xxz_xyyyy, g_0_x_xxz_xyyyz, g_0_x_xxz_xyyz, g_0_x_xxz_xyyzz, g_0_x_xxz_xyzz, g_0_x_xxz_xyzzz, g_0_x_xxz_xzzz, g_0_x_xxz_yyyy, g_0_x_xxz_yyyyy, g_0_x_xxz_yyyyz, g_0_x_xxz_yyyz, g_0_x_xxz_yyyzz, g_0_x_xxz_yyzz, g_0_x_xxz_yyzzz, g_0_x_xxz_yzzz, g_0_x_xxz_yzzzz, g_0_x_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyz_xxxx[k] = -g_0_x_xxz_xxxx[k] * ab_y + g_0_x_xxz_xxxxy[k];

                g_0_x_xxyz_xxxy[k] = -g_0_x_xxz_xxxy[k] * ab_y + g_0_x_xxz_xxxyy[k];

                g_0_x_xxyz_xxxz[k] = -g_0_x_xxz_xxxz[k] * ab_y + g_0_x_xxz_xxxyz[k];

                g_0_x_xxyz_xxyy[k] = -g_0_x_xxz_xxyy[k] * ab_y + g_0_x_xxz_xxyyy[k];

                g_0_x_xxyz_xxyz[k] = -g_0_x_xxz_xxyz[k] * ab_y + g_0_x_xxz_xxyyz[k];

                g_0_x_xxyz_xxzz[k] = -g_0_x_xxz_xxzz[k] * ab_y + g_0_x_xxz_xxyzz[k];

                g_0_x_xxyz_xyyy[k] = -g_0_x_xxz_xyyy[k] * ab_y + g_0_x_xxz_xyyyy[k];

                g_0_x_xxyz_xyyz[k] = -g_0_x_xxz_xyyz[k] * ab_y + g_0_x_xxz_xyyyz[k];

                g_0_x_xxyz_xyzz[k] = -g_0_x_xxz_xyzz[k] * ab_y + g_0_x_xxz_xyyzz[k];

                g_0_x_xxyz_xzzz[k] = -g_0_x_xxz_xzzz[k] * ab_y + g_0_x_xxz_xyzzz[k];

                g_0_x_xxyz_yyyy[k] = -g_0_x_xxz_yyyy[k] * ab_y + g_0_x_xxz_yyyyy[k];

                g_0_x_xxyz_yyyz[k] = -g_0_x_xxz_yyyz[k] * ab_y + g_0_x_xxz_yyyyz[k];

                g_0_x_xxyz_yyzz[k] = -g_0_x_xxz_yyzz[k] * ab_y + g_0_x_xxz_yyyzz[k];

                g_0_x_xxyz_yzzz[k] = -g_0_x_xxz_yzzz[k] * ab_y + g_0_x_xxz_yyzzz[k];

                g_0_x_xxyz_zzzz[k] = -g_0_x_xxz_zzzz[k] * ab_y + g_0_x_xxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzz_xxxx = cbuffer.data(gg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxy = cbuffer.data(gg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxz = cbuffer.data(gg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyy = cbuffer.data(gg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyz = cbuffer.data(gg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxzz_xxzz = cbuffer.data(gg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyy = cbuffer.data(gg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyz = cbuffer.data(gg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxzz_xyzz = cbuffer.data(gg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxzz_xzzz = cbuffer.data(gg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyy = cbuffer.data(gg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyz = cbuffer.data(gg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxzz_yyzz = cbuffer.data(gg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxzz_yzzz = cbuffer.data(gg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxzz_zzzz = cbuffer.data(gg_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxz_xxxx, g_0_x_xxz_xxxxz, g_0_x_xxz_xxxy, g_0_x_xxz_xxxyz, g_0_x_xxz_xxxz, g_0_x_xxz_xxxzz, g_0_x_xxz_xxyy, g_0_x_xxz_xxyyz, g_0_x_xxz_xxyz, g_0_x_xxz_xxyzz, g_0_x_xxz_xxzz, g_0_x_xxz_xxzzz, g_0_x_xxz_xyyy, g_0_x_xxz_xyyyz, g_0_x_xxz_xyyz, g_0_x_xxz_xyyzz, g_0_x_xxz_xyzz, g_0_x_xxz_xyzzz, g_0_x_xxz_xzzz, g_0_x_xxz_xzzzz, g_0_x_xxz_yyyy, g_0_x_xxz_yyyyz, g_0_x_xxz_yyyz, g_0_x_xxz_yyyzz, g_0_x_xxz_yyzz, g_0_x_xxz_yyzzz, g_0_x_xxz_yzzz, g_0_x_xxz_yzzzz, g_0_x_xxz_zzzz, g_0_x_xxz_zzzzz, g_0_x_xxzz_xxxx, g_0_x_xxzz_xxxy, g_0_x_xxzz_xxxz, g_0_x_xxzz_xxyy, g_0_x_xxzz_xxyz, g_0_x_xxzz_xxzz, g_0_x_xxzz_xyyy, g_0_x_xxzz_xyyz, g_0_x_xxzz_xyzz, g_0_x_xxzz_xzzz, g_0_x_xxzz_yyyy, g_0_x_xxzz_yyyz, g_0_x_xxzz_yyzz, g_0_x_xxzz_yzzz, g_0_x_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzz_xxxx[k] = -g_0_x_xxz_xxxx[k] * ab_z + g_0_x_xxz_xxxxz[k];

                g_0_x_xxzz_xxxy[k] = -g_0_x_xxz_xxxy[k] * ab_z + g_0_x_xxz_xxxyz[k];

                g_0_x_xxzz_xxxz[k] = -g_0_x_xxz_xxxz[k] * ab_z + g_0_x_xxz_xxxzz[k];

                g_0_x_xxzz_xxyy[k] = -g_0_x_xxz_xxyy[k] * ab_z + g_0_x_xxz_xxyyz[k];

                g_0_x_xxzz_xxyz[k] = -g_0_x_xxz_xxyz[k] * ab_z + g_0_x_xxz_xxyzz[k];

                g_0_x_xxzz_xxzz[k] = -g_0_x_xxz_xxzz[k] * ab_z + g_0_x_xxz_xxzzz[k];

                g_0_x_xxzz_xyyy[k] = -g_0_x_xxz_xyyy[k] * ab_z + g_0_x_xxz_xyyyz[k];

                g_0_x_xxzz_xyyz[k] = -g_0_x_xxz_xyyz[k] * ab_z + g_0_x_xxz_xyyzz[k];

                g_0_x_xxzz_xyzz[k] = -g_0_x_xxz_xyzz[k] * ab_z + g_0_x_xxz_xyzzz[k];

                g_0_x_xxzz_xzzz[k] = -g_0_x_xxz_xzzz[k] * ab_z + g_0_x_xxz_xzzzz[k];

                g_0_x_xxzz_yyyy[k] = -g_0_x_xxz_yyyy[k] * ab_z + g_0_x_xxz_yyyyz[k];

                g_0_x_xxzz_yyyz[k] = -g_0_x_xxz_yyyz[k] * ab_z + g_0_x_xxz_yyyzz[k];

                g_0_x_xxzz_yyzz[k] = -g_0_x_xxz_yyzz[k] * ab_z + g_0_x_xxz_yyzzz[k];

                g_0_x_xxzz_yzzz[k] = -g_0_x_xxz_yzzz[k] * ab_z + g_0_x_xxz_yzzzz[k];

                g_0_x_xxzz_zzzz[k] = -g_0_x_xxz_zzzz[k] * ab_z + g_0_x_xxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyy_xxxx = cbuffer.data(gg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxy = cbuffer.data(gg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxz = cbuffer.data(gg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyy = cbuffer.data(gg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyz = cbuffer.data(gg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xyyy_xxzz = cbuffer.data(gg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyy = cbuffer.data(gg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyz = cbuffer.data(gg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xyyy_xyzz = cbuffer.data(gg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xyyy_xzzz = cbuffer.data(gg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyy = cbuffer.data(gg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyz = cbuffer.data(gg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xyyy_yyzz = cbuffer.data(gg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xyyy_yzzz = cbuffer.data(gg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xyyy_zzzz = cbuffer.data(gg_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyy_xxxx, g_0_x_xyy_xxxxy, g_0_x_xyy_xxxy, g_0_x_xyy_xxxyy, g_0_x_xyy_xxxyz, g_0_x_xyy_xxxz, g_0_x_xyy_xxyy, g_0_x_xyy_xxyyy, g_0_x_xyy_xxyyz, g_0_x_xyy_xxyz, g_0_x_xyy_xxyzz, g_0_x_xyy_xxzz, g_0_x_xyy_xyyy, g_0_x_xyy_xyyyy, g_0_x_xyy_xyyyz, g_0_x_xyy_xyyz, g_0_x_xyy_xyyzz, g_0_x_xyy_xyzz, g_0_x_xyy_xyzzz, g_0_x_xyy_xzzz, g_0_x_xyy_yyyy, g_0_x_xyy_yyyyy, g_0_x_xyy_yyyyz, g_0_x_xyy_yyyz, g_0_x_xyy_yyyzz, g_0_x_xyy_yyzz, g_0_x_xyy_yyzzz, g_0_x_xyy_yzzz, g_0_x_xyy_yzzzz, g_0_x_xyy_zzzz, g_0_x_xyyy_xxxx, g_0_x_xyyy_xxxy, g_0_x_xyyy_xxxz, g_0_x_xyyy_xxyy, g_0_x_xyyy_xxyz, g_0_x_xyyy_xxzz, g_0_x_xyyy_xyyy, g_0_x_xyyy_xyyz, g_0_x_xyyy_xyzz, g_0_x_xyyy_xzzz, g_0_x_xyyy_yyyy, g_0_x_xyyy_yyyz, g_0_x_xyyy_yyzz, g_0_x_xyyy_yzzz, g_0_x_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyy_xxxx[k] = -g_0_x_xyy_xxxx[k] * ab_y + g_0_x_xyy_xxxxy[k];

                g_0_x_xyyy_xxxy[k] = -g_0_x_xyy_xxxy[k] * ab_y + g_0_x_xyy_xxxyy[k];

                g_0_x_xyyy_xxxz[k] = -g_0_x_xyy_xxxz[k] * ab_y + g_0_x_xyy_xxxyz[k];

                g_0_x_xyyy_xxyy[k] = -g_0_x_xyy_xxyy[k] * ab_y + g_0_x_xyy_xxyyy[k];

                g_0_x_xyyy_xxyz[k] = -g_0_x_xyy_xxyz[k] * ab_y + g_0_x_xyy_xxyyz[k];

                g_0_x_xyyy_xxzz[k] = -g_0_x_xyy_xxzz[k] * ab_y + g_0_x_xyy_xxyzz[k];

                g_0_x_xyyy_xyyy[k] = -g_0_x_xyy_xyyy[k] * ab_y + g_0_x_xyy_xyyyy[k];

                g_0_x_xyyy_xyyz[k] = -g_0_x_xyy_xyyz[k] * ab_y + g_0_x_xyy_xyyyz[k];

                g_0_x_xyyy_xyzz[k] = -g_0_x_xyy_xyzz[k] * ab_y + g_0_x_xyy_xyyzz[k];

                g_0_x_xyyy_xzzz[k] = -g_0_x_xyy_xzzz[k] * ab_y + g_0_x_xyy_xyzzz[k];

                g_0_x_xyyy_yyyy[k] = -g_0_x_xyy_yyyy[k] * ab_y + g_0_x_xyy_yyyyy[k];

                g_0_x_xyyy_yyyz[k] = -g_0_x_xyy_yyyz[k] * ab_y + g_0_x_xyy_yyyyz[k];

                g_0_x_xyyy_yyzz[k] = -g_0_x_xyy_yyzz[k] * ab_y + g_0_x_xyy_yyyzz[k];

                g_0_x_xyyy_yzzz[k] = -g_0_x_xyy_yzzz[k] * ab_y + g_0_x_xyy_yyzzz[k];

                g_0_x_xyyy_zzzz[k] = -g_0_x_xyy_zzzz[k] * ab_y + g_0_x_xyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyz_xxxx = cbuffer.data(gg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxy = cbuffer.data(gg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxz = cbuffer.data(gg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyy = cbuffer.data(gg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyz = cbuffer.data(gg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xyyz_xxzz = cbuffer.data(gg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyy = cbuffer.data(gg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyz = cbuffer.data(gg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xyyz_xyzz = cbuffer.data(gg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xyyz_xzzz = cbuffer.data(gg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyy = cbuffer.data(gg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyz = cbuffer.data(gg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xyyz_yyzz = cbuffer.data(gg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xyyz_yzzz = cbuffer.data(gg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xyyz_zzzz = cbuffer.data(gg_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyz_xxxx, g_0_x_xyyz_xxxy, g_0_x_xyyz_xxxz, g_0_x_xyyz_xxyy, g_0_x_xyyz_xxyz, g_0_x_xyyz_xxzz, g_0_x_xyyz_xyyy, g_0_x_xyyz_xyyz, g_0_x_xyyz_xyzz, g_0_x_xyyz_xzzz, g_0_x_xyyz_yyyy, g_0_x_xyyz_yyyz, g_0_x_xyyz_yyzz, g_0_x_xyyz_yzzz, g_0_x_xyyz_zzzz, g_0_x_xyz_xxxx, g_0_x_xyz_xxxxy, g_0_x_xyz_xxxy, g_0_x_xyz_xxxyy, g_0_x_xyz_xxxyz, g_0_x_xyz_xxxz, g_0_x_xyz_xxyy, g_0_x_xyz_xxyyy, g_0_x_xyz_xxyyz, g_0_x_xyz_xxyz, g_0_x_xyz_xxyzz, g_0_x_xyz_xxzz, g_0_x_xyz_xyyy, g_0_x_xyz_xyyyy, g_0_x_xyz_xyyyz, g_0_x_xyz_xyyz, g_0_x_xyz_xyyzz, g_0_x_xyz_xyzz, g_0_x_xyz_xyzzz, g_0_x_xyz_xzzz, g_0_x_xyz_yyyy, g_0_x_xyz_yyyyy, g_0_x_xyz_yyyyz, g_0_x_xyz_yyyz, g_0_x_xyz_yyyzz, g_0_x_xyz_yyzz, g_0_x_xyz_yyzzz, g_0_x_xyz_yzzz, g_0_x_xyz_yzzzz, g_0_x_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyz_xxxx[k] = -g_0_x_xyz_xxxx[k] * ab_y + g_0_x_xyz_xxxxy[k];

                g_0_x_xyyz_xxxy[k] = -g_0_x_xyz_xxxy[k] * ab_y + g_0_x_xyz_xxxyy[k];

                g_0_x_xyyz_xxxz[k] = -g_0_x_xyz_xxxz[k] * ab_y + g_0_x_xyz_xxxyz[k];

                g_0_x_xyyz_xxyy[k] = -g_0_x_xyz_xxyy[k] * ab_y + g_0_x_xyz_xxyyy[k];

                g_0_x_xyyz_xxyz[k] = -g_0_x_xyz_xxyz[k] * ab_y + g_0_x_xyz_xxyyz[k];

                g_0_x_xyyz_xxzz[k] = -g_0_x_xyz_xxzz[k] * ab_y + g_0_x_xyz_xxyzz[k];

                g_0_x_xyyz_xyyy[k] = -g_0_x_xyz_xyyy[k] * ab_y + g_0_x_xyz_xyyyy[k];

                g_0_x_xyyz_xyyz[k] = -g_0_x_xyz_xyyz[k] * ab_y + g_0_x_xyz_xyyyz[k];

                g_0_x_xyyz_xyzz[k] = -g_0_x_xyz_xyzz[k] * ab_y + g_0_x_xyz_xyyzz[k];

                g_0_x_xyyz_xzzz[k] = -g_0_x_xyz_xzzz[k] * ab_y + g_0_x_xyz_xyzzz[k];

                g_0_x_xyyz_yyyy[k] = -g_0_x_xyz_yyyy[k] * ab_y + g_0_x_xyz_yyyyy[k];

                g_0_x_xyyz_yyyz[k] = -g_0_x_xyz_yyyz[k] * ab_y + g_0_x_xyz_yyyyz[k];

                g_0_x_xyyz_yyzz[k] = -g_0_x_xyz_yyzz[k] * ab_y + g_0_x_xyz_yyyzz[k];

                g_0_x_xyyz_yzzz[k] = -g_0_x_xyz_yzzz[k] * ab_y + g_0_x_xyz_yyzzz[k];

                g_0_x_xyyz_zzzz[k] = -g_0_x_xyz_zzzz[k] * ab_y + g_0_x_xyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzz_xxxx = cbuffer.data(gg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxy = cbuffer.data(gg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxz = cbuffer.data(gg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyy = cbuffer.data(gg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyz = cbuffer.data(gg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xyzz_xxzz = cbuffer.data(gg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyy = cbuffer.data(gg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyz = cbuffer.data(gg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xyzz_xyzz = cbuffer.data(gg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xyzz_xzzz = cbuffer.data(gg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyy = cbuffer.data(gg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyz = cbuffer.data(gg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xyzz_yyzz = cbuffer.data(gg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xyzz_yzzz = cbuffer.data(gg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xyzz_zzzz = cbuffer.data(gg_geom_01_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzz_xxxx, g_0_x_xyzz_xxxy, g_0_x_xyzz_xxxz, g_0_x_xyzz_xxyy, g_0_x_xyzz_xxyz, g_0_x_xyzz_xxzz, g_0_x_xyzz_xyyy, g_0_x_xyzz_xyyz, g_0_x_xyzz_xyzz, g_0_x_xyzz_xzzz, g_0_x_xyzz_yyyy, g_0_x_xyzz_yyyz, g_0_x_xyzz_yyzz, g_0_x_xyzz_yzzz, g_0_x_xyzz_zzzz, g_0_x_xzz_xxxx, g_0_x_xzz_xxxxy, g_0_x_xzz_xxxy, g_0_x_xzz_xxxyy, g_0_x_xzz_xxxyz, g_0_x_xzz_xxxz, g_0_x_xzz_xxyy, g_0_x_xzz_xxyyy, g_0_x_xzz_xxyyz, g_0_x_xzz_xxyz, g_0_x_xzz_xxyzz, g_0_x_xzz_xxzz, g_0_x_xzz_xyyy, g_0_x_xzz_xyyyy, g_0_x_xzz_xyyyz, g_0_x_xzz_xyyz, g_0_x_xzz_xyyzz, g_0_x_xzz_xyzz, g_0_x_xzz_xyzzz, g_0_x_xzz_xzzz, g_0_x_xzz_yyyy, g_0_x_xzz_yyyyy, g_0_x_xzz_yyyyz, g_0_x_xzz_yyyz, g_0_x_xzz_yyyzz, g_0_x_xzz_yyzz, g_0_x_xzz_yyzzz, g_0_x_xzz_yzzz, g_0_x_xzz_yzzzz, g_0_x_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzz_xxxx[k] = -g_0_x_xzz_xxxx[k] * ab_y + g_0_x_xzz_xxxxy[k];

                g_0_x_xyzz_xxxy[k] = -g_0_x_xzz_xxxy[k] * ab_y + g_0_x_xzz_xxxyy[k];

                g_0_x_xyzz_xxxz[k] = -g_0_x_xzz_xxxz[k] * ab_y + g_0_x_xzz_xxxyz[k];

                g_0_x_xyzz_xxyy[k] = -g_0_x_xzz_xxyy[k] * ab_y + g_0_x_xzz_xxyyy[k];

                g_0_x_xyzz_xxyz[k] = -g_0_x_xzz_xxyz[k] * ab_y + g_0_x_xzz_xxyyz[k];

                g_0_x_xyzz_xxzz[k] = -g_0_x_xzz_xxzz[k] * ab_y + g_0_x_xzz_xxyzz[k];

                g_0_x_xyzz_xyyy[k] = -g_0_x_xzz_xyyy[k] * ab_y + g_0_x_xzz_xyyyy[k];

                g_0_x_xyzz_xyyz[k] = -g_0_x_xzz_xyyz[k] * ab_y + g_0_x_xzz_xyyyz[k];

                g_0_x_xyzz_xyzz[k] = -g_0_x_xzz_xyzz[k] * ab_y + g_0_x_xzz_xyyzz[k];

                g_0_x_xyzz_xzzz[k] = -g_0_x_xzz_xzzz[k] * ab_y + g_0_x_xzz_xyzzz[k];

                g_0_x_xyzz_yyyy[k] = -g_0_x_xzz_yyyy[k] * ab_y + g_0_x_xzz_yyyyy[k];

                g_0_x_xyzz_yyyz[k] = -g_0_x_xzz_yyyz[k] * ab_y + g_0_x_xzz_yyyyz[k];

                g_0_x_xyzz_yyzz[k] = -g_0_x_xzz_yyzz[k] * ab_y + g_0_x_xzz_yyyzz[k];

                g_0_x_xyzz_yzzz[k] = -g_0_x_xzz_yzzz[k] * ab_y + g_0_x_xzz_yyzzz[k];

                g_0_x_xyzz_zzzz[k] = -g_0_x_xzz_zzzz[k] * ab_y + g_0_x_xzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzz_xxxx = cbuffer.data(gg_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxy = cbuffer.data(gg_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxz = cbuffer.data(gg_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyy = cbuffer.data(gg_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyz = cbuffer.data(gg_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xzzz_xxzz = cbuffer.data(gg_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyy = cbuffer.data(gg_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyz = cbuffer.data(gg_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xzzz_xyzz = cbuffer.data(gg_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xzzz_xzzz = cbuffer.data(gg_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyy = cbuffer.data(gg_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyz = cbuffer.data(gg_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xzzz_yyzz = cbuffer.data(gg_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xzzz_yzzz = cbuffer.data(gg_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xzzz_zzzz = cbuffer.data(gg_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzz_xxxx, g_0_x_xzz_xxxxz, g_0_x_xzz_xxxy, g_0_x_xzz_xxxyz, g_0_x_xzz_xxxz, g_0_x_xzz_xxxzz, g_0_x_xzz_xxyy, g_0_x_xzz_xxyyz, g_0_x_xzz_xxyz, g_0_x_xzz_xxyzz, g_0_x_xzz_xxzz, g_0_x_xzz_xxzzz, g_0_x_xzz_xyyy, g_0_x_xzz_xyyyz, g_0_x_xzz_xyyz, g_0_x_xzz_xyyzz, g_0_x_xzz_xyzz, g_0_x_xzz_xyzzz, g_0_x_xzz_xzzz, g_0_x_xzz_xzzzz, g_0_x_xzz_yyyy, g_0_x_xzz_yyyyz, g_0_x_xzz_yyyz, g_0_x_xzz_yyyzz, g_0_x_xzz_yyzz, g_0_x_xzz_yyzzz, g_0_x_xzz_yzzz, g_0_x_xzz_yzzzz, g_0_x_xzz_zzzz, g_0_x_xzz_zzzzz, g_0_x_xzzz_xxxx, g_0_x_xzzz_xxxy, g_0_x_xzzz_xxxz, g_0_x_xzzz_xxyy, g_0_x_xzzz_xxyz, g_0_x_xzzz_xxzz, g_0_x_xzzz_xyyy, g_0_x_xzzz_xyyz, g_0_x_xzzz_xyzz, g_0_x_xzzz_xzzz, g_0_x_xzzz_yyyy, g_0_x_xzzz_yyyz, g_0_x_xzzz_yyzz, g_0_x_xzzz_yzzz, g_0_x_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzz_xxxx[k] = -g_0_x_xzz_xxxx[k] * ab_z + g_0_x_xzz_xxxxz[k];

                g_0_x_xzzz_xxxy[k] = -g_0_x_xzz_xxxy[k] * ab_z + g_0_x_xzz_xxxyz[k];

                g_0_x_xzzz_xxxz[k] = -g_0_x_xzz_xxxz[k] * ab_z + g_0_x_xzz_xxxzz[k];

                g_0_x_xzzz_xxyy[k] = -g_0_x_xzz_xxyy[k] * ab_z + g_0_x_xzz_xxyyz[k];

                g_0_x_xzzz_xxyz[k] = -g_0_x_xzz_xxyz[k] * ab_z + g_0_x_xzz_xxyzz[k];

                g_0_x_xzzz_xxzz[k] = -g_0_x_xzz_xxzz[k] * ab_z + g_0_x_xzz_xxzzz[k];

                g_0_x_xzzz_xyyy[k] = -g_0_x_xzz_xyyy[k] * ab_z + g_0_x_xzz_xyyyz[k];

                g_0_x_xzzz_xyyz[k] = -g_0_x_xzz_xyyz[k] * ab_z + g_0_x_xzz_xyyzz[k];

                g_0_x_xzzz_xyzz[k] = -g_0_x_xzz_xyzz[k] * ab_z + g_0_x_xzz_xyzzz[k];

                g_0_x_xzzz_xzzz[k] = -g_0_x_xzz_xzzz[k] * ab_z + g_0_x_xzz_xzzzz[k];

                g_0_x_xzzz_yyyy[k] = -g_0_x_xzz_yyyy[k] * ab_z + g_0_x_xzz_yyyyz[k];

                g_0_x_xzzz_yyyz[k] = -g_0_x_xzz_yyyz[k] * ab_z + g_0_x_xzz_yyyzz[k];

                g_0_x_xzzz_yyzz[k] = -g_0_x_xzz_yyzz[k] * ab_z + g_0_x_xzz_yyzzz[k];

                g_0_x_xzzz_yzzz[k] = -g_0_x_xzz_yzzz[k] * ab_z + g_0_x_xzz_yzzzz[k];

                g_0_x_xzzz_zzzz[k] = -g_0_x_xzz_zzzz[k] * ab_z + g_0_x_xzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyy_xxxx = cbuffer.data(gg_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxy = cbuffer.data(gg_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxz = cbuffer.data(gg_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyy = cbuffer.data(gg_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyz = cbuffer.data(gg_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yyyy_xxzz = cbuffer.data(gg_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyy = cbuffer.data(gg_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyz = cbuffer.data(gg_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yyyy_xyzz = cbuffer.data(gg_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yyyy_xzzz = cbuffer.data(gg_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyy = cbuffer.data(gg_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyz = cbuffer.data(gg_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_yyyy_yyzz = cbuffer.data(gg_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_yyyy_yzzz = cbuffer.data(gg_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_yyyy_zzzz = cbuffer.data(gg_geom_01_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyy_xxxx, g_0_x_yyy_xxxxy, g_0_x_yyy_xxxy, g_0_x_yyy_xxxyy, g_0_x_yyy_xxxyz, g_0_x_yyy_xxxz, g_0_x_yyy_xxyy, g_0_x_yyy_xxyyy, g_0_x_yyy_xxyyz, g_0_x_yyy_xxyz, g_0_x_yyy_xxyzz, g_0_x_yyy_xxzz, g_0_x_yyy_xyyy, g_0_x_yyy_xyyyy, g_0_x_yyy_xyyyz, g_0_x_yyy_xyyz, g_0_x_yyy_xyyzz, g_0_x_yyy_xyzz, g_0_x_yyy_xyzzz, g_0_x_yyy_xzzz, g_0_x_yyy_yyyy, g_0_x_yyy_yyyyy, g_0_x_yyy_yyyyz, g_0_x_yyy_yyyz, g_0_x_yyy_yyyzz, g_0_x_yyy_yyzz, g_0_x_yyy_yyzzz, g_0_x_yyy_yzzz, g_0_x_yyy_yzzzz, g_0_x_yyy_zzzz, g_0_x_yyyy_xxxx, g_0_x_yyyy_xxxy, g_0_x_yyyy_xxxz, g_0_x_yyyy_xxyy, g_0_x_yyyy_xxyz, g_0_x_yyyy_xxzz, g_0_x_yyyy_xyyy, g_0_x_yyyy_xyyz, g_0_x_yyyy_xyzz, g_0_x_yyyy_xzzz, g_0_x_yyyy_yyyy, g_0_x_yyyy_yyyz, g_0_x_yyyy_yyzz, g_0_x_yyyy_yzzz, g_0_x_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyy_xxxx[k] = -g_0_x_yyy_xxxx[k] * ab_y + g_0_x_yyy_xxxxy[k];

                g_0_x_yyyy_xxxy[k] = -g_0_x_yyy_xxxy[k] * ab_y + g_0_x_yyy_xxxyy[k];

                g_0_x_yyyy_xxxz[k] = -g_0_x_yyy_xxxz[k] * ab_y + g_0_x_yyy_xxxyz[k];

                g_0_x_yyyy_xxyy[k] = -g_0_x_yyy_xxyy[k] * ab_y + g_0_x_yyy_xxyyy[k];

                g_0_x_yyyy_xxyz[k] = -g_0_x_yyy_xxyz[k] * ab_y + g_0_x_yyy_xxyyz[k];

                g_0_x_yyyy_xxzz[k] = -g_0_x_yyy_xxzz[k] * ab_y + g_0_x_yyy_xxyzz[k];

                g_0_x_yyyy_xyyy[k] = -g_0_x_yyy_xyyy[k] * ab_y + g_0_x_yyy_xyyyy[k];

                g_0_x_yyyy_xyyz[k] = -g_0_x_yyy_xyyz[k] * ab_y + g_0_x_yyy_xyyyz[k];

                g_0_x_yyyy_xyzz[k] = -g_0_x_yyy_xyzz[k] * ab_y + g_0_x_yyy_xyyzz[k];

                g_0_x_yyyy_xzzz[k] = -g_0_x_yyy_xzzz[k] * ab_y + g_0_x_yyy_xyzzz[k];

                g_0_x_yyyy_yyyy[k] = -g_0_x_yyy_yyyy[k] * ab_y + g_0_x_yyy_yyyyy[k];

                g_0_x_yyyy_yyyz[k] = -g_0_x_yyy_yyyz[k] * ab_y + g_0_x_yyy_yyyyz[k];

                g_0_x_yyyy_yyzz[k] = -g_0_x_yyy_yyzz[k] * ab_y + g_0_x_yyy_yyyzz[k];

                g_0_x_yyyy_yzzz[k] = -g_0_x_yyy_yzzz[k] * ab_y + g_0_x_yyy_yyzzz[k];

                g_0_x_yyyy_zzzz[k] = -g_0_x_yyy_zzzz[k] * ab_y + g_0_x_yyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyz_xxxx = cbuffer.data(gg_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxy = cbuffer.data(gg_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxz = cbuffer.data(gg_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyy = cbuffer.data(gg_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyz = cbuffer.data(gg_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_yyyz_xxzz = cbuffer.data(gg_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyy = cbuffer.data(gg_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyz = cbuffer.data(gg_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_yyyz_xyzz = cbuffer.data(gg_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_yyyz_xzzz = cbuffer.data(gg_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyy = cbuffer.data(gg_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyz = cbuffer.data(gg_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_yyyz_yyzz = cbuffer.data(gg_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_yyyz_yzzz = cbuffer.data(gg_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_yyyz_zzzz = cbuffer.data(gg_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyz_xxxx, g_0_x_yyyz_xxxy, g_0_x_yyyz_xxxz, g_0_x_yyyz_xxyy, g_0_x_yyyz_xxyz, g_0_x_yyyz_xxzz, g_0_x_yyyz_xyyy, g_0_x_yyyz_xyyz, g_0_x_yyyz_xyzz, g_0_x_yyyz_xzzz, g_0_x_yyyz_yyyy, g_0_x_yyyz_yyyz, g_0_x_yyyz_yyzz, g_0_x_yyyz_yzzz, g_0_x_yyyz_zzzz, g_0_x_yyz_xxxx, g_0_x_yyz_xxxxy, g_0_x_yyz_xxxy, g_0_x_yyz_xxxyy, g_0_x_yyz_xxxyz, g_0_x_yyz_xxxz, g_0_x_yyz_xxyy, g_0_x_yyz_xxyyy, g_0_x_yyz_xxyyz, g_0_x_yyz_xxyz, g_0_x_yyz_xxyzz, g_0_x_yyz_xxzz, g_0_x_yyz_xyyy, g_0_x_yyz_xyyyy, g_0_x_yyz_xyyyz, g_0_x_yyz_xyyz, g_0_x_yyz_xyyzz, g_0_x_yyz_xyzz, g_0_x_yyz_xyzzz, g_0_x_yyz_xzzz, g_0_x_yyz_yyyy, g_0_x_yyz_yyyyy, g_0_x_yyz_yyyyz, g_0_x_yyz_yyyz, g_0_x_yyz_yyyzz, g_0_x_yyz_yyzz, g_0_x_yyz_yyzzz, g_0_x_yyz_yzzz, g_0_x_yyz_yzzzz, g_0_x_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyz_xxxx[k] = -g_0_x_yyz_xxxx[k] * ab_y + g_0_x_yyz_xxxxy[k];

                g_0_x_yyyz_xxxy[k] = -g_0_x_yyz_xxxy[k] * ab_y + g_0_x_yyz_xxxyy[k];

                g_0_x_yyyz_xxxz[k] = -g_0_x_yyz_xxxz[k] * ab_y + g_0_x_yyz_xxxyz[k];

                g_0_x_yyyz_xxyy[k] = -g_0_x_yyz_xxyy[k] * ab_y + g_0_x_yyz_xxyyy[k];

                g_0_x_yyyz_xxyz[k] = -g_0_x_yyz_xxyz[k] * ab_y + g_0_x_yyz_xxyyz[k];

                g_0_x_yyyz_xxzz[k] = -g_0_x_yyz_xxzz[k] * ab_y + g_0_x_yyz_xxyzz[k];

                g_0_x_yyyz_xyyy[k] = -g_0_x_yyz_xyyy[k] * ab_y + g_0_x_yyz_xyyyy[k];

                g_0_x_yyyz_xyyz[k] = -g_0_x_yyz_xyyz[k] * ab_y + g_0_x_yyz_xyyyz[k];

                g_0_x_yyyz_xyzz[k] = -g_0_x_yyz_xyzz[k] * ab_y + g_0_x_yyz_xyyzz[k];

                g_0_x_yyyz_xzzz[k] = -g_0_x_yyz_xzzz[k] * ab_y + g_0_x_yyz_xyzzz[k];

                g_0_x_yyyz_yyyy[k] = -g_0_x_yyz_yyyy[k] * ab_y + g_0_x_yyz_yyyyy[k];

                g_0_x_yyyz_yyyz[k] = -g_0_x_yyz_yyyz[k] * ab_y + g_0_x_yyz_yyyyz[k];

                g_0_x_yyyz_yyzz[k] = -g_0_x_yyz_yyzz[k] * ab_y + g_0_x_yyz_yyyzz[k];

                g_0_x_yyyz_yzzz[k] = -g_0_x_yyz_yzzz[k] * ab_y + g_0_x_yyz_yyzzz[k];

                g_0_x_yyyz_zzzz[k] = -g_0_x_yyz_zzzz[k] * ab_y + g_0_x_yyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzz_xxxx = cbuffer.data(gg_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxy = cbuffer.data(gg_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxz = cbuffer.data(gg_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyy = cbuffer.data(gg_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyz = cbuffer.data(gg_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_yyzz_xxzz = cbuffer.data(gg_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyy = cbuffer.data(gg_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyz = cbuffer.data(gg_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_yyzz_xyzz = cbuffer.data(gg_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_yyzz_xzzz = cbuffer.data(gg_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyy = cbuffer.data(gg_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyz = cbuffer.data(gg_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_yyzz_yyzz = cbuffer.data(gg_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_yyzz_yzzz = cbuffer.data(gg_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_yyzz_zzzz = cbuffer.data(gg_geom_01_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzz_xxxx, g_0_x_yyzz_xxxy, g_0_x_yyzz_xxxz, g_0_x_yyzz_xxyy, g_0_x_yyzz_xxyz, g_0_x_yyzz_xxzz, g_0_x_yyzz_xyyy, g_0_x_yyzz_xyyz, g_0_x_yyzz_xyzz, g_0_x_yyzz_xzzz, g_0_x_yyzz_yyyy, g_0_x_yyzz_yyyz, g_0_x_yyzz_yyzz, g_0_x_yyzz_yzzz, g_0_x_yyzz_zzzz, g_0_x_yzz_xxxx, g_0_x_yzz_xxxxy, g_0_x_yzz_xxxy, g_0_x_yzz_xxxyy, g_0_x_yzz_xxxyz, g_0_x_yzz_xxxz, g_0_x_yzz_xxyy, g_0_x_yzz_xxyyy, g_0_x_yzz_xxyyz, g_0_x_yzz_xxyz, g_0_x_yzz_xxyzz, g_0_x_yzz_xxzz, g_0_x_yzz_xyyy, g_0_x_yzz_xyyyy, g_0_x_yzz_xyyyz, g_0_x_yzz_xyyz, g_0_x_yzz_xyyzz, g_0_x_yzz_xyzz, g_0_x_yzz_xyzzz, g_0_x_yzz_xzzz, g_0_x_yzz_yyyy, g_0_x_yzz_yyyyy, g_0_x_yzz_yyyyz, g_0_x_yzz_yyyz, g_0_x_yzz_yyyzz, g_0_x_yzz_yyzz, g_0_x_yzz_yyzzz, g_0_x_yzz_yzzz, g_0_x_yzz_yzzzz, g_0_x_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzz_xxxx[k] = -g_0_x_yzz_xxxx[k] * ab_y + g_0_x_yzz_xxxxy[k];

                g_0_x_yyzz_xxxy[k] = -g_0_x_yzz_xxxy[k] * ab_y + g_0_x_yzz_xxxyy[k];

                g_0_x_yyzz_xxxz[k] = -g_0_x_yzz_xxxz[k] * ab_y + g_0_x_yzz_xxxyz[k];

                g_0_x_yyzz_xxyy[k] = -g_0_x_yzz_xxyy[k] * ab_y + g_0_x_yzz_xxyyy[k];

                g_0_x_yyzz_xxyz[k] = -g_0_x_yzz_xxyz[k] * ab_y + g_0_x_yzz_xxyyz[k];

                g_0_x_yyzz_xxzz[k] = -g_0_x_yzz_xxzz[k] * ab_y + g_0_x_yzz_xxyzz[k];

                g_0_x_yyzz_xyyy[k] = -g_0_x_yzz_xyyy[k] * ab_y + g_0_x_yzz_xyyyy[k];

                g_0_x_yyzz_xyyz[k] = -g_0_x_yzz_xyyz[k] * ab_y + g_0_x_yzz_xyyyz[k];

                g_0_x_yyzz_xyzz[k] = -g_0_x_yzz_xyzz[k] * ab_y + g_0_x_yzz_xyyzz[k];

                g_0_x_yyzz_xzzz[k] = -g_0_x_yzz_xzzz[k] * ab_y + g_0_x_yzz_xyzzz[k];

                g_0_x_yyzz_yyyy[k] = -g_0_x_yzz_yyyy[k] * ab_y + g_0_x_yzz_yyyyy[k];

                g_0_x_yyzz_yyyz[k] = -g_0_x_yzz_yyyz[k] * ab_y + g_0_x_yzz_yyyyz[k];

                g_0_x_yyzz_yyzz[k] = -g_0_x_yzz_yyzz[k] * ab_y + g_0_x_yzz_yyyzz[k];

                g_0_x_yyzz_yzzz[k] = -g_0_x_yzz_yzzz[k] * ab_y + g_0_x_yzz_yyzzz[k];

                g_0_x_yyzz_zzzz[k] = -g_0_x_yzz_zzzz[k] * ab_y + g_0_x_yzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzz_xxxx = cbuffer.data(gg_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxy = cbuffer.data(gg_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxz = cbuffer.data(gg_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyy = cbuffer.data(gg_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyz = cbuffer.data(gg_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_yzzz_xxzz = cbuffer.data(gg_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyy = cbuffer.data(gg_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyz = cbuffer.data(gg_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_yzzz_xyzz = cbuffer.data(gg_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_yzzz_xzzz = cbuffer.data(gg_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyy = cbuffer.data(gg_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyz = cbuffer.data(gg_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_yzzz_yyzz = cbuffer.data(gg_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_yzzz_yzzz = cbuffer.data(gg_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_yzzz_zzzz = cbuffer.data(gg_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzz_xxxx, g_0_x_yzzz_xxxy, g_0_x_yzzz_xxxz, g_0_x_yzzz_xxyy, g_0_x_yzzz_xxyz, g_0_x_yzzz_xxzz, g_0_x_yzzz_xyyy, g_0_x_yzzz_xyyz, g_0_x_yzzz_xyzz, g_0_x_yzzz_xzzz, g_0_x_yzzz_yyyy, g_0_x_yzzz_yyyz, g_0_x_yzzz_yyzz, g_0_x_yzzz_yzzz, g_0_x_yzzz_zzzz, g_0_x_zzz_xxxx, g_0_x_zzz_xxxxy, g_0_x_zzz_xxxy, g_0_x_zzz_xxxyy, g_0_x_zzz_xxxyz, g_0_x_zzz_xxxz, g_0_x_zzz_xxyy, g_0_x_zzz_xxyyy, g_0_x_zzz_xxyyz, g_0_x_zzz_xxyz, g_0_x_zzz_xxyzz, g_0_x_zzz_xxzz, g_0_x_zzz_xyyy, g_0_x_zzz_xyyyy, g_0_x_zzz_xyyyz, g_0_x_zzz_xyyz, g_0_x_zzz_xyyzz, g_0_x_zzz_xyzz, g_0_x_zzz_xyzzz, g_0_x_zzz_xzzz, g_0_x_zzz_yyyy, g_0_x_zzz_yyyyy, g_0_x_zzz_yyyyz, g_0_x_zzz_yyyz, g_0_x_zzz_yyyzz, g_0_x_zzz_yyzz, g_0_x_zzz_yyzzz, g_0_x_zzz_yzzz, g_0_x_zzz_yzzzz, g_0_x_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzz_xxxx[k] = -g_0_x_zzz_xxxx[k] * ab_y + g_0_x_zzz_xxxxy[k];

                g_0_x_yzzz_xxxy[k] = -g_0_x_zzz_xxxy[k] * ab_y + g_0_x_zzz_xxxyy[k];

                g_0_x_yzzz_xxxz[k] = -g_0_x_zzz_xxxz[k] * ab_y + g_0_x_zzz_xxxyz[k];

                g_0_x_yzzz_xxyy[k] = -g_0_x_zzz_xxyy[k] * ab_y + g_0_x_zzz_xxyyy[k];

                g_0_x_yzzz_xxyz[k] = -g_0_x_zzz_xxyz[k] * ab_y + g_0_x_zzz_xxyyz[k];

                g_0_x_yzzz_xxzz[k] = -g_0_x_zzz_xxzz[k] * ab_y + g_0_x_zzz_xxyzz[k];

                g_0_x_yzzz_xyyy[k] = -g_0_x_zzz_xyyy[k] * ab_y + g_0_x_zzz_xyyyy[k];

                g_0_x_yzzz_xyyz[k] = -g_0_x_zzz_xyyz[k] * ab_y + g_0_x_zzz_xyyyz[k];

                g_0_x_yzzz_xyzz[k] = -g_0_x_zzz_xyzz[k] * ab_y + g_0_x_zzz_xyyzz[k];

                g_0_x_yzzz_xzzz[k] = -g_0_x_zzz_xzzz[k] * ab_y + g_0_x_zzz_xyzzz[k];

                g_0_x_yzzz_yyyy[k] = -g_0_x_zzz_yyyy[k] * ab_y + g_0_x_zzz_yyyyy[k];

                g_0_x_yzzz_yyyz[k] = -g_0_x_zzz_yyyz[k] * ab_y + g_0_x_zzz_yyyyz[k];

                g_0_x_yzzz_yyzz[k] = -g_0_x_zzz_yyzz[k] * ab_y + g_0_x_zzz_yyyzz[k];

                g_0_x_yzzz_yzzz[k] = -g_0_x_zzz_yzzz[k] * ab_y + g_0_x_zzz_yyzzz[k];

                g_0_x_yzzz_zzzz[k] = -g_0_x_zzz_zzzz[k] * ab_y + g_0_x_zzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzz_xxxx = cbuffer.data(gg_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxy = cbuffer.data(gg_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxz = cbuffer.data(gg_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyy = cbuffer.data(gg_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyz = cbuffer.data(gg_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_zzzz_xxzz = cbuffer.data(gg_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyy = cbuffer.data(gg_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyz = cbuffer.data(gg_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_zzzz_xyzz = cbuffer.data(gg_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_zzzz_xzzz = cbuffer.data(gg_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyy = cbuffer.data(gg_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyz = cbuffer.data(gg_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_zzzz_yyzz = cbuffer.data(gg_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_zzzz_yzzz = cbuffer.data(gg_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_zzzz_zzzz = cbuffer.data(gg_geom_01_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzz_xxxx, g_0_x_zzz_xxxxz, g_0_x_zzz_xxxy, g_0_x_zzz_xxxyz, g_0_x_zzz_xxxz, g_0_x_zzz_xxxzz, g_0_x_zzz_xxyy, g_0_x_zzz_xxyyz, g_0_x_zzz_xxyz, g_0_x_zzz_xxyzz, g_0_x_zzz_xxzz, g_0_x_zzz_xxzzz, g_0_x_zzz_xyyy, g_0_x_zzz_xyyyz, g_0_x_zzz_xyyz, g_0_x_zzz_xyyzz, g_0_x_zzz_xyzz, g_0_x_zzz_xyzzz, g_0_x_zzz_xzzz, g_0_x_zzz_xzzzz, g_0_x_zzz_yyyy, g_0_x_zzz_yyyyz, g_0_x_zzz_yyyz, g_0_x_zzz_yyyzz, g_0_x_zzz_yyzz, g_0_x_zzz_yyzzz, g_0_x_zzz_yzzz, g_0_x_zzz_yzzzz, g_0_x_zzz_zzzz, g_0_x_zzz_zzzzz, g_0_x_zzzz_xxxx, g_0_x_zzzz_xxxy, g_0_x_zzzz_xxxz, g_0_x_zzzz_xxyy, g_0_x_zzzz_xxyz, g_0_x_zzzz_xxzz, g_0_x_zzzz_xyyy, g_0_x_zzzz_xyyz, g_0_x_zzzz_xyzz, g_0_x_zzzz_xzzz, g_0_x_zzzz_yyyy, g_0_x_zzzz_yyyz, g_0_x_zzzz_yyzz, g_0_x_zzzz_yzzz, g_0_x_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzz_xxxx[k] = -g_0_x_zzz_xxxx[k] * ab_z + g_0_x_zzz_xxxxz[k];

                g_0_x_zzzz_xxxy[k] = -g_0_x_zzz_xxxy[k] * ab_z + g_0_x_zzz_xxxyz[k];

                g_0_x_zzzz_xxxz[k] = -g_0_x_zzz_xxxz[k] * ab_z + g_0_x_zzz_xxxzz[k];

                g_0_x_zzzz_xxyy[k] = -g_0_x_zzz_xxyy[k] * ab_z + g_0_x_zzz_xxyyz[k];

                g_0_x_zzzz_xxyz[k] = -g_0_x_zzz_xxyz[k] * ab_z + g_0_x_zzz_xxyzz[k];

                g_0_x_zzzz_xxzz[k] = -g_0_x_zzz_xxzz[k] * ab_z + g_0_x_zzz_xxzzz[k];

                g_0_x_zzzz_xyyy[k] = -g_0_x_zzz_xyyy[k] * ab_z + g_0_x_zzz_xyyyz[k];

                g_0_x_zzzz_xyyz[k] = -g_0_x_zzz_xyyz[k] * ab_z + g_0_x_zzz_xyyzz[k];

                g_0_x_zzzz_xyzz[k] = -g_0_x_zzz_xyzz[k] * ab_z + g_0_x_zzz_xyzzz[k];

                g_0_x_zzzz_xzzz[k] = -g_0_x_zzz_xzzz[k] * ab_z + g_0_x_zzz_xzzzz[k];

                g_0_x_zzzz_yyyy[k] = -g_0_x_zzz_yyyy[k] * ab_z + g_0_x_zzz_yyyyz[k];

                g_0_x_zzzz_yyyz[k] = -g_0_x_zzz_yyyz[k] * ab_z + g_0_x_zzz_yyyzz[k];

                g_0_x_zzzz_yyzz[k] = -g_0_x_zzz_yyzz[k] * ab_z + g_0_x_zzz_yyzzz[k];

                g_0_x_zzzz_yzzz[k] = -g_0_x_zzz_yzzz[k] * ab_z + g_0_x_zzz_yzzzz[k];

                g_0_x_zzzz_zzzz[k] = -g_0_x_zzz_zzzz[k] * ab_z + g_0_x_zzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxx_xxxx = cbuffer.data(gg_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxy = cbuffer.data(gg_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxz = cbuffer.data(gg_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyy = cbuffer.data(gg_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyz = cbuffer.data(gg_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xxxx_xxzz = cbuffer.data(gg_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyy = cbuffer.data(gg_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyz = cbuffer.data(gg_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xxxx_xyzz = cbuffer.data(gg_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xxxx_xzzz = cbuffer.data(gg_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyy = cbuffer.data(gg_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyz = cbuffer.data(gg_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xxxx_yyzz = cbuffer.data(gg_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xxxx_yzzz = cbuffer.data(gg_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xxxx_zzzz = cbuffer.data(gg_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxx_xxxx, g_0_y_xxx_xxxxx, g_0_y_xxx_xxxxy, g_0_y_xxx_xxxxz, g_0_y_xxx_xxxy, g_0_y_xxx_xxxyy, g_0_y_xxx_xxxyz, g_0_y_xxx_xxxz, g_0_y_xxx_xxxzz, g_0_y_xxx_xxyy, g_0_y_xxx_xxyyy, g_0_y_xxx_xxyyz, g_0_y_xxx_xxyz, g_0_y_xxx_xxyzz, g_0_y_xxx_xxzz, g_0_y_xxx_xxzzz, g_0_y_xxx_xyyy, g_0_y_xxx_xyyyy, g_0_y_xxx_xyyyz, g_0_y_xxx_xyyz, g_0_y_xxx_xyyzz, g_0_y_xxx_xyzz, g_0_y_xxx_xyzzz, g_0_y_xxx_xzzz, g_0_y_xxx_xzzzz, g_0_y_xxx_yyyy, g_0_y_xxx_yyyz, g_0_y_xxx_yyzz, g_0_y_xxx_yzzz, g_0_y_xxx_zzzz, g_0_y_xxxx_xxxx, g_0_y_xxxx_xxxy, g_0_y_xxxx_xxxz, g_0_y_xxxx_xxyy, g_0_y_xxxx_xxyz, g_0_y_xxxx_xxzz, g_0_y_xxxx_xyyy, g_0_y_xxxx_xyyz, g_0_y_xxxx_xyzz, g_0_y_xxxx_xzzz, g_0_y_xxxx_yyyy, g_0_y_xxxx_yyyz, g_0_y_xxxx_yyzz, g_0_y_xxxx_yzzz, g_0_y_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxx_xxxx[k] = -g_0_y_xxx_xxxx[k] * ab_x + g_0_y_xxx_xxxxx[k];

                g_0_y_xxxx_xxxy[k] = -g_0_y_xxx_xxxy[k] * ab_x + g_0_y_xxx_xxxxy[k];

                g_0_y_xxxx_xxxz[k] = -g_0_y_xxx_xxxz[k] * ab_x + g_0_y_xxx_xxxxz[k];

                g_0_y_xxxx_xxyy[k] = -g_0_y_xxx_xxyy[k] * ab_x + g_0_y_xxx_xxxyy[k];

                g_0_y_xxxx_xxyz[k] = -g_0_y_xxx_xxyz[k] * ab_x + g_0_y_xxx_xxxyz[k];

                g_0_y_xxxx_xxzz[k] = -g_0_y_xxx_xxzz[k] * ab_x + g_0_y_xxx_xxxzz[k];

                g_0_y_xxxx_xyyy[k] = -g_0_y_xxx_xyyy[k] * ab_x + g_0_y_xxx_xxyyy[k];

                g_0_y_xxxx_xyyz[k] = -g_0_y_xxx_xyyz[k] * ab_x + g_0_y_xxx_xxyyz[k];

                g_0_y_xxxx_xyzz[k] = -g_0_y_xxx_xyzz[k] * ab_x + g_0_y_xxx_xxyzz[k];

                g_0_y_xxxx_xzzz[k] = -g_0_y_xxx_xzzz[k] * ab_x + g_0_y_xxx_xxzzz[k];

                g_0_y_xxxx_yyyy[k] = -g_0_y_xxx_yyyy[k] * ab_x + g_0_y_xxx_xyyyy[k];

                g_0_y_xxxx_yyyz[k] = -g_0_y_xxx_yyyz[k] * ab_x + g_0_y_xxx_xyyyz[k];

                g_0_y_xxxx_yyzz[k] = -g_0_y_xxx_yyzz[k] * ab_x + g_0_y_xxx_xyyzz[k];

                g_0_y_xxxx_yzzz[k] = -g_0_y_xxx_yzzz[k] * ab_x + g_0_y_xxx_xyzzz[k];

                g_0_y_xxxx_zzzz[k] = -g_0_y_xxx_zzzz[k] * ab_x + g_0_y_xxx_xzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxy_xxxx = cbuffer.data(gg_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxy = cbuffer.data(gg_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxz = cbuffer.data(gg_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyy = cbuffer.data(gg_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyz = cbuffer.data(gg_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xxxy_xxzz = cbuffer.data(gg_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyy = cbuffer.data(gg_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyz = cbuffer.data(gg_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xxxy_xyzz = cbuffer.data(gg_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xxxy_xzzz = cbuffer.data(gg_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyy = cbuffer.data(gg_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyz = cbuffer.data(gg_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_xxxy_yyzz = cbuffer.data(gg_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xxxy_yzzz = cbuffer.data(gg_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xxxy_zzzz = cbuffer.data(gg_geom_01_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxy_xxxx, g_0_y_xxxy_xxxy, g_0_y_xxxy_xxxz, g_0_y_xxxy_xxyy, g_0_y_xxxy_xxyz, g_0_y_xxxy_xxzz, g_0_y_xxxy_xyyy, g_0_y_xxxy_xyyz, g_0_y_xxxy_xyzz, g_0_y_xxxy_xzzz, g_0_y_xxxy_yyyy, g_0_y_xxxy_yyyz, g_0_y_xxxy_yyzz, g_0_y_xxxy_yzzz, g_0_y_xxxy_zzzz, g_0_y_xxy_xxxx, g_0_y_xxy_xxxxx, g_0_y_xxy_xxxxy, g_0_y_xxy_xxxxz, g_0_y_xxy_xxxy, g_0_y_xxy_xxxyy, g_0_y_xxy_xxxyz, g_0_y_xxy_xxxz, g_0_y_xxy_xxxzz, g_0_y_xxy_xxyy, g_0_y_xxy_xxyyy, g_0_y_xxy_xxyyz, g_0_y_xxy_xxyz, g_0_y_xxy_xxyzz, g_0_y_xxy_xxzz, g_0_y_xxy_xxzzz, g_0_y_xxy_xyyy, g_0_y_xxy_xyyyy, g_0_y_xxy_xyyyz, g_0_y_xxy_xyyz, g_0_y_xxy_xyyzz, g_0_y_xxy_xyzz, g_0_y_xxy_xyzzz, g_0_y_xxy_xzzz, g_0_y_xxy_xzzzz, g_0_y_xxy_yyyy, g_0_y_xxy_yyyz, g_0_y_xxy_yyzz, g_0_y_xxy_yzzz, g_0_y_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxy_xxxx[k] = -g_0_y_xxy_xxxx[k] * ab_x + g_0_y_xxy_xxxxx[k];

                g_0_y_xxxy_xxxy[k] = -g_0_y_xxy_xxxy[k] * ab_x + g_0_y_xxy_xxxxy[k];

                g_0_y_xxxy_xxxz[k] = -g_0_y_xxy_xxxz[k] * ab_x + g_0_y_xxy_xxxxz[k];

                g_0_y_xxxy_xxyy[k] = -g_0_y_xxy_xxyy[k] * ab_x + g_0_y_xxy_xxxyy[k];

                g_0_y_xxxy_xxyz[k] = -g_0_y_xxy_xxyz[k] * ab_x + g_0_y_xxy_xxxyz[k];

                g_0_y_xxxy_xxzz[k] = -g_0_y_xxy_xxzz[k] * ab_x + g_0_y_xxy_xxxzz[k];

                g_0_y_xxxy_xyyy[k] = -g_0_y_xxy_xyyy[k] * ab_x + g_0_y_xxy_xxyyy[k];

                g_0_y_xxxy_xyyz[k] = -g_0_y_xxy_xyyz[k] * ab_x + g_0_y_xxy_xxyyz[k];

                g_0_y_xxxy_xyzz[k] = -g_0_y_xxy_xyzz[k] * ab_x + g_0_y_xxy_xxyzz[k];

                g_0_y_xxxy_xzzz[k] = -g_0_y_xxy_xzzz[k] * ab_x + g_0_y_xxy_xxzzz[k];

                g_0_y_xxxy_yyyy[k] = -g_0_y_xxy_yyyy[k] * ab_x + g_0_y_xxy_xyyyy[k];

                g_0_y_xxxy_yyyz[k] = -g_0_y_xxy_yyyz[k] * ab_x + g_0_y_xxy_xyyyz[k];

                g_0_y_xxxy_yyzz[k] = -g_0_y_xxy_yyzz[k] * ab_x + g_0_y_xxy_xyyzz[k];

                g_0_y_xxxy_yzzz[k] = -g_0_y_xxy_yzzz[k] * ab_x + g_0_y_xxy_xyzzz[k];

                g_0_y_xxxy_zzzz[k] = -g_0_y_xxy_zzzz[k] * ab_x + g_0_y_xxy_xzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxz_xxxx = cbuffer.data(gg_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxy = cbuffer.data(gg_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxz = cbuffer.data(gg_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyy = cbuffer.data(gg_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyz = cbuffer.data(gg_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xxxz_xxzz = cbuffer.data(gg_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyy = cbuffer.data(gg_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyz = cbuffer.data(gg_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xxxz_xyzz = cbuffer.data(gg_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_xxxz_xzzz = cbuffer.data(gg_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyy = cbuffer.data(gg_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyz = cbuffer.data(gg_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xxxz_yyzz = cbuffer.data(gg_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xxxz_yzzz = cbuffer.data(gg_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xxxz_zzzz = cbuffer.data(gg_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxz_xxxx, g_0_y_xxxz_xxxy, g_0_y_xxxz_xxxz, g_0_y_xxxz_xxyy, g_0_y_xxxz_xxyz, g_0_y_xxxz_xxzz, g_0_y_xxxz_xyyy, g_0_y_xxxz_xyyz, g_0_y_xxxz_xyzz, g_0_y_xxxz_xzzz, g_0_y_xxxz_yyyy, g_0_y_xxxz_yyyz, g_0_y_xxxz_yyzz, g_0_y_xxxz_yzzz, g_0_y_xxxz_zzzz, g_0_y_xxz_xxxx, g_0_y_xxz_xxxxx, g_0_y_xxz_xxxxy, g_0_y_xxz_xxxxz, g_0_y_xxz_xxxy, g_0_y_xxz_xxxyy, g_0_y_xxz_xxxyz, g_0_y_xxz_xxxz, g_0_y_xxz_xxxzz, g_0_y_xxz_xxyy, g_0_y_xxz_xxyyy, g_0_y_xxz_xxyyz, g_0_y_xxz_xxyz, g_0_y_xxz_xxyzz, g_0_y_xxz_xxzz, g_0_y_xxz_xxzzz, g_0_y_xxz_xyyy, g_0_y_xxz_xyyyy, g_0_y_xxz_xyyyz, g_0_y_xxz_xyyz, g_0_y_xxz_xyyzz, g_0_y_xxz_xyzz, g_0_y_xxz_xyzzz, g_0_y_xxz_xzzz, g_0_y_xxz_xzzzz, g_0_y_xxz_yyyy, g_0_y_xxz_yyyz, g_0_y_xxz_yyzz, g_0_y_xxz_yzzz, g_0_y_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxz_xxxx[k] = -g_0_y_xxz_xxxx[k] * ab_x + g_0_y_xxz_xxxxx[k];

                g_0_y_xxxz_xxxy[k] = -g_0_y_xxz_xxxy[k] * ab_x + g_0_y_xxz_xxxxy[k];

                g_0_y_xxxz_xxxz[k] = -g_0_y_xxz_xxxz[k] * ab_x + g_0_y_xxz_xxxxz[k];

                g_0_y_xxxz_xxyy[k] = -g_0_y_xxz_xxyy[k] * ab_x + g_0_y_xxz_xxxyy[k];

                g_0_y_xxxz_xxyz[k] = -g_0_y_xxz_xxyz[k] * ab_x + g_0_y_xxz_xxxyz[k];

                g_0_y_xxxz_xxzz[k] = -g_0_y_xxz_xxzz[k] * ab_x + g_0_y_xxz_xxxzz[k];

                g_0_y_xxxz_xyyy[k] = -g_0_y_xxz_xyyy[k] * ab_x + g_0_y_xxz_xxyyy[k];

                g_0_y_xxxz_xyyz[k] = -g_0_y_xxz_xyyz[k] * ab_x + g_0_y_xxz_xxyyz[k];

                g_0_y_xxxz_xyzz[k] = -g_0_y_xxz_xyzz[k] * ab_x + g_0_y_xxz_xxyzz[k];

                g_0_y_xxxz_xzzz[k] = -g_0_y_xxz_xzzz[k] * ab_x + g_0_y_xxz_xxzzz[k];

                g_0_y_xxxz_yyyy[k] = -g_0_y_xxz_yyyy[k] * ab_x + g_0_y_xxz_xyyyy[k];

                g_0_y_xxxz_yyyz[k] = -g_0_y_xxz_yyyz[k] * ab_x + g_0_y_xxz_xyyyz[k];

                g_0_y_xxxz_yyzz[k] = -g_0_y_xxz_yyzz[k] * ab_x + g_0_y_xxz_xyyzz[k];

                g_0_y_xxxz_yzzz[k] = -g_0_y_xxz_yzzz[k] * ab_x + g_0_y_xxz_xyzzz[k];

                g_0_y_xxxz_zzzz[k] = -g_0_y_xxz_zzzz[k] * ab_x + g_0_y_xxz_xzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyy_xxxx = cbuffer.data(gg_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxy = cbuffer.data(gg_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxz = cbuffer.data(gg_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyy = cbuffer.data(gg_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyz = cbuffer.data(gg_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xxyy_xxzz = cbuffer.data(gg_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyy = cbuffer.data(gg_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyz = cbuffer.data(gg_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xxyy_xyzz = cbuffer.data(gg_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xxyy_xzzz = cbuffer.data(gg_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyy = cbuffer.data(gg_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyz = cbuffer.data(gg_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xxyy_yyzz = cbuffer.data(gg_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xxyy_yzzz = cbuffer.data(gg_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xxyy_zzzz = cbuffer.data(gg_geom_01_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyy_xxxx, g_0_y_xxyy_xxxy, g_0_y_xxyy_xxxz, g_0_y_xxyy_xxyy, g_0_y_xxyy_xxyz, g_0_y_xxyy_xxzz, g_0_y_xxyy_xyyy, g_0_y_xxyy_xyyz, g_0_y_xxyy_xyzz, g_0_y_xxyy_xzzz, g_0_y_xxyy_yyyy, g_0_y_xxyy_yyyz, g_0_y_xxyy_yyzz, g_0_y_xxyy_yzzz, g_0_y_xxyy_zzzz, g_0_y_xyy_xxxx, g_0_y_xyy_xxxxx, g_0_y_xyy_xxxxy, g_0_y_xyy_xxxxz, g_0_y_xyy_xxxy, g_0_y_xyy_xxxyy, g_0_y_xyy_xxxyz, g_0_y_xyy_xxxz, g_0_y_xyy_xxxzz, g_0_y_xyy_xxyy, g_0_y_xyy_xxyyy, g_0_y_xyy_xxyyz, g_0_y_xyy_xxyz, g_0_y_xyy_xxyzz, g_0_y_xyy_xxzz, g_0_y_xyy_xxzzz, g_0_y_xyy_xyyy, g_0_y_xyy_xyyyy, g_0_y_xyy_xyyyz, g_0_y_xyy_xyyz, g_0_y_xyy_xyyzz, g_0_y_xyy_xyzz, g_0_y_xyy_xyzzz, g_0_y_xyy_xzzz, g_0_y_xyy_xzzzz, g_0_y_xyy_yyyy, g_0_y_xyy_yyyz, g_0_y_xyy_yyzz, g_0_y_xyy_yzzz, g_0_y_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyy_xxxx[k] = -g_0_y_xyy_xxxx[k] * ab_x + g_0_y_xyy_xxxxx[k];

                g_0_y_xxyy_xxxy[k] = -g_0_y_xyy_xxxy[k] * ab_x + g_0_y_xyy_xxxxy[k];

                g_0_y_xxyy_xxxz[k] = -g_0_y_xyy_xxxz[k] * ab_x + g_0_y_xyy_xxxxz[k];

                g_0_y_xxyy_xxyy[k] = -g_0_y_xyy_xxyy[k] * ab_x + g_0_y_xyy_xxxyy[k];

                g_0_y_xxyy_xxyz[k] = -g_0_y_xyy_xxyz[k] * ab_x + g_0_y_xyy_xxxyz[k];

                g_0_y_xxyy_xxzz[k] = -g_0_y_xyy_xxzz[k] * ab_x + g_0_y_xyy_xxxzz[k];

                g_0_y_xxyy_xyyy[k] = -g_0_y_xyy_xyyy[k] * ab_x + g_0_y_xyy_xxyyy[k];

                g_0_y_xxyy_xyyz[k] = -g_0_y_xyy_xyyz[k] * ab_x + g_0_y_xyy_xxyyz[k];

                g_0_y_xxyy_xyzz[k] = -g_0_y_xyy_xyzz[k] * ab_x + g_0_y_xyy_xxyzz[k];

                g_0_y_xxyy_xzzz[k] = -g_0_y_xyy_xzzz[k] * ab_x + g_0_y_xyy_xxzzz[k];

                g_0_y_xxyy_yyyy[k] = -g_0_y_xyy_yyyy[k] * ab_x + g_0_y_xyy_xyyyy[k];

                g_0_y_xxyy_yyyz[k] = -g_0_y_xyy_yyyz[k] * ab_x + g_0_y_xyy_xyyyz[k];

                g_0_y_xxyy_yyzz[k] = -g_0_y_xyy_yyzz[k] * ab_x + g_0_y_xyy_xyyzz[k];

                g_0_y_xxyy_yzzz[k] = -g_0_y_xyy_yzzz[k] * ab_x + g_0_y_xyy_xyzzz[k];

                g_0_y_xxyy_zzzz[k] = -g_0_y_xyy_zzzz[k] * ab_x + g_0_y_xyy_xzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyz_xxxx = cbuffer.data(gg_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxy = cbuffer.data(gg_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxz = cbuffer.data(gg_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyy = cbuffer.data(gg_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyz = cbuffer.data(gg_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xxyz_xxzz = cbuffer.data(gg_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyy = cbuffer.data(gg_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyz = cbuffer.data(gg_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xxyz_xyzz = cbuffer.data(gg_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xxyz_xzzz = cbuffer.data(gg_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyy = cbuffer.data(gg_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyz = cbuffer.data(gg_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xxyz_yyzz = cbuffer.data(gg_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xxyz_yzzz = cbuffer.data(gg_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xxyz_zzzz = cbuffer.data(gg_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyz_xxxx, g_0_y_xxyz_xxxy, g_0_y_xxyz_xxxz, g_0_y_xxyz_xxyy, g_0_y_xxyz_xxyz, g_0_y_xxyz_xxzz, g_0_y_xxyz_xyyy, g_0_y_xxyz_xyyz, g_0_y_xxyz_xyzz, g_0_y_xxyz_xzzz, g_0_y_xxyz_yyyy, g_0_y_xxyz_yyyz, g_0_y_xxyz_yyzz, g_0_y_xxyz_yzzz, g_0_y_xxyz_zzzz, g_0_y_xyz_xxxx, g_0_y_xyz_xxxxx, g_0_y_xyz_xxxxy, g_0_y_xyz_xxxxz, g_0_y_xyz_xxxy, g_0_y_xyz_xxxyy, g_0_y_xyz_xxxyz, g_0_y_xyz_xxxz, g_0_y_xyz_xxxzz, g_0_y_xyz_xxyy, g_0_y_xyz_xxyyy, g_0_y_xyz_xxyyz, g_0_y_xyz_xxyz, g_0_y_xyz_xxyzz, g_0_y_xyz_xxzz, g_0_y_xyz_xxzzz, g_0_y_xyz_xyyy, g_0_y_xyz_xyyyy, g_0_y_xyz_xyyyz, g_0_y_xyz_xyyz, g_0_y_xyz_xyyzz, g_0_y_xyz_xyzz, g_0_y_xyz_xyzzz, g_0_y_xyz_xzzz, g_0_y_xyz_xzzzz, g_0_y_xyz_yyyy, g_0_y_xyz_yyyz, g_0_y_xyz_yyzz, g_0_y_xyz_yzzz, g_0_y_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyz_xxxx[k] = -g_0_y_xyz_xxxx[k] * ab_x + g_0_y_xyz_xxxxx[k];

                g_0_y_xxyz_xxxy[k] = -g_0_y_xyz_xxxy[k] * ab_x + g_0_y_xyz_xxxxy[k];

                g_0_y_xxyz_xxxz[k] = -g_0_y_xyz_xxxz[k] * ab_x + g_0_y_xyz_xxxxz[k];

                g_0_y_xxyz_xxyy[k] = -g_0_y_xyz_xxyy[k] * ab_x + g_0_y_xyz_xxxyy[k];

                g_0_y_xxyz_xxyz[k] = -g_0_y_xyz_xxyz[k] * ab_x + g_0_y_xyz_xxxyz[k];

                g_0_y_xxyz_xxzz[k] = -g_0_y_xyz_xxzz[k] * ab_x + g_0_y_xyz_xxxzz[k];

                g_0_y_xxyz_xyyy[k] = -g_0_y_xyz_xyyy[k] * ab_x + g_0_y_xyz_xxyyy[k];

                g_0_y_xxyz_xyyz[k] = -g_0_y_xyz_xyyz[k] * ab_x + g_0_y_xyz_xxyyz[k];

                g_0_y_xxyz_xyzz[k] = -g_0_y_xyz_xyzz[k] * ab_x + g_0_y_xyz_xxyzz[k];

                g_0_y_xxyz_xzzz[k] = -g_0_y_xyz_xzzz[k] * ab_x + g_0_y_xyz_xxzzz[k];

                g_0_y_xxyz_yyyy[k] = -g_0_y_xyz_yyyy[k] * ab_x + g_0_y_xyz_xyyyy[k];

                g_0_y_xxyz_yyyz[k] = -g_0_y_xyz_yyyz[k] * ab_x + g_0_y_xyz_xyyyz[k];

                g_0_y_xxyz_yyzz[k] = -g_0_y_xyz_yyzz[k] * ab_x + g_0_y_xyz_xyyzz[k];

                g_0_y_xxyz_yzzz[k] = -g_0_y_xyz_yzzz[k] * ab_x + g_0_y_xyz_xyzzz[k];

                g_0_y_xxyz_zzzz[k] = -g_0_y_xyz_zzzz[k] * ab_x + g_0_y_xyz_xzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzz_xxxx = cbuffer.data(gg_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxy = cbuffer.data(gg_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxz = cbuffer.data(gg_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyy = cbuffer.data(gg_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyz = cbuffer.data(gg_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xxzz_xxzz = cbuffer.data(gg_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyy = cbuffer.data(gg_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyz = cbuffer.data(gg_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xxzz_xyzz = cbuffer.data(gg_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xxzz_xzzz = cbuffer.data(gg_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyy = cbuffer.data(gg_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyz = cbuffer.data(gg_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xxzz_yyzz = cbuffer.data(gg_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xxzz_yzzz = cbuffer.data(gg_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xxzz_zzzz = cbuffer.data(gg_geom_01_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzz_xxxx, g_0_y_xxzz_xxxy, g_0_y_xxzz_xxxz, g_0_y_xxzz_xxyy, g_0_y_xxzz_xxyz, g_0_y_xxzz_xxzz, g_0_y_xxzz_xyyy, g_0_y_xxzz_xyyz, g_0_y_xxzz_xyzz, g_0_y_xxzz_xzzz, g_0_y_xxzz_yyyy, g_0_y_xxzz_yyyz, g_0_y_xxzz_yyzz, g_0_y_xxzz_yzzz, g_0_y_xxzz_zzzz, g_0_y_xzz_xxxx, g_0_y_xzz_xxxxx, g_0_y_xzz_xxxxy, g_0_y_xzz_xxxxz, g_0_y_xzz_xxxy, g_0_y_xzz_xxxyy, g_0_y_xzz_xxxyz, g_0_y_xzz_xxxz, g_0_y_xzz_xxxzz, g_0_y_xzz_xxyy, g_0_y_xzz_xxyyy, g_0_y_xzz_xxyyz, g_0_y_xzz_xxyz, g_0_y_xzz_xxyzz, g_0_y_xzz_xxzz, g_0_y_xzz_xxzzz, g_0_y_xzz_xyyy, g_0_y_xzz_xyyyy, g_0_y_xzz_xyyyz, g_0_y_xzz_xyyz, g_0_y_xzz_xyyzz, g_0_y_xzz_xyzz, g_0_y_xzz_xyzzz, g_0_y_xzz_xzzz, g_0_y_xzz_xzzzz, g_0_y_xzz_yyyy, g_0_y_xzz_yyyz, g_0_y_xzz_yyzz, g_0_y_xzz_yzzz, g_0_y_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzz_xxxx[k] = -g_0_y_xzz_xxxx[k] * ab_x + g_0_y_xzz_xxxxx[k];

                g_0_y_xxzz_xxxy[k] = -g_0_y_xzz_xxxy[k] * ab_x + g_0_y_xzz_xxxxy[k];

                g_0_y_xxzz_xxxz[k] = -g_0_y_xzz_xxxz[k] * ab_x + g_0_y_xzz_xxxxz[k];

                g_0_y_xxzz_xxyy[k] = -g_0_y_xzz_xxyy[k] * ab_x + g_0_y_xzz_xxxyy[k];

                g_0_y_xxzz_xxyz[k] = -g_0_y_xzz_xxyz[k] * ab_x + g_0_y_xzz_xxxyz[k];

                g_0_y_xxzz_xxzz[k] = -g_0_y_xzz_xxzz[k] * ab_x + g_0_y_xzz_xxxzz[k];

                g_0_y_xxzz_xyyy[k] = -g_0_y_xzz_xyyy[k] * ab_x + g_0_y_xzz_xxyyy[k];

                g_0_y_xxzz_xyyz[k] = -g_0_y_xzz_xyyz[k] * ab_x + g_0_y_xzz_xxyyz[k];

                g_0_y_xxzz_xyzz[k] = -g_0_y_xzz_xyzz[k] * ab_x + g_0_y_xzz_xxyzz[k];

                g_0_y_xxzz_xzzz[k] = -g_0_y_xzz_xzzz[k] * ab_x + g_0_y_xzz_xxzzz[k];

                g_0_y_xxzz_yyyy[k] = -g_0_y_xzz_yyyy[k] * ab_x + g_0_y_xzz_xyyyy[k];

                g_0_y_xxzz_yyyz[k] = -g_0_y_xzz_yyyz[k] * ab_x + g_0_y_xzz_xyyyz[k];

                g_0_y_xxzz_yyzz[k] = -g_0_y_xzz_yyzz[k] * ab_x + g_0_y_xzz_xyyzz[k];

                g_0_y_xxzz_yzzz[k] = -g_0_y_xzz_yzzz[k] * ab_x + g_0_y_xzz_xyzzz[k];

                g_0_y_xxzz_zzzz[k] = -g_0_y_xzz_zzzz[k] * ab_x + g_0_y_xzz_xzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyy_xxxx = cbuffer.data(gg_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxy = cbuffer.data(gg_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxz = cbuffer.data(gg_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyy = cbuffer.data(gg_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyz = cbuffer.data(gg_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xyyy_xxzz = cbuffer.data(gg_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyy = cbuffer.data(gg_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyz = cbuffer.data(gg_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xyyy_xyzz = cbuffer.data(gg_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xyyy_xzzz = cbuffer.data(gg_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyy = cbuffer.data(gg_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyz = cbuffer.data(gg_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xyyy_yyzz = cbuffer.data(gg_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xyyy_yzzz = cbuffer.data(gg_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xyyy_zzzz = cbuffer.data(gg_geom_01_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyy_xxxx, g_0_y_xyyy_xxxy, g_0_y_xyyy_xxxz, g_0_y_xyyy_xxyy, g_0_y_xyyy_xxyz, g_0_y_xyyy_xxzz, g_0_y_xyyy_xyyy, g_0_y_xyyy_xyyz, g_0_y_xyyy_xyzz, g_0_y_xyyy_xzzz, g_0_y_xyyy_yyyy, g_0_y_xyyy_yyyz, g_0_y_xyyy_yyzz, g_0_y_xyyy_yzzz, g_0_y_xyyy_zzzz, g_0_y_yyy_xxxx, g_0_y_yyy_xxxxx, g_0_y_yyy_xxxxy, g_0_y_yyy_xxxxz, g_0_y_yyy_xxxy, g_0_y_yyy_xxxyy, g_0_y_yyy_xxxyz, g_0_y_yyy_xxxz, g_0_y_yyy_xxxzz, g_0_y_yyy_xxyy, g_0_y_yyy_xxyyy, g_0_y_yyy_xxyyz, g_0_y_yyy_xxyz, g_0_y_yyy_xxyzz, g_0_y_yyy_xxzz, g_0_y_yyy_xxzzz, g_0_y_yyy_xyyy, g_0_y_yyy_xyyyy, g_0_y_yyy_xyyyz, g_0_y_yyy_xyyz, g_0_y_yyy_xyyzz, g_0_y_yyy_xyzz, g_0_y_yyy_xyzzz, g_0_y_yyy_xzzz, g_0_y_yyy_xzzzz, g_0_y_yyy_yyyy, g_0_y_yyy_yyyz, g_0_y_yyy_yyzz, g_0_y_yyy_yzzz, g_0_y_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyy_xxxx[k] = -g_0_y_yyy_xxxx[k] * ab_x + g_0_y_yyy_xxxxx[k];

                g_0_y_xyyy_xxxy[k] = -g_0_y_yyy_xxxy[k] * ab_x + g_0_y_yyy_xxxxy[k];

                g_0_y_xyyy_xxxz[k] = -g_0_y_yyy_xxxz[k] * ab_x + g_0_y_yyy_xxxxz[k];

                g_0_y_xyyy_xxyy[k] = -g_0_y_yyy_xxyy[k] * ab_x + g_0_y_yyy_xxxyy[k];

                g_0_y_xyyy_xxyz[k] = -g_0_y_yyy_xxyz[k] * ab_x + g_0_y_yyy_xxxyz[k];

                g_0_y_xyyy_xxzz[k] = -g_0_y_yyy_xxzz[k] * ab_x + g_0_y_yyy_xxxzz[k];

                g_0_y_xyyy_xyyy[k] = -g_0_y_yyy_xyyy[k] * ab_x + g_0_y_yyy_xxyyy[k];

                g_0_y_xyyy_xyyz[k] = -g_0_y_yyy_xyyz[k] * ab_x + g_0_y_yyy_xxyyz[k];

                g_0_y_xyyy_xyzz[k] = -g_0_y_yyy_xyzz[k] * ab_x + g_0_y_yyy_xxyzz[k];

                g_0_y_xyyy_xzzz[k] = -g_0_y_yyy_xzzz[k] * ab_x + g_0_y_yyy_xxzzz[k];

                g_0_y_xyyy_yyyy[k] = -g_0_y_yyy_yyyy[k] * ab_x + g_0_y_yyy_xyyyy[k];

                g_0_y_xyyy_yyyz[k] = -g_0_y_yyy_yyyz[k] * ab_x + g_0_y_yyy_xyyyz[k];

                g_0_y_xyyy_yyzz[k] = -g_0_y_yyy_yyzz[k] * ab_x + g_0_y_yyy_xyyzz[k];

                g_0_y_xyyy_yzzz[k] = -g_0_y_yyy_yzzz[k] * ab_x + g_0_y_yyy_xyzzz[k];

                g_0_y_xyyy_zzzz[k] = -g_0_y_yyy_zzzz[k] * ab_x + g_0_y_yyy_xzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyz_xxxx = cbuffer.data(gg_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxy = cbuffer.data(gg_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxz = cbuffer.data(gg_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyy = cbuffer.data(gg_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyz = cbuffer.data(gg_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xyyz_xxzz = cbuffer.data(gg_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyy = cbuffer.data(gg_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyz = cbuffer.data(gg_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_xyyz_xyzz = cbuffer.data(gg_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_xyyz_xzzz = cbuffer.data(gg_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyy = cbuffer.data(gg_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyz = cbuffer.data(gg_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_xyyz_yyzz = cbuffer.data(gg_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_xyyz_yzzz = cbuffer.data(gg_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_xyyz_zzzz = cbuffer.data(gg_geom_01_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyz_xxxx, g_0_y_xyyz_xxxy, g_0_y_xyyz_xxxz, g_0_y_xyyz_xxyy, g_0_y_xyyz_xxyz, g_0_y_xyyz_xxzz, g_0_y_xyyz_xyyy, g_0_y_xyyz_xyyz, g_0_y_xyyz_xyzz, g_0_y_xyyz_xzzz, g_0_y_xyyz_yyyy, g_0_y_xyyz_yyyz, g_0_y_xyyz_yyzz, g_0_y_xyyz_yzzz, g_0_y_xyyz_zzzz, g_0_y_yyz_xxxx, g_0_y_yyz_xxxxx, g_0_y_yyz_xxxxy, g_0_y_yyz_xxxxz, g_0_y_yyz_xxxy, g_0_y_yyz_xxxyy, g_0_y_yyz_xxxyz, g_0_y_yyz_xxxz, g_0_y_yyz_xxxzz, g_0_y_yyz_xxyy, g_0_y_yyz_xxyyy, g_0_y_yyz_xxyyz, g_0_y_yyz_xxyz, g_0_y_yyz_xxyzz, g_0_y_yyz_xxzz, g_0_y_yyz_xxzzz, g_0_y_yyz_xyyy, g_0_y_yyz_xyyyy, g_0_y_yyz_xyyyz, g_0_y_yyz_xyyz, g_0_y_yyz_xyyzz, g_0_y_yyz_xyzz, g_0_y_yyz_xyzzz, g_0_y_yyz_xzzz, g_0_y_yyz_xzzzz, g_0_y_yyz_yyyy, g_0_y_yyz_yyyz, g_0_y_yyz_yyzz, g_0_y_yyz_yzzz, g_0_y_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyz_xxxx[k] = -g_0_y_yyz_xxxx[k] * ab_x + g_0_y_yyz_xxxxx[k];

                g_0_y_xyyz_xxxy[k] = -g_0_y_yyz_xxxy[k] * ab_x + g_0_y_yyz_xxxxy[k];

                g_0_y_xyyz_xxxz[k] = -g_0_y_yyz_xxxz[k] * ab_x + g_0_y_yyz_xxxxz[k];

                g_0_y_xyyz_xxyy[k] = -g_0_y_yyz_xxyy[k] * ab_x + g_0_y_yyz_xxxyy[k];

                g_0_y_xyyz_xxyz[k] = -g_0_y_yyz_xxyz[k] * ab_x + g_0_y_yyz_xxxyz[k];

                g_0_y_xyyz_xxzz[k] = -g_0_y_yyz_xxzz[k] * ab_x + g_0_y_yyz_xxxzz[k];

                g_0_y_xyyz_xyyy[k] = -g_0_y_yyz_xyyy[k] * ab_x + g_0_y_yyz_xxyyy[k];

                g_0_y_xyyz_xyyz[k] = -g_0_y_yyz_xyyz[k] * ab_x + g_0_y_yyz_xxyyz[k];

                g_0_y_xyyz_xyzz[k] = -g_0_y_yyz_xyzz[k] * ab_x + g_0_y_yyz_xxyzz[k];

                g_0_y_xyyz_xzzz[k] = -g_0_y_yyz_xzzz[k] * ab_x + g_0_y_yyz_xxzzz[k];

                g_0_y_xyyz_yyyy[k] = -g_0_y_yyz_yyyy[k] * ab_x + g_0_y_yyz_xyyyy[k];

                g_0_y_xyyz_yyyz[k] = -g_0_y_yyz_yyyz[k] * ab_x + g_0_y_yyz_xyyyz[k];

                g_0_y_xyyz_yyzz[k] = -g_0_y_yyz_yyzz[k] * ab_x + g_0_y_yyz_xyyzz[k];

                g_0_y_xyyz_yzzz[k] = -g_0_y_yyz_yzzz[k] * ab_x + g_0_y_yyz_xyzzz[k];

                g_0_y_xyyz_zzzz[k] = -g_0_y_yyz_zzzz[k] * ab_x + g_0_y_yyz_xzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzz_xxxx = cbuffer.data(gg_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxy = cbuffer.data(gg_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxz = cbuffer.data(gg_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyy = cbuffer.data(gg_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyz = cbuffer.data(gg_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_xyzz_xxzz = cbuffer.data(gg_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyy = cbuffer.data(gg_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyz = cbuffer.data(gg_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_xyzz_xyzz = cbuffer.data(gg_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_xyzz_xzzz = cbuffer.data(gg_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyy = cbuffer.data(gg_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyz = cbuffer.data(gg_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_xyzz_yyzz = cbuffer.data(gg_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_xyzz_yzzz = cbuffer.data(gg_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_xyzz_zzzz = cbuffer.data(gg_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzz_xxxx, g_0_y_xyzz_xxxy, g_0_y_xyzz_xxxz, g_0_y_xyzz_xxyy, g_0_y_xyzz_xxyz, g_0_y_xyzz_xxzz, g_0_y_xyzz_xyyy, g_0_y_xyzz_xyyz, g_0_y_xyzz_xyzz, g_0_y_xyzz_xzzz, g_0_y_xyzz_yyyy, g_0_y_xyzz_yyyz, g_0_y_xyzz_yyzz, g_0_y_xyzz_yzzz, g_0_y_xyzz_zzzz, g_0_y_yzz_xxxx, g_0_y_yzz_xxxxx, g_0_y_yzz_xxxxy, g_0_y_yzz_xxxxz, g_0_y_yzz_xxxy, g_0_y_yzz_xxxyy, g_0_y_yzz_xxxyz, g_0_y_yzz_xxxz, g_0_y_yzz_xxxzz, g_0_y_yzz_xxyy, g_0_y_yzz_xxyyy, g_0_y_yzz_xxyyz, g_0_y_yzz_xxyz, g_0_y_yzz_xxyzz, g_0_y_yzz_xxzz, g_0_y_yzz_xxzzz, g_0_y_yzz_xyyy, g_0_y_yzz_xyyyy, g_0_y_yzz_xyyyz, g_0_y_yzz_xyyz, g_0_y_yzz_xyyzz, g_0_y_yzz_xyzz, g_0_y_yzz_xyzzz, g_0_y_yzz_xzzz, g_0_y_yzz_xzzzz, g_0_y_yzz_yyyy, g_0_y_yzz_yyyz, g_0_y_yzz_yyzz, g_0_y_yzz_yzzz, g_0_y_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzz_xxxx[k] = -g_0_y_yzz_xxxx[k] * ab_x + g_0_y_yzz_xxxxx[k];

                g_0_y_xyzz_xxxy[k] = -g_0_y_yzz_xxxy[k] * ab_x + g_0_y_yzz_xxxxy[k];

                g_0_y_xyzz_xxxz[k] = -g_0_y_yzz_xxxz[k] * ab_x + g_0_y_yzz_xxxxz[k];

                g_0_y_xyzz_xxyy[k] = -g_0_y_yzz_xxyy[k] * ab_x + g_0_y_yzz_xxxyy[k];

                g_0_y_xyzz_xxyz[k] = -g_0_y_yzz_xxyz[k] * ab_x + g_0_y_yzz_xxxyz[k];

                g_0_y_xyzz_xxzz[k] = -g_0_y_yzz_xxzz[k] * ab_x + g_0_y_yzz_xxxzz[k];

                g_0_y_xyzz_xyyy[k] = -g_0_y_yzz_xyyy[k] * ab_x + g_0_y_yzz_xxyyy[k];

                g_0_y_xyzz_xyyz[k] = -g_0_y_yzz_xyyz[k] * ab_x + g_0_y_yzz_xxyyz[k];

                g_0_y_xyzz_xyzz[k] = -g_0_y_yzz_xyzz[k] * ab_x + g_0_y_yzz_xxyzz[k];

                g_0_y_xyzz_xzzz[k] = -g_0_y_yzz_xzzz[k] * ab_x + g_0_y_yzz_xxzzz[k];

                g_0_y_xyzz_yyyy[k] = -g_0_y_yzz_yyyy[k] * ab_x + g_0_y_yzz_xyyyy[k];

                g_0_y_xyzz_yyyz[k] = -g_0_y_yzz_yyyz[k] * ab_x + g_0_y_yzz_xyyyz[k];

                g_0_y_xyzz_yyzz[k] = -g_0_y_yzz_yyzz[k] * ab_x + g_0_y_yzz_xyyzz[k];

                g_0_y_xyzz_yzzz[k] = -g_0_y_yzz_yzzz[k] * ab_x + g_0_y_yzz_xyzzz[k];

                g_0_y_xyzz_zzzz[k] = -g_0_y_yzz_zzzz[k] * ab_x + g_0_y_yzz_xzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzz_xxxx = cbuffer.data(gg_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxy = cbuffer.data(gg_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxz = cbuffer.data(gg_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyy = cbuffer.data(gg_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyz = cbuffer.data(gg_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xzzz_xxzz = cbuffer.data(gg_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyy = cbuffer.data(gg_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyz = cbuffer.data(gg_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xzzz_xyzz = cbuffer.data(gg_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xzzz_xzzz = cbuffer.data(gg_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyy = cbuffer.data(gg_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyz = cbuffer.data(gg_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_xzzz_yyzz = cbuffer.data(gg_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xzzz_yzzz = cbuffer.data(gg_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xzzz_zzzz = cbuffer.data(gg_geom_01_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzz_xxxx, g_0_y_xzzz_xxxy, g_0_y_xzzz_xxxz, g_0_y_xzzz_xxyy, g_0_y_xzzz_xxyz, g_0_y_xzzz_xxzz, g_0_y_xzzz_xyyy, g_0_y_xzzz_xyyz, g_0_y_xzzz_xyzz, g_0_y_xzzz_xzzz, g_0_y_xzzz_yyyy, g_0_y_xzzz_yyyz, g_0_y_xzzz_yyzz, g_0_y_xzzz_yzzz, g_0_y_xzzz_zzzz, g_0_y_zzz_xxxx, g_0_y_zzz_xxxxx, g_0_y_zzz_xxxxy, g_0_y_zzz_xxxxz, g_0_y_zzz_xxxy, g_0_y_zzz_xxxyy, g_0_y_zzz_xxxyz, g_0_y_zzz_xxxz, g_0_y_zzz_xxxzz, g_0_y_zzz_xxyy, g_0_y_zzz_xxyyy, g_0_y_zzz_xxyyz, g_0_y_zzz_xxyz, g_0_y_zzz_xxyzz, g_0_y_zzz_xxzz, g_0_y_zzz_xxzzz, g_0_y_zzz_xyyy, g_0_y_zzz_xyyyy, g_0_y_zzz_xyyyz, g_0_y_zzz_xyyz, g_0_y_zzz_xyyzz, g_0_y_zzz_xyzz, g_0_y_zzz_xyzzz, g_0_y_zzz_xzzz, g_0_y_zzz_xzzzz, g_0_y_zzz_yyyy, g_0_y_zzz_yyyz, g_0_y_zzz_yyzz, g_0_y_zzz_yzzz, g_0_y_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzz_xxxx[k] = -g_0_y_zzz_xxxx[k] * ab_x + g_0_y_zzz_xxxxx[k];

                g_0_y_xzzz_xxxy[k] = -g_0_y_zzz_xxxy[k] * ab_x + g_0_y_zzz_xxxxy[k];

                g_0_y_xzzz_xxxz[k] = -g_0_y_zzz_xxxz[k] * ab_x + g_0_y_zzz_xxxxz[k];

                g_0_y_xzzz_xxyy[k] = -g_0_y_zzz_xxyy[k] * ab_x + g_0_y_zzz_xxxyy[k];

                g_0_y_xzzz_xxyz[k] = -g_0_y_zzz_xxyz[k] * ab_x + g_0_y_zzz_xxxyz[k];

                g_0_y_xzzz_xxzz[k] = -g_0_y_zzz_xxzz[k] * ab_x + g_0_y_zzz_xxxzz[k];

                g_0_y_xzzz_xyyy[k] = -g_0_y_zzz_xyyy[k] * ab_x + g_0_y_zzz_xxyyy[k];

                g_0_y_xzzz_xyyz[k] = -g_0_y_zzz_xyyz[k] * ab_x + g_0_y_zzz_xxyyz[k];

                g_0_y_xzzz_xyzz[k] = -g_0_y_zzz_xyzz[k] * ab_x + g_0_y_zzz_xxyzz[k];

                g_0_y_xzzz_xzzz[k] = -g_0_y_zzz_xzzz[k] * ab_x + g_0_y_zzz_xxzzz[k];

                g_0_y_xzzz_yyyy[k] = -g_0_y_zzz_yyyy[k] * ab_x + g_0_y_zzz_xyyyy[k];

                g_0_y_xzzz_yyyz[k] = -g_0_y_zzz_yyyz[k] * ab_x + g_0_y_zzz_xyyyz[k];

                g_0_y_xzzz_yyzz[k] = -g_0_y_zzz_yyzz[k] * ab_x + g_0_y_zzz_xyyzz[k];

                g_0_y_xzzz_yzzz[k] = -g_0_y_zzz_yzzz[k] * ab_x + g_0_y_zzz_xyzzz[k];

                g_0_y_xzzz_zzzz[k] = -g_0_y_zzz_zzzz[k] * ab_x + g_0_y_zzz_xzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyy_xxxx = cbuffer.data(gg_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxy = cbuffer.data(gg_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxz = cbuffer.data(gg_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyy = cbuffer.data(gg_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyz = cbuffer.data(gg_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_yyyy_xxzz = cbuffer.data(gg_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyy = cbuffer.data(gg_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyz = cbuffer.data(gg_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_yyyy_xyzz = cbuffer.data(gg_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_yyyy_xzzz = cbuffer.data(gg_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyy = cbuffer.data(gg_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyz = cbuffer.data(gg_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_yyyy_yyzz = cbuffer.data(gg_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_yyyy_yzzz = cbuffer.data(gg_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_yyyy_zzzz = cbuffer.data(gg_geom_01_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_xxxx, g_0_y_yyy_xxxxy, g_0_y_yyy_xxxy, g_0_y_yyy_xxxyy, g_0_y_yyy_xxxyz, g_0_y_yyy_xxxz, g_0_y_yyy_xxyy, g_0_y_yyy_xxyyy, g_0_y_yyy_xxyyz, g_0_y_yyy_xxyz, g_0_y_yyy_xxyzz, g_0_y_yyy_xxzz, g_0_y_yyy_xyyy, g_0_y_yyy_xyyyy, g_0_y_yyy_xyyyz, g_0_y_yyy_xyyz, g_0_y_yyy_xyyzz, g_0_y_yyy_xyzz, g_0_y_yyy_xyzzz, g_0_y_yyy_xzzz, g_0_y_yyy_yyyy, g_0_y_yyy_yyyyy, g_0_y_yyy_yyyyz, g_0_y_yyy_yyyz, g_0_y_yyy_yyyzz, g_0_y_yyy_yyzz, g_0_y_yyy_yyzzz, g_0_y_yyy_yzzz, g_0_y_yyy_yzzzz, g_0_y_yyy_zzzz, g_0_y_yyyy_xxxx, g_0_y_yyyy_xxxy, g_0_y_yyyy_xxxz, g_0_y_yyyy_xxyy, g_0_y_yyyy_xxyz, g_0_y_yyyy_xxzz, g_0_y_yyyy_xyyy, g_0_y_yyyy_xyyz, g_0_y_yyyy_xyzz, g_0_y_yyyy_xzzz, g_0_y_yyyy_yyyy, g_0_y_yyyy_yyyz, g_0_y_yyyy_yyzz, g_0_y_yyyy_yzzz, g_0_y_yyyy_zzzz, g_yyy_xxxx, g_yyy_xxxy, g_yyy_xxxz, g_yyy_xxyy, g_yyy_xxyz, g_yyy_xxzz, g_yyy_xyyy, g_yyy_xyyz, g_yyy_xyzz, g_yyy_xzzz, g_yyy_yyyy, g_yyy_yyyz, g_yyy_yyzz, g_yyy_yzzz, g_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyy_xxxx[k] = g_yyy_xxxx[k] - g_0_y_yyy_xxxx[k] * ab_y + g_0_y_yyy_xxxxy[k];

                g_0_y_yyyy_xxxy[k] = g_yyy_xxxy[k] - g_0_y_yyy_xxxy[k] * ab_y + g_0_y_yyy_xxxyy[k];

                g_0_y_yyyy_xxxz[k] = g_yyy_xxxz[k] - g_0_y_yyy_xxxz[k] * ab_y + g_0_y_yyy_xxxyz[k];

                g_0_y_yyyy_xxyy[k] = g_yyy_xxyy[k] - g_0_y_yyy_xxyy[k] * ab_y + g_0_y_yyy_xxyyy[k];

                g_0_y_yyyy_xxyz[k] = g_yyy_xxyz[k] - g_0_y_yyy_xxyz[k] * ab_y + g_0_y_yyy_xxyyz[k];

                g_0_y_yyyy_xxzz[k] = g_yyy_xxzz[k] - g_0_y_yyy_xxzz[k] * ab_y + g_0_y_yyy_xxyzz[k];

                g_0_y_yyyy_xyyy[k] = g_yyy_xyyy[k] - g_0_y_yyy_xyyy[k] * ab_y + g_0_y_yyy_xyyyy[k];

                g_0_y_yyyy_xyyz[k] = g_yyy_xyyz[k] - g_0_y_yyy_xyyz[k] * ab_y + g_0_y_yyy_xyyyz[k];

                g_0_y_yyyy_xyzz[k] = g_yyy_xyzz[k] - g_0_y_yyy_xyzz[k] * ab_y + g_0_y_yyy_xyyzz[k];

                g_0_y_yyyy_xzzz[k] = g_yyy_xzzz[k] - g_0_y_yyy_xzzz[k] * ab_y + g_0_y_yyy_xyzzz[k];

                g_0_y_yyyy_yyyy[k] = g_yyy_yyyy[k] - g_0_y_yyy_yyyy[k] * ab_y + g_0_y_yyy_yyyyy[k];

                g_0_y_yyyy_yyyz[k] = g_yyy_yyyz[k] - g_0_y_yyy_yyyz[k] * ab_y + g_0_y_yyy_yyyyz[k];

                g_0_y_yyyy_yyzz[k] = g_yyy_yyzz[k] - g_0_y_yyy_yyzz[k] * ab_y + g_0_y_yyy_yyyzz[k];

                g_0_y_yyyy_yzzz[k] = g_yyy_yzzz[k] - g_0_y_yyy_yzzz[k] * ab_y + g_0_y_yyy_yyzzz[k];

                g_0_y_yyyy_zzzz[k] = g_yyy_zzzz[k] - g_0_y_yyy_zzzz[k] * ab_y + g_0_y_yyy_yzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyz_xxxx = cbuffer.data(gg_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxy = cbuffer.data(gg_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxz = cbuffer.data(gg_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyy = cbuffer.data(gg_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyz = cbuffer.data(gg_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_yyyz_xxzz = cbuffer.data(gg_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyy = cbuffer.data(gg_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyz = cbuffer.data(gg_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_yyyz_xyzz = cbuffer.data(gg_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_yyyz_xzzz = cbuffer.data(gg_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyy = cbuffer.data(gg_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyz = cbuffer.data(gg_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_yyyz_yyzz = cbuffer.data(gg_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_yyyz_yzzz = cbuffer.data(gg_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_yyyz_zzzz = cbuffer.data(gg_geom_01_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_xxxx, g_0_y_yyy_xxxxz, g_0_y_yyy_xxxy, g_0_y_yyy_xxxyz, g_0_y_yyy_xxxz, g_0_y_yyy_xxxzz, g_0_y_yyy_xxyy, g_0_y_yyy_xxyyz, g_0_y_yyy_xxyz, g_0_y_yyy_xxyzz, g_0_y_yyy_xxzz, g_0_y_yyy_xxzzz, g_0_y_yyy_xyyy, g_0_y_yyy_xyyyz, g_0_y_yyy_xyyz, g_0_y_yyy_xyyzz, g_0_y_yyy_xyzz, g_0_y_yyy_xyzzz, g_0_y_yyy_xzzz, g_0_y_yyy_xzzzz, g_0_y_yyy_yyyy, g_0_y_yyy_yyyyz, g_0_y_yyy_yyyz, g_0_y_yyy_yyyzz, g_0_y_yyy_yyzz, g_0_y_yyy_yyzzz, g_0_y_yyy_yzzz, g_0_y_yyy_yzzzz, g_0_y_yyy_zzzz, g_0_y_yyy_zzzzz, g_0_y_yyyz_xxxx, g_0_y_yyyz_xxxy, g_0_y_yyyz_xxxz, g_0_y_yyyz_xxyy, g_0_y_yyyz_xxyz, g_0_y_yyyz_xxzz, g_0_y_yyyz_xyyy, g_0_y_yyyz_xyyz, g_0_y_yyyz_xyzz, g_0_y_yyyz_xzzz, g_0_y_yyyz_yyyy, g_0_y_yyyz_yyyz, g_0_y_yyyz_yyzz, g_0_y_yyyz_yzzz, g_0_y_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyz_xxxx[k] = -g_0_y_yyy_xxxx[k] * ab_z + g_0_y_yyy_xxxxz[k];

                g_0_y_yyyz_xxxy[k] = -g_0_y_yyy_xxxy[k] * ab_z + g_0_y_yyy_xxxyz[k];

                g_0_y_yyyz_xxxz[k] = -g_0_y_yyy_xxxz[k] * ab_z + g_0_y_yyy_xxxzz[k];

                g_0_y_yyyz_xxyy[k] = -g_0_y_yyy_xxyy[k] * ab_z + g_0_y_yyy_xxyyz[k];

                g_0_y_yyyz_xxyz[k] = -g_0_y_yyy_xxyz[k] * ab_z + g_0_y_yyy_xxyzz[k];

                g_0_y_yyyz_xxzz[k] = -g_0_y_yyy_xxzz[k] * ab_z + g_0_y_yyy_xxzzz[k];

                g_0_y_yyyz_xyyy[k] = -g_0_y_yyy_xyyy[k] * ab_z + g_0_y_yyy_xyyyz[k];

                g_0_y_yyyz_xyyz[k] = -g_0_y_yyy_xyyz[k] * ab_z + g_0_y_yyy_xyyzz[k];

                g_0_y_yyyz_xyzz[k] = -g_0_y_yyy_xyzz[k] * ab_z + g_0_y_yyy_xyzzz[k];

                g_0_y_yyyz_xzzz[k] = -g_0_y_yyy_xzzz[k] * ab_z + g_0_y_yyy_xzzzz[k];

                g_0_y_yyyz_yyyy[k] = -g_0_y_yyy_yyyy[k] * ab_z + g_0_y_yyy_yyyyz[k];

                g_0_y_yyyz_yyyz[k] = -g_0_y_yyy_yyyz[k] * ab_z + g_0_y_yyy_yyyzz[k];

                g_0_y_yyyz_yyzz[k] = -g_0_y_yyy_yyzz[k] * ab_z + g_0_y_yyy_yyzzz[k];

                g_0_y_yyyz_yzzz[k] = -g_0_y_yyy_yzzz[k] * ab_z + g_0_y_yyy_yzzzz[k];

                g_0_y_yyyz_zzzz[k] = -g_0_y_yyy_zzzz[k] * ab_z + g_0_y_yyy_zzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzz_xxxx = cbuffer.data(gg_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxy = cbuffer.data(gg_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxz = cbuffer.data(gg_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyy = cbuffer.data(gg_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyz = cbuffer.data(gg_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_yyzz_xxzz = cbuffer.data(gg_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyy = cbuffer.data(gg_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyz = cbuffer.data(gg_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_yyzz_xyzz = cbuffer.data(gg_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_yyzz_xzzz = cbuffer.data(gg_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyy = cbuffer.data(gg_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyz = cbuffer.data(gg_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_yyzz_yyzz = cbuffer.data(gg_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_yyzz_yzzz = cbuffer.data(gg_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_yyzz_zzzz = cbuffer.data(gg_geom_01_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyz_xxxx, g_0_y_yyz_xxxxz, g_0_y_yyz_xxxy, g_0_y_yyz_xxxyz, g_0_y_yyz_xxxz, g_0_y_yyz_xxxzz, g_0_y_yyz_xxyy, g_0_y_yyz_xxyyz, g_0_y_yyz_xxyz, g_0_y_yyz_xxyzz, g_0_y_yyz_xxzz, g_0_y_yyz_xxzzz, g_0_y_yyz_xyyy, g_0_y_yyz_xyyyz, g_0_y_yyz_xyyz, g_0_y_yyz_xyyzz, g_0_y_yyz_xyzz, g_0_y_yyz_xyzzz, g_0_y_yyz_xzzz, g_0_y_yyz_xzzzz, g_0_y_yyz_yyyy, g_0_y_yyz_yyyyz, g_0_y_yyz_yyyz, g_0_y_yyz_yyyzz, g_0_y_yyz_yyzz, g_0_y_yyz_yyzzz, g_0_y_yyz_yzzz, g_0_y_yyz_yzzzz, g_0_y_yyz_zzzz, g_0_y_yyz_zzzzz, g_0_y_yyzz_xxxx, g_0_y_yyzz_xxxy, g_0_y_yyzz_xxxz, g_0_y_yyzz_xxyy, g_0_y_yyzz_xxyz, g_0_y_yyzz_xxzz, g_0_y_yyzz_xyyy, g_0_y_yyzz_xyyz, g_0_y_yyzz_xyzz, g_0_y_yyzz_xzzz, g_0_y_yyzz_yyyy, g_0_y_yyzz_yyyz, g_0_y_yyzz_yyzz, g_0_y_yyzz_yzzz, g_0_y_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzz_xxxx[k] = -g_0_y_yyz_xxxx[k] * ab_z + g_0_y_yyz_xxxxz[k];

                g_0_y_yyzz_xxxy[k] = -g_0_y_yyz_xxxy[k] * ab_z + g_0_y_yyz_xxxyz[k];

                g_0_y_yyzz_xxxz[k] = -g_0_y_yyz_xxxz[k] * ab_z + g_0_y_yyz_xxxzz[k];

                g_0_y_yyzz_xxyy[k] = -g_0_y_yyz_xxyy[k] * ab_z + g_0_y_yyz_xxyyz[k];

                g_0_y_yyzz_xxyz[k] = -g_0_y_yyz_xxyz[k] * ab_z + g_0_y_yyz_xxyzz[k];

                g_0_y_yyzz_xxzz[k] = -g_0_y_yyz_xxzz[k] * ab_z + g_0_y_yyz_xxzzz[k];

                g_0_y_yyzz_xyyy[k] = -g_0_y_yyz_xyyy[k] * ab_z + g_0_y_yyz_xyyyz[k];

                g_0_y_yyzz_xyyz[k] = -g_0_y_yyz_xyyz[k] * ab_z + g_0_y_yyz_xyyzz[k];

                g_0_y_yyzz_xyzz[k] = -g_0_y_yyz_xyzz[k] * ab_z + g_0_y_yyz_xyzzz[k];

                g_0_y_yyzz_xzzz[k] = -g_0_y_yyz_xzzz[k] * ab_z + g_0_y_yyz_xzzzz[k];

                g_0_y_yyzz_yyyy[k] = -g_0_y_yyz_yyyy[k] * ab_z + g_0_y_yyz_yyyyz[k];

                g_0_y_yyzz_yyyz[k] = -g_0_y_yyz_yyyz[k] * ab_z + g_0_y_yyz_yyyzz[k];

                g_0_y_yyzz_yyzz[k] = -g_0_y_yyz_yyzz[k] * ab_z + g_0_y_yyz_yyzzz[k];

                g_0_y_yyzz_yzzz[k] = -g_0_y_yyz_yzzz[k] * ab_z + g_0_y_yyz_yzzzz[k];

                g_0_y_yyzz_zzzz[k] = -g_0_y_yyz_zzzz[k] * ab_z + g_0_y_yyz_zzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzz_xxxx = cbuffer.data(gg_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxy = cbuffer.data(gg_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxz = cbuffer.data(gg_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyy = cbuffer.data(gg_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyz = cbuffer.data(gg_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_yzzz_xxzz = cbuffer.data(gg_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyy = cbuffer.data(gg_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyz = cbuffer.data(gg_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_yzzz_xyzz = cbuffer.data(gg_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_yzzz_xzzz = cbuffer.data(gg_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyy = cbuffer.data(gg_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyz = cbuffer.data(gg_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_yzzz_yyzz = cbuffer.data(gg_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_yzzz_yzzz = cbuffer.data(gg_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_yzzz_zzzz = cbuffer.data(gg_geom_01_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzz_xxxx, g_0_y_yzz_xxxxz, g_0_y_yzz_xxxy, g_0_y_yzz_xxxyz, g_0_y_yzz_xxxz, g_0_y_yzz_xxxzz, g_0_y_yzz_xxyy, g_0_y_yzz_xxyyz, g_0_y_yzz_xxyz, g_0_y_yzz_xxyzz, g_0_y_yzz_xxzz, g_0_y_yzz_xxzzz, g_0_y_yzz_xyyy, g_0_y_yzz_xyyyz, g_0_y_yzz_xyyz, g_0_y_yzz_xyyzz, g_0_y_yzz_xyzz, g_0_y_yzz_xyzzz, g_0_y_yzz_xzzz, g_0_y_yzz_xzzzz, g_0_y_yzz_yyyy, g_0_y_yzz_yyyyz, g_0_y_yzz_yyyz, g_0_y_yzz_yyyzz, g_0_y_yzz_yyzz, g_0_y_yzz_yyzzz, g_0_y_yzz_yzzz, g_0_y_yzz_yzzzz, g_0_y_yzz_zzzz, g_0_y_yzz_zzzzz, g_0_y_yzzz_xxxx, g_0_y_yzzz_xxxy, g_0_y_yzzz_xxxz, g_0_y_yzzz_xxyy, g_0_y_yzzz_xxyz, g_0_y_yzzz_xxzz, g_0_y_yzzz_xyyy, g_0_y_yzzz_xyyz, g_0_y_yzzz_xyzz, g_0_y_yzzz_xzzz, g_0_y_yzzz_yyyy, g_0_y_yzzz_yyyz, g_0_y_yzzz_yyzz, g_0_y_yzzz_yzzz, g_0_y_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzz_xxxx[k] = -g_0_y_yzz_xxxx[k] * ab_z + g_0_y_yzz_xxxxz[k];

                g_0_y_yzzz_xxxy[k] = -g_0_y_yzz_xxxy[k] * ab_z + g_0_y_yzz_xxxyz[k];

                g_0_y_yzzz_xxxz[k] = -g_0_y_yzz_xxxz[k] * ab_z + g_0_y_yzz_xxxzz[k];

                g_0_y_yzzz_xxyy[k] = -g_0_y_yzz_xxyy[k] * ab_z + g_0_y_yzz_xxyyz[k];

                g_0_y_yzzz_xxyz[k] = -g_0_y_yzz_xxyz[k] * ab_z + g_0_y_yzz_xxyzz[k];

                g_0_y_yzzz_xxzz[k] = -g_0_y_yzz_xxzz[k] * ab_z + g_0_y_yzz_xxzzz[k];

                g_0_y_yzzz_xyyy[k] = -g_0_y_yzz_xyyy[k] * ab_z + g_0_y_yzz_xyyyz[k];

                g_0_y_yzzz_xyyz[k] = -g_0_y_yzz_xyyz[k] * ab_z + g_0_y_yzz_xyyzz[k];

                g_0_y_yzzz_xyzz[k] = -g_0_y_yzz_xyzz[k] * ab_z + g_0_y_yzz_xyzzz[k];

                g_0_y_yzzz_xzzz[k] = -g_0_y_yzz_xzzz[k] * ab_z + g_0_y_yzz_xzzzz[k];

                g_0_y_yzzz_yyyy[k] = -g_0_y_yzz_yyyy[k] * ab_z + g_0_y_yzz_yyyyz[k];

                g_0_y_yzzz_yyyz[k] = -g_0_y_yzz_yyyz[k] * ab_z + g_0_y_yzz_yyyzz[k];

                g_0_y_yzzz_yyzz[k] = -g_0_y_yzz_yyzz[k] * ab_z + g_0_y_yzz_yyzzz[k];

                g_0_y_yzzz_yzzz[k] = -g_0_y_yzz_yzzz[k] * ab_z + g_0_y_yzz_yzzzz[k];

                g_0_y_yzzz_zzzz[k] = -g_0_y_yzz_zzzz[k] * ab_z + g_0_y_yzz_zzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzz_xxxx = cbuffer.data(gg_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxy = cbuffer.data(gg_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxz = cbuffer.data(gg_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyy = cbuffer.data(gg_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyz = cbuffer.data(gg_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_zzzz_xxzz = cbuffer.data(gg_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyy = cbuffer.data(gg_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyz = cbuffer.data(gg_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_zzzz_xyzz = cbuffer.data(gg_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_zzzz_xzzz = cbuffer.data(gg_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyy = cbuffer.data(gg_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyz = cbuffer.data(gg_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_zzzz_yyzz = cbuffer.data(gg_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_zzzz_yzzz = cbuffer.data(gg_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_zzzz_zzzz = cbuffer.data(gg_geom_01_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzz_xxxx, g_0_y_zzz_xxxxz, g_0_y_zzz_xxxy, g_0_y_zzz_xxxyz, g_0_y_zzz_xxxz, g_0_y_zzz_xxxzz, g_0_y_zzz_xxyy, g_0_y_zzz_xxyyz, g_0_y_zzz_xxyz, g_0_y_zzz_xxyzz, g_0_y_zzz_xxzz, g_0_y_zzz_xxzzz, g_0_y_zzz_xyyy, g_0_y_zzz_xyyyz, g_0_y_zzz_xyyz, g_0_y_zzz_xyyzz, g_0_y_zzz_xyzz, g_0_y_zzz_xyzzz, g_0_y_zzz_xzzz, g_0_y_zzz_xzzzz, g_0_y_zzz_yyyy, g_0_y_zzz_yyyyz, g_0_y_zzz_yyyz, g_0_y_zzz_yyyzz, g_0_y_zzz_yyzz, g_0_y_zzz_yyzzz, g_0_y_zzz_yzzz, g_0_y_zzz_yzzzz, g_0_y_zzz_zzzz, g_0_y_zzz_zzzzz, g_0_y_zzzz_xxxx, g_0_y_zzzz_xxxy, g_0_y_zzzz_xxxz, g_0_y_zzzz_xxyy, g_0_y_zzzz_xxyz, g_0_y_zzzz_xxzz, g_0_y_zzzz_xyyy, g_0_y_zzzz_xyyz, g_0_y_zzzz_xyzz, g_0_y_zzzz_xzzz, g_0_y_zzzz_yyyy, g_0_y_zzzz_yyyz, g_0_y_zzzz_yyzz, g_0_y_zzzz_yzzz, g_0_y_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzz_xxxx[k] = -g_0_y_zzz_xxxx[k] * ab_z + g_0_y_zzz_xxxxz[k];

                g_0_y_zzzz_xxxy[k] = -g_0_y_zzz_xxxy[k] * ab_z + g_0_y_zzz_xxxyz[k];

                g_0_y_zzzz_xxxz[k] = -g_0_y_zzz_xxxz[k] * ab_z + g_0_y_zzz_xxxzz[k];

                g_0_y_zzzz_xxyy[k] = -g_0_y_zzz_xxyy[k] * ab_z + g_0_y_zzz_xxyyz[k];

                g_0_y_zzzz_xxyz[k] = -g_0_y_zzz_xxyz[k] * ab_z + g_0_y_zzz_xxyzz[k];

                g_0_y_zzzz_xxzz[k] = -g_0_y_zzz_xxzz[k] * ab_z + g_0_y_zzz_xxzzz[k];

                g_0_y_zzzz_xyyy[k] = -g_0_y_zzz_xyyy[k] * ab_z + g_0_y_zzz_xyyyz[k];

                g_0_y_zzzz_xyyz[k] = -g_0_y_zzz_xyyz[k] * ab_z + g_0_y_zzz_xyyzz[k];

                g_0_y_zzzz_xyzz[k] = -g_0_y_zzz_xyzz[k] * ab_z + g_0_y_zzz_xyzzz[k];

                g_0_y_zzzz_xzzz[k] = -g_0_y_zzz_xzzz[k] * ab_z + g_0_y_zzz_xzzzz[k];

                g_0_y_zzzz_yyyy[k] = -g_0_y_zzz_yyyy[k] * ab_z + g_0_y_zzz_yyyyz[k];

                g_0_y_zzzz_yyyz[k] = -g_0_y_zzz_yyyz[k] * ab_z + g_0_y_zzz_yyyzz[k];

                g_0_y_zzzz_yyzz[k] = -g_0_y_zzz_yyzz[k] * ab_z + g_0_y_zzz_yyzzz[k];

                g_0_y_zzzz_yzzz[k] = -g_0_y_zzz_yzzz[k] * ab_z + g_0_y_zzz_yzzzz[k];

                g_0_y_zzzz_zzzz[k] = -g_0_y_zzz_zzzz[k] * ab_z + g_0_y_zzz_zzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxx_xxxx = cbuffer.data(gg_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxy = cbuffer.data(gg_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxz = cbuffer.data(gg_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyy = cbuffer.data(gg_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyz = cbuffer.data(gg_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xxxx_xxzz = cbuffer.data(gg_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyy = cbuffer.data(gg_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyz = cbuffer.data(gg_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xxxx_xyzz = cbuffer.data(gg_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xxxx_xzzz = cbuffer.data(gg_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyy = cbuffer.data(gg_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyz = cbuffer.data(gg_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_z_xxxx_yyzz = cbuffer.data(gg_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_xxxx_yzzz = cbuffer.data(gg_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_xxxx_zzzz = cbuffer.data(gg_geom_01_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxx_xxxx, g_0_z_xxx_xxxxx, g_0_z_xxx_xxxxy, g_0_z_xxx_xxxxz, g_0_z_xxx_xxxy, g_0_z_xxx_xxxyy, g_0_z_xxx_xxxyz, g_0_z_xxx_xxxz, g_0_z_xxx_xxxzz, g_0_z_xxx_xxyy, g_0_z_xxx_xxyyy, g_0_z_xxx_xxyyz, g_0_z_xxx_xxyz, g_0_z_xxx_xxyzz, g_0_z_xxx_xxzz, g_0_z_xxx_xxzzz, g_0_z_xxx_xyyy, g_0_z_xxx_xyyyy, g_0_z_xxx_xyyyz, g_0_z_xxx_xyyz, g_0_z_xxx_xyyzz, g_0_z_xxx_xyzz, g_0_z_xxx_xyzzz, g_0_z_xxx_xzzz, g_0_z_xxx_xzzzz, g_0_z_xxx_yyyy, g_0_z_xxx_yyyz, g_0_z_xxx_yyzz, g_0_z_xxx_yzzz, g_0_z_xxx_zzzz, g_0_z_xxxx_xxxx, g_0_z_xxxx_xxxy, g_0_z_xxxx_xxxz, g_0_z_xxxx_xxyy, g_0_z_xxxx_xxyz, g_0_z_xxxx_xxzz, g_0_z_xxxx_xyyy, g_0_z_xxxx_xyyz, g_0_z_xxxx_xyzz, g_0_z_xxxx_xzzz, g_0_z_xxxx_yyyy, g_0_z_xxxx_yyyz, g_0_z_xxxx_yyzz, g_0_z_xxxx_yzzz, g_0_z_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxx_xxxx[k] = -g_0_z_xxx_xxxx[k] * ab_x + g_0_z_xxx_xxxxx[k];

                g_0_z_xxxx_xxxy[k] = -g_0_z_xxx_xxxy[k] * ab_x + g_0_z_xxx_xxxxy[k];

                g_0_z_xxxx_xxxz[k] = -g_0_z_xxx_xxxz[k] * ab_x + g_0_z_xxx_xxxxz[k];

                g_0_z_xxxx_xxyy[k] = -g_0_z_xxx_xxyy[k] * ab_x + g_0_z_xxx_xxxyy[k];

                g_0_z_xxxx_xxyz[k] = -g_0_z_xxx_xxyz[k] * ab_x + g_0_z_xxx_xxxyz[k];

                g_0_z_xxxx_xxzz[k] = -g_0_z_xxx_xxzz[k] * ab_x + g_0_z_xxx_xxxzz[k];

                g_0_z_xxxx_xyyy[k] = -g_0_z_xxx_xyyy[k] * ab_x + g_0_z_xxx_xxyyy[k];

                g_0_z_xxxx_xyyz[k] = -g_0_z_xxx_xyyz[k] * ab_x + g_0_z_xxx_xxyyz[k];

                g_0_z_xxxx_xyzz[k] = -g_0_z_xxx_xyzz[k] * ab_x + g_0_z_xxx_xxyzz[k];

                g_0_z_xxxx_xzzz[k] = -g_0_z_xxx_xzzz[k] * ab_x + g_0_z_xxx_xxzzz[k];

                g_0_z_xxxx_yyyy[k] = -g_0_z_xxx_yyyy[k] * ab_x + g_0_z_xxx_xyyyy[k];

                g_0_z_xxxx_yyyz[k] = -g_0_z_xxx_yyyz[k] * ab_x + g_0_z_xxx_xyyyz[k];

                g_0_z_xxxx_yyzz[k] = -g_0_z_xxx_yyzz[k] * ab_x + g_0_z_xxx_xyyzz[k];

                g_0_z_xxxx_yzzz[k] = -g_0_z_xxx_yzzz[k] * ab_x + g_0_z_xxx_xyzzz[k];

                g_0_z_xxxx_zzzz[k] = -g_0_z_xxx_zzzz[k] * ab_x + g_0_z_xxx_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxy_xxxx = cbuffer.data(gg_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxy = cbuffer.data(gg_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxz = cbuffer.data(gg_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyy = cbuffer.data(gg_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyz = cbuffer.data(gg_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_xxxy_xxzz = cbuffer.data(gg_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyy = cbuffer.data(gg_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyz = cbuffer.data(gg_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_xxxy_xyzz = cbuffer.data(gg_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_z_xxxy_xzzz = cbuffer.data(gg_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyy = cbuffer.data(gg_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyz = cbuffer.data(gg_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_xxxy_yyzz = cbuffer.data(gg_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_xxxy_yzzz = cbuffer.data(gg_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_xxxy_zzzz = cbuffer.data(gg_geom_01_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxy_xxxx, g_0_z_xxxy_xxxy, g_0_z_xxxy_xxxz, g_0_z_xxxy_xxyy, g_0_z_xxxy_xxyz, g_0_z_xxxy_xxzz, g_0_z_xxxy_xyyy, g_0_z_xxxy_xyyz, g_0_z_xxxy_xyzz, g_0_z_xxxy_xzzz, g_0_z_xxxy_yyyy, g_0_z_xxxy_yyyz, g_0_z_xxxy_yyzz, g_0_z_xxxy_yzzz, g_0_z_xxxy_zzzz, g_0_z_xxy_xxxx, g_0_z_xxy_xxxxx, g_0_z_xxy_xxxxy, g_0_z_xxy_xxxxz, g_0_z_xxy_xxxy, g_0_z_xxy_xxxyy, g_0_z_xxy_xxxyz, g_0_z_xxy_xxxz, g_0_z_xxy_xxxzz, g_0_z_xxy_xxyy, g_0_z_xxy_xxyyy, g_0_z_xxy_xxyyz, g_0_z_xxy_xxyz, g_0_z_xxy_xxyzz, g_0_z_xxy_xxzz, g_0_z_xxy_xxzzz, g_0_z_xxy_xyyy, g_0_z_xxy_xyyyy, g_0_z_xxy_xyyyz, g_0_z_xxy_xyyz, g_0_z_xxy_xyyzz, g_0_z_xxy_xyzz, g_0_z_xxy_xyzzz, g_0_z_xxy_xzzz, g_0_z_xxy_xzzzz, g_0_z_xxy_yyyy, g_0_z_xxy_yyyz, g_0_z_xxy_yyzz, g_0_z_xxy_yzzz, g_0_z_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxy_xxxx[k] = -g_0_z_xxy_xxxx[k] * ab_x + g_0_z_xxy_xxxxx[k];

                g_0_z_xxxy_xxxy[k] = -g_0_z_xxy_xxxy[k] * ab_x + g_0_z_xxy_xxxxy[k];

                g_0_z_xxxy_xxxz[k] = -g_0_z_xxy_xxxz[k] * ab_x + g_0_z_xxy_xxxxz[k];

                g_0_z_xxxy_xxyy[k] = -g_0_z_xxy_xxyy[k] * ab_x + g_0_z_xxy_xxxyy[k];

                g_0_z_xxxy_xxyz[k] = -g_0_z_xxy_xxyz[k] * ab_x + g_0_z_xxy_xxxyz[k];

                g_0_z_xxxy_xxzz[k] = -g_0_z_xxy_xxzz[k] * ab_x + g_0_z_xxy_xxxzz[k];

                g_0_z_xxxy_xyyy[k] = -g_0_z_xxy_xyyy[k] * ab_x + g_0_z_xxy_xxyyy[k];

                g_0_z_xxxy_xyyz[k] = -g_0_z_xxy_xyyz[k] * ab_x + g_0_z_xxy_xxyyz[k];

                g_0_z_xxxy_xyzz[k] = -g_0_z_xxy_xyzz[k] * ab_x + g_0_z_xxy_xxyzz[k];

                g_0_z_xxxy_xzzz[k] = -g_0_z_xxy_xzzz[k] * ab_x + g_0_z_xxy_xxzzz[k];

                g_0_z_xxxy_yyyy[k] = -g_0_z_xxy_yyyy[k] * ab_x + g_0_z_xxy_xyyyy[k];

                g_0_z_xxxy_yyyz[k] = -g_0_z_xxy_yyyz[k] * ab_x + g_0_z_xxy_xyyyz[k];

                g_0_z_xxxy_yyzz[k] = -g_0_z_xxy_yyzz[k] * ab_x + g_0_z_xxy_xyyzz[k];

                g_0_z_xxxy_yzzz[k] = -g_0_z_xxy_yzzz[k] * ab_x + g_0_z_xxy_xyzzz[k];

                g_0_z_xxxy_zzzz[k] = -g_0_z_xxy_zzzz[k] * ab_x + g_0_z_xxy_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxz_xxxx = cbuffer.data(gg_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxy = cbuffer.data(gg_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxz = cbuffer.data(gg_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyy = cbuffer.data(gg_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyz = cbuffer.data(gg_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_xxxz_xxzz = cbuffer.data(gg_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyy = cbuffer.data(gg_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyz = cbuffer.data(gg_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_xxxz_xyzz = cbuffer.data(gg_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_xxxz_xzzz = cbuffer.data(gg_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyy = cbuffer.data(gg_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyz = cbuffer.data(gg_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_z_xxxz_yyzz = cbuffer.data(gg_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_xxxz_yzzz = cbuffer.data(gg_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_xxxz_zzzz = cbuffer.data(gg_geom_01_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxz_xxxx, g_0_z_xxxz_xxxy, g_0_z_xxxz_xxxz, g_0_z_xxxz_xxyy, g_0_z_xxxz_xxyz, g_0_z_xxxz_xxzz, g_0_z_xxxz_xyyy, g_0_z_xxxz_xyyz, g_0_z_xxxz_xyzz, g_0_z_xxxz_xzzz, g_0_z_xxxz_yyyy, g_0_z_xxxz_yyyz, g_0_z_xxxz_yyzz, g_0_z_xxxz_yzzz, g_0_z_xxxz_zzzz, g_0_z_xxz_xxxx, g_0_z_xxz_xxxxx, g_0_z_xxz_xxxxy, g_0_z_xxz_xxxxz, g_0_z_xxz_xxxy, g_0_z_xxz_xxxyy, g_0_z_xxz_xxxyz, g_0_z_xxz_xxxz, g_0_z_xxz_xxxzz, g_0_z_xxz_xxyy, g_0_z_xxz_xxyyy, g_0_z_xxz_xxyyz, g_0_z_xxz_xxyz, g_0_z_xxz_xxyzz, g_0_z_xxz_xxzz, g_0_z_xxz_xxzzz, g_0_z_xxz_xyyy, g_0_z_xxz_xyyyy, g_0_z_xxz_xyyyz, g_0_z_xxz_xyyz, g_0_z_xxz_xyyzz, g_0_z_xxz_xyzz, g_0_z_xxz_xyzzz, g_0_z_xxz_xzzz, g_0_z_xxz_xzzzz, g_0_z_xxz_yyyy, g_0_z_xxz_yyyz, g_0_z_xxz_yyzz, g_0_z_xxz_yzzz, g_0_z_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxz_xxxx[k] = -g_0_z_xxz_xxxx[k] * ab_x + g_0_z_xxz_xxxxx[k];

                g_0_z_xxxz_xxxy[k] = -g_0_z_xxz_xxxy[k] * ab_x + g_0_z_xxz_xxxxy[k];

                g_0_z_xxxz_xxxz[k] = -g_0_z_xxz_xxxz[k] * ab_x + g_0_z_xxz_xxxxz[k];

                g_0_z_xxxz_xxyy[k] = -g_0_z_xxz_xxyy[k] * ab_x + g_0_z_xxz_xxxyy[k];

                g_0_z_xxxz_xxyz[k] = -g_0_z_xxz_xxyz[k] * ab_x + g_0_z_xxz_xxxyz[k];

                g_0_z_xxxz_xxzz[k] = -g_0_z_xxz_xxzz[k] * ab_x + g_0_z_xxz_xxxzz[k];

                g_0_z_xxxz_xyyy[k] = -g_0_z_xxz_xyyy[k] * ab_x + g_0_z_xxz_xxyyy[k];

                g_0_z_xxxz_xyyz[k] = -g_0_z_xxz_xyyz[k] * ab_x + g_0_z_xxz_xxyyz[k];

                g_0_z_xxxz_xyzz[k] = -g_0_z_xxz_xyzz[k] * ab_x + g_0_z_xxz_xxyzz[k];

                g_0_z_xxxz_xzzz[k] = -g_0_z_xxz_xzzz[k] * ab_x + g_0_z_xxz_xxzzz[k];

                g_0_z_xxxz_yyyy[k] = -g_0_z_xxz_yyyy[k] * ab_x + g_0_z_xxz_xyyyy[k];

                g_0_z_xxxz_yyyz[k] = -g_0_z_xxz_yyyz[k] * ab_x + g_0_z_xxz_xyyyz[k];

                g_0_z_xxxz_yyzz[k] = -g_0_z_xxz_yyzz[k] * ab_x + g_0_z_xxz_xyyzz[k];

                g_0_z_xxxz_yzzz[k] = -g_0_z_xxz_yzzz[k] * ab_x + g_0_z_xxz_xyzzz[k];

                g_0_z_xxxz_zzzz[k] = -g_0_z_xxz_zzzz[k] * ab_x + g_0_z_xxz_xzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyy_xxxx = cbuffer.data(gg_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxy = cbuffer.data(gg_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxz = cbuffer.data(gg_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyy = cbuffer.data(gg_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyz = cbuffer.data(gg_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_xxyy_xxzz = cbuffer.data(gg_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyy = cbuffer.data(gg_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyz = cbuffer.data(gg_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_xxyy_xyzz = cbuffer.data(gg_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_z_xxyy_xzzz = cbuffer.data(gg_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyy = cbuffer.data(gg_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyz = cbuffer.data(gg_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_z_xxyy_yyzz = cbuffer.data(gg_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_z_xxyy_yzzz = cbuffer.data(gg_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_z_xxyy_zzzz = cbuffer.data(gg_geom_01_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyy_xxxx, g_0_z_xxyy_xxxy, g_0_z_xxyy_xxxz, g_0_z_xxyy_xxyy, g_0_z_xxyy_xxyz, g_0_z_xxyy_xxzz, g_0_z_xxyy_xyyy, g_0_z_xxyy_xyyz, g_0_z_xxyy_xyzz, g_0_z_xxyy_xzzz, g_0_z_xxyy_yyyy, g_0_z_xxyy_yyyz, g_0_z_xxyy_yyzz, g_0_z_xxyy_yzzz, g_0_z_xxyy_zzzz, g_0_z_xyy_xxxx, g_0_z_xyy_xxxxx, g_0_z_xyy_xxxxy, g_0_z_xyy_xxxxz, g_0_z_xyy_xxxy, g_0_z_xyy_xxxyy, g_0_z_xyy_xxxyz, g_0_z_xyy_xxxz, g_0_z_xyy_xxxzz, g_0_z_xyy_xxyy, g_0_z_xyy_xxyyy, g_0_z_xyy_xxyyz, g_0_z_xyy_xxyz, g_0_z_xyy_xxyzz, g_0_z_xyy_xxzz, g_0_z_xyy_xxzzz, g_0_z_xyy_xyyy, g_0_z_xyy_xyyyy, g_0_z_xyy_xyyyz, g_0_z_xyy_xyyz, g_0_z_xyy_xyyzz, g_0_z_xyy_xyzz, g_0_z_xyy_xyzzz, g_0_z_xyy_xzzz, g_0_z_xyy_xzzzz, g_0_z_xyy_yyyy, g_0_z_xyy_yyyz, g_0_z_xyy_yyzz, g_0_z_xyy_yzzz, g_0_z_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyy_xxxx[k] = -g_0_z_xyy_xxxx[k] * ab_x + g_0_z_xyy_xxxxx[k];

                g_0_z_xxyy_xxxy[k] = -g_0_z_xyy_xxxy[k] * ab_x + g_0_z_xyy_xxxxy[k];

                g_0_z_xxyy_xxxz[k] = -g_0_z_xyy_xxxz[k] * ab_x + g_0_z_xyy_xxxxz[k];

                g_0_z_xxyy_xxyy[k] = -g_0_z_xyy_xxyy[k] * ab_x + g_0_z_xyy_xxxyy[k];

                g_0_z_xxyy_xxyz[k] = -g_0_z_xyy_xxyz[k] * ab_x + g_0_z_xyy_xxxyz[k];

                g_0_z_xxyy_xxzz[k] = -g_0_z_xyy_xxzz[k] * ab_x + g_0_z_xyy_xxxzz[k];

                g_0_z_xxyy_xyyy[k] = -g_0_z_xyy_xyyy[k] * ab_x + g_0_z_xyy_xxyyy[k];

                g_0_z_xxyy_xyyz[k] = -g_0_z_xyy_xyyz[k] * ab_x + g_0_z_xyy_xxyyz[k];

                g_0_z_xxyy_xyzz[k] = -g_0_z_xyy_xyzz[k] * ab_x + g_0_z_xyy_xxyzz[k];

                g_0_z_xxyy_xzzz[k] = -g_0_z_xyy_xzzz[k] * ab_x + g_0_z_xyy_xxzzz[k];

                g_0_z_xxyy_yyyy[k] = -g_0_z_xyy_yyyy[k] * ab_x + g_0_z_xyy_xyyyy[k];

                g_0_z_xxyy_yyyz[k] = -g_0_z_xyy_yyyz[k] * ab_x + g_0_z_xyy_xyyyz[k];

                g_0_z_xxyy_yyzz[k] = -g_0_z_xyy_yyzz[k] * ab_x + g_0_z_xyy_xyyzz[k];

                g_0_z_xxyy_yzzz[k] = -g_0_z_xyy_yzzz[k] * ab_x + g_0_z_xyy_xyzzz[k];

                g_0_z_xxyy_zzzz[k] = -g_0_z_xyy_zzzz[k] * ab_x + g_0_z_xyy_xzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyz_xxxx = cbuffer.data(gg_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxy = cbuffer.data(gg_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxz = cbuffer.data(gg_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyy = cbuffer.data(gg_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyz = cbuffer.data(gg_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_z_xxyz_xxzz = cbuffer.data(gg_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyy = cbuffer.data(gg_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyz = cbuffer.data(gg_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_z_xxyz_xyzz = cbuffer.data(gg_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_z_xxyz_xzzz = cbuffer.data(gg_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyy = cbuffer.data(gg_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyz = cbuffer.data(gg_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_z_xxyz_yyzz = cbuffer.data(gg_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_z_xxyz_yzzz = cbuffer.data(gg_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_z_xxyz_zzzz = cbuffer.data(gg_geom_01_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyz_xxxx, g_0_z_xxyz_xxxy, g_0_z_xxyz_xxxz, g_0_z_xxyz_xxyy, g_0_z_xxyz_xxyz, g_0_z_xxyz_xxzz, g_0_z_xxyz_xyyy, g_0_z_xxyz_xyyz, g_0_z_xxyz_xyzz, g_0_z_xxyz_xzzz, g_0_z_xxyz_yyyy, g_0_z_xxyz_yyyz, g_0_z_xxyz_yyzz, g_0_z_xxyz_yzzz, g_0_z_xxyz_zzzz, g_0_z_xyz_xxxx, g_0_z_xyz_xxxxx, g_0_z_xyz_xxxxy, g_0_z_xyz_xxxxz, g_0_z_xyz_xxxy, g_0_z_xyz_xxxyy, g_0_z_xyz_xxxyz, g_0_z_xyz_xxxz, g_0_z_xyz_xxxzz, g_0_z_xyz_xxyy, g_0_z_xyz_xxyyy, g_0_z_xyz_xxyyz, g_0_z_xyz_xxyz, g_0_z_xyz_xxyzz, g_0_z_xyz_xxzz, g_0_z_xyz_xxzzz, g_0_z_xyz_xyyy, g_0_z_xyz_xyyyy, g_0_z_xyz_xyyyz, g_0_z_xyz_xyyz, g_0_z_xyz_xyyzz, g_0_z_xyz_xyzz, g_0_z_xyz_xyzzz, g_0_z_xyz_xzzz, g_0_z_xyz_xzzzz, g_0_z_xyz_yyyy, g_0_z_xyz_yyyz, g_0_z_xyz_yyzz, g_0_z_xyz_yzzz, g_0_z_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyz_xxxx[k] = -g_0_z_xyz_xxxx[k] * ab_x + g_0_z_xyz_xxxxx[k];

                g_0_z_xxyz_xxxy[k] = -g_0_z_xyz_xxxy[k] * ab_x + g_0_z_xyz_xxxxy[k];

                g_0_z_xxyz_xxxz[k] = -g_0_z_xyz_xxxz[k] * ab_x + g_0_z_xyz_xxxxz[k];

                g_0_z_xxyz_xxyy[k] = -g_0_z_xyz_xxyy[k] * ab_x + g_0_z_xyz_xxxyy[k];

                g_0_z_xxyz_xxyz[k] = -g_0_z_xyz_xxyz[k] * ab_x + g_0_z_xyz_xxxyz[k];

                g_0_z_xxyz_xxzz[k] = -g_0_z_xyz_xxzz[k] * ab_x + g_0_z_xyz_xxxzz[k];

                g_0_z_xxyz_xyyy[k] = -g_0_z_xyz_xyyy[k] * ab_x + g_0_z_xyz_xxyyy[k];

                g_0_z_xxyz_xyyz[k] = -g_0_z_xyz_xyyz[k] * ab_x + g_0_z_xyz_xxyyz[k];

                g_0_z_xxyz_xyzz[k] = -g_0_z_xyz_xyzz[k] * ab_x + g_0_z_xyz_xxyzz[k];

                g_0_z_xxyz_xzzz[k] = -g_0_z_xyz_xzzz[k] * ab_x + g_0_z_xyz_xxzzz[k];

                g_0_z_xxyz_yyyy[k] = -g_0_z_xyz_yyyy[k] * ab_x + g_0_z_xyz_xyyyy[k];

                g_0_z_xxyz_yyyz[k] = -g_0_z_xyz_yyyz[k] * ab_x + g_0_z_xyz_xyyyz[k];

                g_0_z_xxyz_yyzz[k] = -g_0_z_xyz_yyzz[k] * ab_x + g_0_z_xyz_xyyzz[k];

                g_0_z_xxyz_yzzz[k] = -g_0_z_xyz_yzzz[k] * ab_x + g_0_z_xyz_xyzzz[k];

                g_0_z_xxyz_zzzz[k] = -g_0_z_xyz_zzzz[k] * ab_x + g_0_z_xyz_xzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzz_xxxx = cbuffer.data(gg_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxy = cbuffer.data(gg_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxz = cbuffer.data(gg_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyy = cbuffer.data(gg_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyz = cbuffer.data(gg_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_z_xxzz_xxzz = cbuffer.data(gg_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyy = cbuffer.data(gg_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyz = cbuffer.data(gg_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_z_xxzz_xyzz = cbuffer.data(gg_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_z_xxzz_xzzz = cbuffer.data(gg_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyy = cbuffer.data(gg_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyz = cbuffer.data(gg_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_z_xxzz_yyzz = cbuffer.data(gg_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_z_xxzz_yzzz = cbuffer.data(gg_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_z_xxzz_zzzz = cbuffer.data(gg_geom_01_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzz_xxxx, g_0_z_xxzz_xxxy, g_0_z_xxzz_xxxz, g_0_z_xxzz_xxyy, g_0_z_xxzz_xxyz, g_0_z_xxzz_xxzz, g_0_z_xxzz_xyyy, g_0_z_xxzz_xyyz, g_0_z_xxzz_xyzz, g_0_z_xxzz_xzzz, g_0_z_xxzz_yyyy, g_0_z_xxzz_yyyz, g_0_z_xxzz_yyzz, g_0_z_xxzz_yzzz, g_0_z_xxzz_zzzz, g_0_z_xzz_xxxx, g_0_z_xzz_xxxxx, g_0_z_xzz_xxxxy, g_0_z_xzz_xxxxz, g_0_z_xzz_xxxy, g_0_z_xzz_xxxyy, g_0_z_xzz_xxxyz, g_0_z_xzz_xxxz, g_0_z_xzz_xxxzz, g_0_z_xzz_xxyy, g_0_z_xzz_xxyyy, g_0_z_xzz_xxyyz, g_0_z_xzz_xxyz, g_0_z_xzz_xxyzz, g_0_z_xzz_xxzz, g_0_z_xzz_xxzzz, g_0_z_xzz_xyyy, g_0_z_xzz_xyyyy, g_0_z_xzz_xyyyz, g_0_z_xzz_xyyz, g_0_z_xzz_xyyzz, g_0_z_xzz_xyzz, g_0_z_xzz_xyzzz, g_0_z_xzz_xzzz, g_0_z_xzz_xzzzz, g_0_z_xzz_yyyy, g_0_z_xzz_yyyz, g_0_z_xzz_yyzz, g_0_z_xzz_yzzz, g_0_z_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzz_xxxx[k] = -g_0_z_xzz_xxxx[k] * ab_x + g_0_z_xzz_xxxxx[k];

                g_0_z_xxzz_xxxy[k] = -g_0_z_xzz_xxxy[k] * ab_x + g_0_z_xzz_xxxxy[k];

                g_0_z_xxzz_xxxz[k] = -g_0_z_xzz_xxxz[k] * ab_x + g_0_z_xzz_xxxxz[k];

                g_0_z_xxzz_xxyy[k] = -g_0_z_xzz_xxyy[k] * ab_x + g_0_z_xzz_xxxyy[k];

                g_0_z_xxzz_xxyz[k] = -g_0_z_xzz_xxyz[k] * ab_x + g_0_z_xzz_xxxyz[k];

                g_0_z_xxzz_xxzz[k] = -g_0_z_xzz_xxzz[k] * ab_x + g_0_z_xzz_xxxzz[k];

                g_0_z_xxzz_xyyy[k] = -g_0_z_xzz_xyyy[k] * ab_x + g_0_z_xzz_xxyyy[k];

                g_0_z_xxzz_xyyz[k] = -g_0_z_xzz_xyyz[k] * ab_x + g_0_z_xzz_xxyyz[k];

                g_0_z_xxzz_xyzz[k] = -g_0_z_xzz_xyzz[k] * ab_x + g_0_z_xzz_xxyzz[k];

                g_0_z_xxzz_xzzz[k] = -g_0_z_xzz_xzzz[k] * ab_x + g_0_z_xzz_xxzzz[k];

                g_0_z_xxzz_yyyy[k] = -g_0_z_xzz_yyyy[k] * ab_x + g_0_z_xzz_xyyyy[k];

                g_0_z_xxzz_yyyz[k] = -g_0_z_xzz_yyyz[k] * ab_x + g_0_z_xzz_xyyyz[k];

                g_0_z_xxzz_yyzz[k] = -g_0_z_xzz_yyzz[k] * ab_x + g_0_z_xzz_xyyzz[k];

                g_0_z_xxzz_yzzz[k] = -g_0_z_xzz_yzzz[k] * ab_x + g_0_z_xzz_xyzzz[k];

                g_0_z_xxzz_zzzz[k] = -g_0_z_xzz_zzzz[k] * ab_x + g_0_z_xzz_xzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyy_xxxx = cbuffer.data(gg_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxy = cbuffer.data(gg_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxz = cbuffer.data(gg_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyy = cbuffer.data(gg_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyz = cbuffer.data(gg_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_z_xyyy_xxzz = cbuffer.data(gg_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyy = cbuffer.data(gg_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyz = cbuffer.data(gg_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_z_xyyy_xyzz = cbuffer.data(gg_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_z_xyyy_xzzz = cbuffer.data(gg_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyy = cbuffer.data(gg_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyz = cbuffer.data(gg_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_z_xyyy_yyzz = cbuffer.data(gg_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_z_xyyy_yzzz = cbuffer.data(gg_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_z_xyyy_zzzz = cbuffer.data(gg_geom_01_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyy_xxxx, g_0_z_xyyy_xxxy, g_0_z_xyyy_xxxz, g_0_z_xyyy_xxyy, g_0_z_xyyy_xxyz, g_0_z_xyyy_xxzz, g_0_z_xyyy_xyyy, g_0_z_xyyy_xyyz, g_0_z_xyyy_xyzz, g_0_z_xyyy_xzzz, g_0_z_xyyy_yyyy, g_0_z_xyyy_yyyz, g_0_z_xyyy_yyzz, g_0_z_xyyy_yzzz, g_0_z_xyyy_zzzz, g_0_z_yyy_xxxx, g_0_z_yyy_xxxxx, g_0_z_yyy_xxxxy, g_0_z_yyy_xxxxz, g_0_z_yyy_xxxy, g_0_z_yyy_xxxyy, g_0_z_yyy_xxxyz, g_0_z_yyy_xxxz, g_0_z_yyy_xxxzz, g_0_z_yyy_xxyy, g_0_z_yyy_xxyyy, g_0_z_yyy_xxyyz, g_0_z_yyy_xxyz, g_0_z_yyy_xxyzz, g_0_z_yyy_xxzz, g_0_z_yyy_xxzzz, g_0_z_yyy_xyyy, g_0_z_yyy_xyyyy, g_0_z_yyy_xyyyz, g_0_z_yyy_xyyz, g_0_z_yyy_xyyzz, g_0_z_yyy_xyzz, g_0_z_yyy_xyzzz, g_0_z_yyy_xzzz, g_0_z_yyy_xzzzz, g_0_z_yyy_yyyy, g_0_z_yyy_yyyz, g_0_z_yyy_yyzz, g_0_z_yyy_yzzz, g_0_z_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyy_xxxx[k] = -g_0_z_yyy_xxxx[k] * ab_x + g_0_z_yyy_xxxxx[k];

                g_0_z_xyyy_xxxy[k] = -g_0_z_yyy_xxxy[k] * ab_x + g_0_z_yyy_xxxxy[k];

                g_0_z_xyyy_xxxz[k] = -g_0_z_yyy_xxxz[k] * ab_x + g_0_z_yyy_xxxxz[k];

                g_0_z_xyyy_xxyy[k] = -g_0_z_yyy_xxyy[k] * ab_x + g_0_z_yyy_xxxyy[k];

                g_0_z_xyyy_xxyz[k] = -g_0_z_yyy_xxyz[k] * ab_x + g_0_z_yyy_xxxyz[k];

                g_0_z_xyyy_xxzz[k] = -g_0_z_yyy_xxzz[k] * ab_x + g_0_z_yyy_xxxzz[k];

                g_0_z_xyyy_xyyy[k] = -g_0_z_yyy_xyyy[k] * ab_x + g_0_z_yyy_xxyyy[k];

                g_0_z_xyyy_xyyz[k] = -g_0_z_yyy_xyyz[k] * ab_x + g_0_z_yyy_xxyyz[k];

                g_0_z_xyyy_xyzz[k] = -g_0_z_yyy_xyzz[k] * ab_x + g_0_z_yyy_xxyzz[k];

                g_0_z_xyyy_xzzz[k] = -g_0_z_yyy_xzzz[k] * ab_x + g_0_z_yyy_xxzzz[k];

                g_0_z_xyyy_yyyy[k] = -g_0_z_yyy_yyyy[k] * ab_x + g_0_z_yyy_xyyyy[k];

                g_0_z_xyyy_yyyz[k] = -g_0_z_yyy_yyyz[k] * ab_x + g_0_z_yyy_xyyyz[k];

                g_0_z_xyyy_yyzz[k] = -g_0_z_yyy_yyzz[k] * ab_x + g_0_z_yyy_xyyzz[k];

                g_0_z_xyyy_yzzz[k] = -g_0_z_yyy_yzzz[k] * ab_x + g_0_z_yyy_xyzzz[k];

                g_0_z_xyyy_zzzz[k] = -g_0_z_yyy_zzzz[k] * ab_x + g_0_z_yyy_xzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyz_xxxx = cbuffer.data(gg_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxy = cbuffer.data(gg_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxz = cbuffer.data(gg_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyy = cbuffer.data(gg_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyz = cbuffer.data(gg_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_xyyz_xxzz = cbuffer.data(gg_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyy = cbuffer.data(gg_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyz = cbuffer.data(gg_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_xyyz_xyzz = cbuffer.data(gg_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_xyyz_xzzz = cbuffer.data(gg_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyy = cbuffer.data(gg_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyz = cbuffer.data(gg_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_xyyz_yyzz = cbuffer.data(gg_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_xyyz_yzzz = cbuffer.data(gg_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_xyyz_zzzz = cbuffer.data(gg_geom_01_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyz_xxxx, g_0_z_xyyz_xxxy, g_0_z_xyyz_xxxz, g_0_z_xyyz_xxyy, g_0_z_xyyz_xxyz, g_0_z_xyyz_xxzz, g_0_z_xyyz_xyyy, g_0_z_xyyz_xyyz, g_0_z_xyyz_xyzz, g_0_z_xyyz_xzzz, g_0_z_xyyz_yyyy, g_0_z_xyyz_yyyz, g_0_z_xyyz_yyzz, g_0_z_xyyz_yzzz, g_0_z_xyyz_zzzz, g_0_z_yyz_xxxx, g_0_z_yyz_xxxxx, g_0_z_yyz_xxxxy, g_0_z_yyz_xxxxz, g_0_z_yyz_xxxy, g_0_z_yyz_xxxyy, g_0_z_yyz_xxxyz, g_0_z_yyz_xxxz, g_0_z_yyz_xxxzz, g_0_z_yyz_xxyy, g_0_z_yyz_xxyyy, g_0_z_yyz_xxyyz, g_0_z_yyz_xxyz, g_0_z_yyz_xxyzz, g_0_z_yyz_xxzz, g_0_z_yyz_xxzzz, g_0_z_yyz_xyyy, g_0_z_yyz_xyyyy, g_0_z_yyz_xyyyz, g_0_z_yyz_xyyz, g_0_z_yyz_xyyzz, g_0_z_yyz_xyzz, g_0_z_yyz_xyzzz, g_0_z_yyz_xzzz, g_0_z_yyz_xzzzz, g_0_z_yyz_yyyy, g_0_z_yyz_yyyz, g_0_z_yyz_yyzz, g_0_z_yyz_yzzz, g_0_z_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyz_xxxx[k] = -g_0_z_yyz_xxxx[k] * ab_x + g_0_z_yyz_xxxxx[k];

                g_0_z_xyyz_xxxy[k] = -g_0_z_yyz_xxxy[k] * ab_x + g_0_z_yyz_xxxxy[k];

                g_0_z_xyyz_xxxz[k] = -g_0_z_yyz_xxxz[k] * ab_x + g_0_z_yyz_xxxxz[k];

                g_0_z_xyyz_xxyy[k] = -g_0_z_yyz_xxyy[k] * ab_x + g_0_z_yyz_xxxyy[k];

                g_0_z_xyyz_xxyz[k] = -g_0_z_yyz_xxyz[k] * ab_x + g_0_z_yyz_xxxyz[k];

                g_0_z_xyyz_xxzz[k] = -g_0_z_yyz_xxzz[k] * ab_x + g_0_z_yyz_xxxzz[k];

                g_0_z_xyyz_xyyy[k] = -g_0_z_yyz_xyyy[k] * ab_x + g_0_z_yyz_xxyyy[k];

                g_0_z_xyyz_xyyz[k] = -g_0_z_yyz_xyyz[k] * ab_x + g_0_z_yyz_xxyyz[k];

                g_0_z_xyyz_xyzz[k] = -g_0_z_yyz_xyzz[k] * ab_x + g_0_z_yyz_xxyzz[k];

                g_0_z_xyyz_xzzz[k] = -g_0_z_yyz_xzzz[k] * ab_x + g_0_z_yyz_xxzzz[k];

                g_0_z_xyyz_yyyy[k] = -g_0_z_yyz_yyyy[k] * ab_x + g_0_z_yyz_xyyyy[k];

                g_0_z_xyyz_yyyz[k] = -g_0_z_yyz_yyyz[k] * ab_x + g_0_z_yyz_xyyyz[k];

                g_0_z_xyyz_yyzz[k] = -g_0_z_yyz_yyzz[k] * ab_x + g_0_z_yyz_xyyzz[k];

                g_0_z_xyyz_yzzz[k] = -g_0_z_yyz_yzzz[k] * ab_x + g_0_z_yyz_xyzzz[k];

                g_0_z_xyyz_zzzz[k] = -g_0_z_yyz_zzzz[k] * ab_x + g_0_z_yyz_xzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzz_xxxx = cbuffer.data(gg_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxy = cbuffer.data(gg_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxz = cbuffer.data(gg_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyy = cbuffer.data(gg_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyz = cbuffer.data(gg_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_xyzz_xxzz = cbuffer.data(gg_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyy = cbuffer.data(gg_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyz = cbuffer.data(gg_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_xyzz_xyzz = cbuffer.data(gg_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_xyzz_xzzz = cbuffer.data(gg_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyy = cbuffer.data(gg_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyz = cbuffer.data(gg_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_xyzz_yyzz = cbuffer.data(gg_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_xyzz_yzzz = cbuffer.data(gg_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_xyzz_zzzz = cbuffer.data(gg_geom_01_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzz_xxxx, g_0_z_xyzz_xxxy, g_0_z_xyzz_xxxz, g_0_z_xyzz_xxyy, g_0_z_xyzz_xxyz, g_0_z_xyzz_xxzz, g_0_z_xyzz_xyyy, g_0_z_xyzz_xyyz, g_0_z_xyzz_xyzz, g_0_z_xyzz_xzzz, g_0_z_xyzz_yyyy, g_0_z_xyzz_yyyz, g_0_z_xyzz_yyzz, g_0_z_xyzz_yzzz, g_0_z_xyzz_zzzz, g_0_z_yzz_xxxx, g_0_z_yzz_xxxxx, g_0_z_yzz_xxxxy, g_0_z_yzz_xxxxz, g_0_z_yzz_xxxy, g_0_z_yzz_xxxyy, g_0_z_yzz_xxxyz, g_0_z_yzz_xxxz, g_0_z_yzz_xxxzz, g_0_z_yzz_xxyy, g_0_z_yzz_xxyyy, g_0_z_yzz_xxyyz, g_0_z_yzz_xxyz, g_0_z_yzz_xxyzz, g_0_z_yzz_xxzz, g_0_z_yzz_xxzzz, g_0_z_yzz_xyyy, g_0_z_yzz_xyyyy, g_0_z_yzz_xyyyz, g_0_z_yzz_xyyz, g_0_z_yzz_xyyzz, g_0_z_yzz_xyzz, g_0_z_yzz_xyzzz, g_0_z_yzz_xzzz, g_0_z_yzz_xzzzz, g_0_z_yzz_yyyy, g_0_z_yzz_yyyz, g_0_z_yzz_yyzz, g_0_z_yzz_yzzz, g_0_z_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzz_xxxx[k] = -g_0_z_yzz_xxxx[k] * ab_x + g_0_z_yzz_xxxxx[k];

                g_0_z_xyzz_xxxy[k] = -g_0_z_yzz_xxxy[k] * ab_x + g_0_z_yzz_xxxxy[k];

                g_0_z_xyzz_xxxz[k] = -g_0_z_yzz_xxxz[k] * ab_x + g_0_z_yzz_xxxxz[k];

                g_0_z_xyzz_xxyy[k] = -g_0_z_yzz_xxyy[k] * ab_x + g_0_z_yzz_xxxyy[k];

                g_0_z_xyzz_xxyz[k] = -g_0_z_yzz_xxyz[k] * ab_x + g_0_z_yzz_xxxyz[k];

                g_0_z_xyzz_xxzz[k] = -g_0_z_yzz_xxzz[k] * ab_x + g_0_z_yzz_xxxzz[k];

                g_0_z_xyzz_xyyy[k] = -g_0_z_yzz_xyyy[k] * ab_x + g_0_z_yzz_xxyyy[k];

                g_0_z_xyzz_xyyz[k] = -g_0_z_yzz_xyyz[k] * ab_x + g_0_z_yzz_xxyyz[k];

                g_0_z_xyzz_xyzz[k] = -g_0_z_yzz_xyzz[k] * ab_x + g_0_z_yzz_xxyzz[k];

                g_0_z_xyzz_xzzz[k] = -g_0_z_yzz_xzzz[k] * ab_x + g_0_z_yzz_xxzzz[k];

                g_0_z_xyzz_yyyy[k] = -g_0_z_yzz_yyyy[k] * ab_x + g_0_z_yzz_xyyyy[k];

                g_0_z_xyzz_yyyz[k] = -g_0_z_yzz_yyyz[k] * ab_x + g_0_z_yzz_xyyyz[k];

                g_0_z_xyzz_yyzz[k] = -g_0_z_yzz_yyzz[k] * ab_x + g_0_z_yzz_xyyzz[k];

                g_0_z_xyzz_yzzz[k] = -g_0_z_yzz_yzzz[k] * ab_x + g_0_z_yzz_xyzzz[k];

                g_0_z_xyzz_zzzz[k] = -g_0_z_yzz_zzzz[k] * ab_x + g_0_z_yzz_xzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzz_xxxx = cbuffer.data(gg_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxy = cbuffer.data(gg_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxz = cbuffer.data(gg_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyy = cbuffer.data(gg_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyz = cbuffer.data(gg_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_xzzz_xxzz = cbuffer.data(gg_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyy = cbuffer.data(gg_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyz = cbuffer.data(gg_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_xzzz_xyzz = cbuffer.data(gg_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_xzzz_xzzz = cbuffer.data(gg_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyy = cbuffer.data(gg_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyz = cbuffer.data(gg_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_xzzz_yyzz = cbuffer.data(gg_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_xzzz_yzzz = cbuffer.data(gg_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_xzzz_zzzz = cbuffer.data(gg_geom_01_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzz_xxxx, g_0_z_xzzz_xxxy, g_0_z_xzzz_xxxz, g_0_z_xzzz_xxyy, g_0_z_xzzz_xxyz, g_0_z_xzzz_xxzz, g_0_z_xzzz_xyyy, g_0_z_xzzz_xyyz, g_0_z_xzzz_xyzz, g_0_z_xzzz_xzzz, g_0_z_xzzz_yyyy, g_0_z_xzzz_yyyz, g_0_z_xzzz_yyzz, g_0_z_xzzz_yzzz, g_0_z_xzzz_zzzz, g_0_z_zzz_xxxx, g_0_z_zzz_xxxxx, g_0_z_zzz_xxxxy, g_0_z_zzz_xxxxz, g_0_z_zzz_xxxy, g_0_z_zzz_xxxyy, g_0_z_zzz_xxxyz, g_0_z_zzz_xxxz, g_0_z_zzz_xxxzz, g_0_z_zzz_xxyy, g_0_z_zzz_xxyyy, g_0_z_zzz_xxyyz, g_0_z_zzz_xxyz, g_0_z_zzz_xxyzz, g_0_z_zzz_xxzz, g_0_z_zzz_xxzzz, g_0_z_zzz_xyyy, g_0_z_zzz_xyyyy, g_0_z_zzz_xyyyz, g_0_z_zzz_xyyz, g_0_z_zzz_xyyzz, g_0_z_zzz_xyzz, g_0_z_zzz_xyzzz, g_0_z_zzz_xzzz, g_0_z_zzz_xzzzz, g_0_z_zzz_yyyy, g_0_z_zzz_yyyz, g_0_z_zzz_yyzz, g_0_z_zzz_yzzz, g_0_z_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzz_xxxx[k] = -g_0_z_zzz_xxxx[k] * ab_x + g_0_z_zzz_xxxxx[k];

                g_0_z_xzzz_xxxy[k] = -g_0_z_zzz_xxxy[k] * ab_x + g_0_z_zzz_xxxxy[k];

                g_0_z_xzzz_xxxz[k] = -g_0_z_zzz_xxxz[k] * ab_x + g_0_z_zzz_xxxxz[k];

                g_0_z_xzzz_xxyy[k] = -g_0_z_zzz_xxyy[k] * ab_x + g_0_z_zzz_xxxyy[k];

                g_0_z_xzzz_xxyz[k] = -g_0_z_zzz_xxyz[k] * ab_x + g_0_z_zzz_xxxyz[k];

                g_0_z_xzzz_xxzz[k] = -g_0_z_zzz_xxzz[k] * ab_x + g_0_z_zzz_xxxzz[k];

                g_0_z_xzzz_xyyy[k] = -g_0_z_zzz_xyyy[k] * ab_x + g_0_z_zzz_xxyyy[k];

                g_0_z_xzzz_xyyz[k] = -g_0_z_zzz_xyyz[k] * ab_x + g_0_z_zzz_xxyyz[k];

                g_0_z_xzzz_xyzz[k] = -g_0_z_zzz_xyzz[k] * ab_x + g_0_z_zzz_xxyzz[k];

                g_0_z_xzzz_xzzz[k] = -g_0_z_zzz_xzzz[k] * ab_x + g_0_z_zzz_xxzzz[k];

                g_0_z_xzzz_yyyy[k] = -g_0_z_zzz_yyyy[k] * ab_x + g_0_z_zzz_xyyyy[k];

                g_0_z_xzzz_yyyz[k] = -g_0_z_zzz_yyyz[k] * ab_x + g_0_z_zzz_xyyyz[k];

                g_0_z_xzzz_yyzz[k] = -g_0_z_zzz_yyzz[k] * ab_x + g_0_z_zzz_xyyzz[k];

                g_0_z_xzzz_yzzz[k] = -g_0_z_zzz_yzzz[k] * ab_x + g_0_z_zzz_xyzzz[k];

                g_0_z_xzzz_zzzz[k] = -g_0_z_zzz_zzzz[k] * ab_x + g_0_z_zzz_xzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyy_xxxx = cbuffer.data(gg_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxy = cbuffer.data(gg_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxz = cbuffer.data(gg_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyy = cbuffer.data(gg_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyz = cbuffer.data(gg_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_yyyy_xxzz = cbuffer.data(gg_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyy = cbuffer.data(gg_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyz = cbuffer.data(gg_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_yyyy_xyzz = cbuffer.data(gg_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_yyyy_xzzz = cbuffer.data(gg_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyy = cbuffer.data(gg_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyz = cbuffer.data(gg_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_z_yyyy_yyzz = cbuffer.data(gg_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_yyyy_yzzz = cbuffer.data(gg_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_yyyy_zzzz = cbuffer.data(gg_geom_01_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyy_xxxx, g_0_z_yyy_xxxxy, g_0_z_yyy_xxxy, g_0_z_yyy_xxxyy, g_0_z_yyy_xxxyz, g_0_z_yyy_xxxz, g_0_z_yyy_xxyy, g_0_z_yyy_xxyyy, g_0_z_yyy_xxyyz, g_0_z_yyy_xxyz, g_0_z_yyy_xxyzz, g_0_z_yyy_xxzz, g_0_z_yyy_xyyy, g_0_z_yyy_xyyyy, g_0_z_yyy_xyyyz, g_0_z_yyy_xyyz, g_0_z_yyy_xyyzz, g_0_z_yyy_xyzz, g_0_z_yyy_xyzzz, g_0_z_yyy_xzzz, g_0_z_yyy_yyyy, g_0_z_yyy_yyyyy, g_0_z_yyy_yyyyz, g_0_z_yyy_yyyz, g_0_z_yyy_yyyzz, g_0_z_yyy_yyzz, g_0_z_yyy_yyzzz, g_0_z_yyy_yzzz, g_0_z_yyy_yzzzz, g_0_z_yyy_zzzz, g_0_z_yyyy_xxxx, g_0_z_yyyy_xxxy, g_0_z_yyyy_xxxz, g_0_z_yyyy_xxyy, g_0_z_yyyy_xxyz, g_0_z_yyyy_xxzz, g_0_z_yyyy_xyyy, g_0_z_yyyy_xyyz, g_0_z_yyyy_xyzz, g_0_z_yyyy_xzzz, g_0_z_yyyy_yyyy, g_0_z_yyyy_yyyz, g_0_z_yyyy_yyzz, g_0_z_yyyy_yzzz, g_0_z_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyy_xxxx[k] = -g_0_z_yyy_xxxx[k] * ab_y + g_0_z_yyy_xxxxy[k];

                g_0_z_yyyy_xxxy[k] = -g_0_z_yyy_xxxy[k] * ab_y + g_0_z_yyy_xxxyy[k];

                g_0_z_yyyy_xxxz[k] = -g_0_z_yyy_xxxz[k] * ab_y + g_0_z_yyy_xxxyz[k];

                g_0_z_yyyy_xxyy[k] = -g_0_z_yyy_xxyy[k] * ab_y + g_0_z_yyy_xxyyy[k];

                g_0_z_yyyy_xxyz[k] = -g_0_z_yyy_xxyz[k] * ab_y + g_0_z_yyy_xxyyz[k];

                g_0_z_yyyy_xxzz[k] = -g_0_z_yyy_xxzz[k] * ab_y + g_0_z_yyy_xxyzz[k];

                g_0_z_yyyy_xyyy[k] = -g_0_z_yyy_xyyy[k] * ab_y + g_0_z_yyy_xyyyy[k];

                g_0_z_yyyy_xyyz[k] = -g_0_z_yyy_xyyz[k] * ab_y + g_0_z_yyy_xyyyz[k];

                g_0_z_yyyy_xyzz[k] = -g_0_z_yyy_xyzz[k] * ab_y + g_0_z_yyy_xyyzz[k];

                g_0_z_yyyy_xzzz[k] = -g_0_z_yyy_xzzz[k] * ab_y + g_0_z_yyy_xyzzz[k];

                g_0_z_yyyy_yyyy[k] = -g_0_z_yyy_yyyy[k] * ab_y + g_0_z_yyy_yyyyy[k];

                g_0_z_yyyy_yyyz[k] = -g_0_z_yyy_yyyz[k] * ab_y + g_0_z_yyy_yyyyz[k];

                g_0_z_yyyy_yyzz[k] = -g_0_z_yyy_yyzz[k] * ab_y + g_0_z_yyy_yyyzz[k];

                g_0_z_yyyy_yzzz[k] = -g_0_z_yyy_yzzz[k] * ab_y + g_0_z_yyy_yyzzz[k];

                g_0_z_yyyy_zzzz[k] = -g_0_z_yyy_zzzz[k] * ab_y + g_0_z_yyy_yzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyz_xxxx = cbuffer.data(gg_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxy = cbuffer.data(gg_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxz = cbuffer.data(gg_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyy = cbuffer.data(gg_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyz = cbuffer.data(gg_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_yyyz_xxzz = cbuffer.data(gg_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyy = cbuffer.data(gg_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyz = cbuffer.data(gg_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_yyyz_xyzz = cbuffer.data(gg_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_yyyz_xzzz = cbuffer.data(gg_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyy = cbuffer.data(gg_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyz = cbuffer.data(gg_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_yyyz_yyzz = cbuffer.data(gg_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_yyyz_yzzz = cbuffer.data(gg_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_yyyz_zzzz = cbuffer.data(gg_geom_01_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyz_xxxx, g_0_z_yyyz_xxxy, g_0_z_yyyz_xxxz, g_0_z_yyyz_xxyy, g_0_z_yyyz_xxyz, g_0_z_yyyz_xxzz, g_0_z_yyyz_xyyy, g_0_z_yyyz_xyyz, g_0_z_yyyz_xyzz, g_0_z_yyyz_xzzz, g_0_z_yyyz_yyyy, g_0_z_yyyz_yyyz, g_0_z_yyyz_yyzz, g_0_z_yyyz_yzzz, g_0_z_yyyz_zzzz, g_0_z_yyz_xxxx, g_0_z_yyz_xxxxy, g_0_z_yyz_xxxy, g_0_z_yyz_xxxyy, g_0_z_yyz_xxxyz, g_0_z_yyz_xxxz, g_0_z_yyz_xxyy, g_0_z_yyz_xxyyy, g_0_z_yyz_xxyyz, g_0_z_yyz_xxyz, g_0_z_yyz_xxyzz, g_0_z_yyz_xxzz, g_0_z_yyz_xyyy, g_0_z_yyz_xyyyy, g_0_z_yyz_xyyyz, g_0_z_yyz_xyyz, g_0_z_yyz_xyyzz, g_0_z_yyz_xyzz, g_0_z_yyz_xyzzz, g_0_z_yyz_xzzz, g_0_z_yyz_yyyy, g_0_z_yyz_yyyyy, g_0_z_yyz_yyyyz, g_0_z_yyz_yyyz, g_0_z_yyz_yyyzz, g_0_z_yyz_yyzz, g_0_z_yyz_yyzzz, g_0_z_yyz_yzzz, g_0_z_yyz_yzzzz, g_0_z_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyz_xxxx[k] = -g_0_z_yyz_xxxx[k] * ab_y + g_0_z_yyz_xxxxy[k];

                g_0_z_yyyz_xxxy[k] = -g_0_z_yyz_xxxy[k] * ab_y + g_0_z_yyz_xxxyy[k];

                g_0_z_yyyz_xxxz[k] = -g_0_z_yyz_xxxz[k] * ab_y + g_0_z_yyz_xxxyz[k];

                g_0_z_yyyz_xxyy[k] = -g_0_z_yyz_xxyy[k] * ab_y + g_0_z_yyz_xxyyy[k];

                g_0_z_yyyz_xxyz[k] = -g_0_z_yyz_xxyz[k] * ab_y + g_0_z_yyz_xxyyz[k];

                g_0_z_yyyz_xxzz[k] = -g_0_z_yyz_xxzz[k] * ab_y + g_0_z_yyz_xxyzz[k];

                g_0_z_yyyz_xyyy[k] = -g_0_z_yyz_xyyy[k] * ab_y + g_0_z_yyz_xyyyy[k];

                g_0_z_yyyz_xyyz[k] = -g_0_z_yyz_xyyz[k] * ab_y + g_0_z_yyz_xyyyz[k];

                g_0_z_yyyz_xyzz[k] = -g_0_z_yyz_xyzz[k] * ab_y + g_0_z_yyz_xyyzz[k];

                g_0_z_yyyz_xzzz[k] = -g_0_z_yyz_xzzz[k] * ab_y + g_0_z_yyz_xyzzz[k];

                g_0_z_yyyz_yyyy[k] = -g_0_z_yyz_yyyy[k] * ab_y + g_0_z_yyz_yyyyy[k];

                g_0_z_yyyz_yyyz[k] = -g_0_z_yyz_yyyz[k] * ab_y + g_0_z_yyz_yyyyz[k];

                g_0_z_yyyz_yyzz[k] = -g_0_z_yyz_yyzz[k] * ab_y + g_0_z_yyz_yyyzz[k];

                g_0_z_yyyz_yzzz[k] = -g_0_z_yyz_yzzz[k] * ab_y + g_0_z_yyz_yyzzz[k];

                g_0_z_yyyz_zzzz[k] = -g_0_z_yyz_zzzz[k] * ab_y + g_0_z_yyz_yzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzz_xxxx = cbuffer.data(gg_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxy = cbuffer.data(gg_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxz = cbuffer.data(gg_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyy = cbuffer.data(gg_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyz = cbuffer.data(gg_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_yyzz_xxzz = cbuffer.data(gg_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyy = cbuffer.data(gg_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyz = cbuffer.data(gg_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_yyzz_xyzz = cbuffer.data(gg_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_yyzz_xzzz = cbuffer.data(gg_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyy = cbuffer.data(gg_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyz = cbuffer.data(gg_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_z_yyzz_yyzz = cbuffer.data(gg_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_yyzz_yzzz = cbuffer.data(gg_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_yyzz_zzzz = cbuffer.data(gg_geom_01_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzz_xxxx, g_0_z_yyzz_xxxy, g_0_z_yyzz_xxxz, g_0_z_yyzz_xxyy, g_0_z_yyzz_xxyz, g_0_z_yyzz_xxzz, g_0_z_yyzz_xyyy, g_0_z_yyzz_xyyz, g_0_z_yyzz_xyzz, g_0_z_yyzz_xzzz, g_0_z_yyzz_yyyy, g_0_z_yyzz_yyyz, g_0_z_yyzz_yyzz, g_0_z_yyzz_yzzz, g_0_z_yyzz_zzzz, g_0_z_yzz_xxxx, g_0_z_yzz_xxxxy, g_0_z_yzz_xxxy, g_0_z_yzz_xxxyy, g_0_z_yzz_xxxyz, g_0_z_yzz_xxxz, g_0_z_yzz_xxyy, g_0_z_yzz_xxyyy, g_0_z_yzz_xxyyz, g_0_z_yzz_xxyz, g_0_z_yzz_xxyzz, g_0_z_yzz_xxzz, g_0_z_yzz_xyyy, g_0_z_yzz_xyyyy, g_0_z_yzz_xyyyz, g_0_z_yzz_xyyz, g_0_z_yzz_xyyzz, g_0_z_yzz_xyzz, g_0_z_yzz_xyzzz, g_0_z_yzz_xzzz, g_0_z_yzz_yyyy, g_0_z_yzz_yyyyy, g_0_z_yzz_yyyyz, g_0_z_yzz_yyyz, g_0_z_yzz_yyyzz, g_0_z_yzz_yyzz, g_0_z_yzz_yyzzz, g_0_z_yzz_yzzz, g_0_z_yzz_yzzzz, g_0_z_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzz_xxxx[k] = -g_0_z_yzz_xxxx[k] * ab_y + g_0_z_yzz_xxxxy[k];

                g_0_z_yyzz_xxxy[k] = -g_0_z_yzz_xxxy[k] * ab_y + g_0_z_yzz_xxxyy[k];

                g_0_z_yyzz_xxxz[k] = -g_0_z_yzz_xxxz[k] * ab_y + g_0_z_yzz_xxxyz[k];

                g_0_z_yyzz_xxyy[k] = -g_0_z_yzz_xxyy[k] * ab_y + g_0_z_yzz_xxyyy[k];

                g_0_z_yyzz_xxyz[k] = -g_0_z_yzz_xxyz[k] * ab_y + g_0_z_yzz_xxyyz[k];

                g_0_z_yyzz_xxzz[k] = -g_0_z_yzz_xxzz[k] * ab_y + g_0_z_yzz_xxyzz[k];

                g_0_z_yyzz_xyyy[k] = -g_0_z_yzz_xyyy[k] * ab_y + g_0_z_yzz_xyyyy[k];

                g_0_z_yyzz_xyyz[k] = -g_0_z_yzz_xyyz[k] * ab_y + g_0_z_yzz_xyyyz[k];

                g_0_z_yyzz_xyzz[k] = -g_0_z_yzz_xyzz[k] * ab_y + g_0_z_yzz_xyyzz[k];

                g_0_z_yyzz_xzzz[k] = -g_0_z_yzz_xzzz[k] * ab_y + g_0_z_yzz_xyzzz[k];

                g_0_z_yyzz_yyyy[k] = -g_0_z_yzz_yyyy[k] * ab_y + g_0_z_yzz_yyyyy[k];

                g_0_z_yyzz_yyyz[k] = -g_0_z_yzz_yyyz[k] * ab_y + g_0_z_yzz_yyyyz[k];

                g_0_z_yyzz_yyzz[k] = -g_0_z_yzz_yyzz[k] * ab_y + g_0_z_yzz_yyyzz[k];

                g_0_z_yyzz_yzzz[k] = -g_0_z_yzz_yzzz[k] * ab_y + g_0_z_yzz_yyzzz[k];

                g_0_z_yyzz_zzzz[k] = -g_0_z_yzz_zzzz[k] * ab_y + g_0_z_yzz_yzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzz_xxxx = cbuffer.data(gg_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxy = cbuffer.data(gg_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxz = cbuffer.data(gg_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyy = cbuffer.data(gg_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyz = cbuffer.data(gg_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_z_yzzz_xxzz = cbuffer.data(gg_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyy = cbuffer.data(gg_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyz = cbuffer.data(gg_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_z_yzzz_xyzz = cbuffer.data(gg_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_z_yzzz_xzzz = cbuffer.data(gg_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyy = cbuffer.data(gg_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyz = cbuffer.data(gg_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_z_yzzz_yyzz = cbuffer.data(gg_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_z_yzzz_yzzz = cbuffer.data(gg_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_z_yzzz_zzzz = cbuffer.data(gg_geom_01_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzz_xxxx, g_0_z_yzzz_xxxy, g_0_z_yzzz_xxxz, g_0_z_yzzz_xxyy, g_0_z_yzzz_xxyz, g_0_z_yzzz_xxzz, g_0_z_yzzz_xyyy, g_0_z_yzzz_xyyz, g_0_z_yzzz_xyzz, g_0_z_yzzz_xzzz, g_0_z_yzzz_yyyy, g_0_z_yzzz_yyyz, g_0_z_yzzz_yyzz, g_0_z_yzzz_yzzz, g_0_z_yzzz_zzzz, g_0_z_zzz_xxxx, g_0_z_zzz_xxxxy, g_0_z_zzz_xxxy, g_0_z_zzz_xxxyy, g_0_z_zzz_xxxyz, g_0_z_zzz_xxxz, g_0_z_zzz_xxyy, g_0_z_zzz_xxyyy, g_0_z_zzz_xxyyz, g_0_z_zzz_xxyz, g_0_z_zzz_xxyzz, g_0_z_zzz_xxzz, g_0_z_zzz_xyyy, g_0_z_zzz_xyyyy, g_0_z_zzz_xyyyz, g_0_z_zzz_xyyz, g_0_z_zzz_xyyzz, g_0_z_zzz_xyzz, g_0_z_zzz_xyzzz, g_0_z_zzz_xzzz, g_0_z_zzz_yyyy, g_0_z_zzz_yyyyy, g_0_z_zzz_yyyyz, g_0_z_zzz_yyyz, g_0_z_zzz_yyyzz, g_0_z_zzz_yyzz, g_0_z_zzz_yyzzz, g_0_z_zzz_yzzz, g_0_z_zzz_yzzzz, g_0_z_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzz_xxxx[k] = -g_0_z_zzz_xxxx[k] * ab_y + g_0_z_zzz_xxxxy[k];

                g_0_z_yzzz_xxxy[k] = -g_0_z_zzz_xxxy[k] * ab_y + g_0_z_zzz_xxxyy[k];

                g_0_z_yzzz_xxxz[k] = -g_0_z_zzz_xxxz[k] * ab_y + g_0_z_zzz_xxxyz[k];

                g_0_z_yzzz_xxyy[k] = -g_0_z_zzz_xxyy[k] * ab_y + g_0_z_zzz_xxyyy[k];

                g_0_z_yzzz_xxyz[k] = -g_0_z_zzz_xxyz[k] * ab_y + g_0_z_zzz_xxyyz[k];

                g_0_z_yzzz_xxzz[k] = -g_0_z_zzz_xxzz[k] * ab_y + g_0_z_zzz_xxyzz[k];

                g_0_z_yzzz_xyyy[k] = -g_0_z_zzz_xyyy[k] * ab_y + g_0_z_zzz_xyyyy[k];

                g_0_z_yzzz_xyyz[k] = -g_0_z_zzz_xyyz[k] * ab_y + g_0_z_zzz_xyyyz[k];

                g_0_z_yzzz_xyzz[k] = -g_0_z_zzz_xyzz[k] * ab_y + g_0_z_zzz_xyyzz[k];

                g_0_z_yzzz_xzzz[k] = -g_0_z_zzz_xzzz[k] * ab_y + g_0_z_zzz_xyzzz[k];

                g_0_z_yzzz_yyyy[k] = -g_0_z_zzz_yyyy[k] * ab_y + g_0_z_zzz_yyyyy[k];

                g_0_z_yzzz_yyyz[k] = -g_0_z_zzz_yyyz[k] * ab_y + g_0_z_zzz_yyyyz[k];

                g_0_z_yzzz_yyzz[k] = -g_0_z_zzz_yyzz[k] * ab_y + g_0_z_zzz_yyyzz[k];

                g_0_z_yzzz_yzzz[k] = -g_0_z_zzz_yzzz[k] * ab_y + g_0_z_zzz_yyzzz[k];

                g_0_z_yzzz_zzzz[k] = -g_0_z_zzz_zzzz[k] * ab_y + g_0_z_zzz_yzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzz_xxxx = cbuffer.data(gg_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxy = cbuffer.data(gg_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxz = cbuffer.data(gg_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyy = cbuffer.data(gg_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyz = cbuffer.data(gg_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_z_zzzz_xxzz = cbuffer.data(gg_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyy = cbuffer.data(gg_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyz = cbuffer.data(gg_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_z_zzzz_xyzz = cbuffer.data(gg_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_z_zzzz_xzzz = cbuffer.data(gg_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyy = cbuffer.data(gg_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyz = cbuffer.data(gg_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_z_zzzz_yyzz = cbuffer.data(gg_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_z_zzzz_yzzz = cbuffer.data(gg_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_z_zzzz_zzzz = cbuffer.data(gg_geom_01_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_xxxx, g_0_z_zzz_xxxxz, g_0_z_zzz_xxxy, g_0_z_zzz_xxxyz, g_0_z_zzz_xxxz, g_0_z_zzz_xxxzz, g_0_z_zzz_xxyy, g_0_z_zzz_xxyyz, g_0_z_zzz_xxyz, g_0_z_zzz_xxyzz, g_0_z_zzz_xxzz, g_0_z_zzz_xxzzz, g_0_z_zzz_xyyy, g_0_z_zzz_xyyyz, g_0_z_zzz_xyyz, g_0_z_zzz_xyyzz, g_0_z_zzz_xyzz, g_0_z_zzz_xyzzz, g_0_z_zzz_xzzz, g_0_z_zzz_xzzzz, g_0_z_zzz_yyyy, g_0_z_zzz_yyyyz, g_0_z_zzz_yyyz, g_0_z_zzz_yyyzz, g_0_z_zzz_yyzz, g_0_z_zzz_yyzzz, g_0_z_zzz_yzzz, g_0_z_zzz_yzzzz, g_0_z_zzz_zzzz, g_0_z_zzz_zzzzz, g_0_z_zzzz_xxxx, g_0_z_zzzz_xxxy, g_0_z_zzzz_xxxz, g_0_z_zzzz_xxyy, g_0_z_zzzz_xxyz, g_0_z_zzzz_xxzz, g_0_z_zzzz_xyyy, g_0_z_zzzz_xyyz, g_0_z_zzzz_xyzz, g_0_z_zzzz_xzzz, g_0_z_zzzz_yyyy, g_0_z_zzzz_yyyz, g_0_z_zzzz_yyzz, g_0_z_zzzz_yzzz, g_0_z_zzzz_zzzz, g_zzz_xxxx, g_zzz_xxxy, g_zzz_xxxz, g_zzz_xxyy, g_zzz_xxyz, g_zzz_xxzz, g_zzz_xyyy, g_zzz_xyyz, g_zzz_xyzz, g_zzz_xzzz, g_zzz_yyyy, g_zzz_yyyz, g_zzz_yyzz, g_zzz_yzzz, g_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzz_xxxx[k] = g_zzz_xxxx[k] - g_0_z_zzz_xxxx[k] * ab_z + g_0_z_zzz_xxxxz[k];

                g_0_z_zzzz_xxxy[k] = g_zzz_xxxy[k] - g_0_z_zzz_xxxy[k] * ab_z + g_0_z_zzz_xxxyz[k];

                g_0_z_zzzz_xxxz[k] = g_zzz_xxxz[k] - g_0_z_zzz_xxxz[k] * ab_z + g_0_z_zzz_xxxzz[k];

                g_0_z_zzzz_xxyy[k] = g_zzz_xxyy[k] - g_0_z_zzz_xxyy[k] * ab_z + g_0_z_zzz_xxyyz[k];

                g_0_z_zzzz_xxyz[k] = g_zzz_xxyz[k] - g_0_z_zzz_xxyz[k] * ab_z + g_0_z_zzz_xxyzz[k];

                g_0_z_zzzz_xxzz[k] = g_zzz_xxzz[k] - g_0_z_zzz_xxzz[k] * ab_z + g_0_z_zzz_xxzzz[k];

                g_0_z_zzzz_xyyy[k] = g_zzz_xyyy[k] - g_0_z_zzz_xyyy[k] * ab_z + g_0_z_zzz_xyyyz[k];

                g_0_z_zzzz_xyyz[k] = g_zzz_xyyz[k] - g_0_z_zzz_xyyz[k] * ab_z + g_0_z_zzz_xyyzz[k];

                g_0_z_zzzz_xyzz[k] = g_zzz_xyzz[k] - g_0_z_zzz_xyzz[k] * ab_z + g_0_z_zzz_xyzzz[k];

                g_0_z_zzzz_xzzz[k] = g_zzz_xzzz[k] - g_0_z_zzz_xzzz[k] * ab_z + g_0_z_zzz_xzzzz[k];

                g_0_z_zzzz_yyyy[k] = g_zzz_yyyy[k] - g_0_z_zzz_yyyy[k] * ab_z + g_0_z_zzz_yyyyz[k];

                g_0_z_zzzz_yyyz[k] = g_zzz_yyyz[k] - g_0_z_zzz_yyyz[k] * ab_z + g_0_z_zzz_yyyzz[k];

                g_0_z_zzzz_yyzz[k] = g_zzz_yyzz[k] - g_0_z_zzz_yyzz[k] * ab_z + g_0_z_zzz_yyzzz[k];

                g_0_z_zzzz_yzzz[k] = g_zzz_yzzz[k] - g_0_z_zzz_yzzz[k] * ab_z + g_0_z_zzz_yzzzz[k];

                g_0_z_zzzz_zzzz[k] = g_zzz_zzzz[k] - g_0_z_zzz_zzzz[k] * ab_z + g_0_z_zzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

