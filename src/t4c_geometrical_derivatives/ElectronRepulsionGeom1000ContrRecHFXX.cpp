#include "ElectronRepulsionGeom1000ContrRecHFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_hfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_hfxx,
                                            const size_t idx_gfxx,
                                            const size_t idx_geom_10_gfxx,
                                            const size_t idx_geom_10_ggxx,
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
            /// Set up components of auxilary buffer : GFSS

            const auto gf_off = idx_gfxx + i * dcomps + j;

            auto g_xxxx_xxx = cbuffer.data(gf_off + 0 * ccomps * dcomps);

            auto g_xxxx_xxy = cbuffer.data(gf_off + 1 * ccomps * dcomps);

            auto g_xxxx_xxz = cbuffer.data(gf_off + 2 * ccomps * dcomps);

            auto g_xxxx_xyy = cbuffer.data(gf_off + 3 * ccomps * dcomps);

            auto g_xxxx_xyz = cbuffer.data(gf_off + 4 * ccomps * dcomps);

            auto g_xxxx_xzz = cbuffer.data(gf_off + 5 * ccomps * dcomps);

            auto g_xxxx_yyy = cbuffer.data(gf_off + 6 * ccomps * dcomps);

            auto g_xxxx_yyz = cbuffer.data(gf_off + 7 * ccomps * dcomps);

            auto g_xxxx_yzz = cbuffer.data(gf_off + 8 * ccomps * dcomps);

            auto g_xxxx_zzz = cbuffer.data(gf_off + 9 * ccomps * dcomps);

            auto g_xxxy_xxx = cbuffer.data(gf_off + 10 * ccomps * dcomps);

            auto g_xxxy_xxy = cbuffer.data(gf_off + 11 * ccomps * dcomps);

            auto g_xxxy_xxz = cbuffer.data(gf_off + 12 * ccomps * dcomps);

            auto g_xxxy_xyy = cbuffer.data(gf_off + 13 * ccomps * dcomps);

            auto g_xxxy_xyz = cbuffer.data(gf_off + 14 * ccomps * dcomps);

            auto g_xxxy_xzz = cbuffer.data(gf_off + 15 * ccomps * dcomps);

            auto g_xxxy_yyy = cbuffer.data(gf_off + 16 * ccomps * dcomps);

            auto g_xxxy_yyz = cbuffer.data(gf_off + 17 * ccomps * dcomps);

            auto g_xxxy_yzz = cbuffer.data(gf_off + 18 * ccomps * dcomps);

            auto g_xxxy_zzz = cbuffer.data(gf_off + 19 * ccomps * dcomps);

            auto g_xxxz_xxx = cbuffer.data(gf_off + 20 * ccomps * dcomps);

            auto g_xxxz_xxy = cbuffer.data(gf_off + 21 * ccomps * dcomps);

            auto g_xxxz_xxz = cbuffer.data(gf_off + 22 * ccomps * dcomps);

            auto g_xxxz_xyy = cbuffer.data(gf_off + 23 * ccomps * dcomps);

            auto g_xxxz_xyz = cbuffer.data(gf_off + 24 * ccomps * dcomps);

            auto g_xxxz_xzz = cbuffer.data(gf_off + 25 * ccomps * dcomps);

            auto g_xxxz_yyy = cbuffer.data(gf_off + 26 * ccomps * dcomps);

            auto g_xxxz_yyz = cbuffer.data(gf_off + 27 * ccomps * dcomps);

            auto g_xxxz_yzz = cbuffer.data(gf_off + 28 * ccomps * dcomps);

            auto g_xxxz_zzz = cbuffer.data(gf_off + 29 * ccomps * dcomps);

            auto g_xxyy_xxx = cbuffer.data(gf_off + 30 * ccomps * dcomps);

            auto g_xxyy_xxy = cbuffer.data(gf_off + 31 * ccomps * dcomps);

            auto g_xxyy_xxz = cbuffer.data(gf_off + 32 * ccomps * dcomps);

            auto g_xxyy_xyy = cbuffer.data(gf_off + 33 * ccomps * dcomps);

            auto g_xxyy_xyz = cbuffer.data(gf_off + 34 * ccomps * dcomps);

            auto g_xxyy_xzz = cbuffer.data(gf_off + 35 * ccomps * dcomps);

            auto g_xxyy_yyy = cbuffer.data(gf_off + 36 * ccomps * dcomps);

            auto g_xxyy_yyz = cbuffer.data(gf_off + 37 * ccomps * dcomps);

            auto g_xxyy_yzz = cbuffer.data(gf_off + 38 * ccomps * dcomps);

            auto g_xxyy_zzz = cbuffer.data(gf_off + 39 * ccomps * dcomps);

            auto g_xxyz_xxx = cbuffer.data(gf_off + 40 * ccomps * dcomps);

            auto g_xxyz_xxy = cbuffer.data(gf_off + 41 * ccomps * dcomps);

            auto g_xxyz_xxz = cbuffer.data(gf_off + 42 * ccomps * dcomps);

            auto g_xxyz_xyy = cbuffer.data(gf_off + 43 * ccomps * dcomps);

            auto g_xxyz_xyz = cbuffer.data(gf_off + 44 * ccomps * dcomps);

            auto g_xxyz_xzz = cbuffer.data(gf_off + 45 * ccomps * dcomps);

            auto g_xxyz_yyy = cbuffer.data(gf_off + 46 * ccomps * dcomps);

            auto g_xxyz_yyz = cbuffer.data(gf_off + 47 * ccomps * dcomps);

            auto g_xxyz_yzz = cbuffer.data(gf_off + 48 * ccomps * dcomps);

            auto g_xxyz_zzz = cbuffer.data(gf_off + 49 * ccomps * dcomps);

            auto g_xxzz_xxx = cbuffer.data(gf_off + 50 * ccomps * dcomps);

            auto g_xxzz_xxy = cbuffer.data(gf_off + 51 * ccomps * dcomps);

            auto g_xxzz_xxz = cbuffer.data(gf_off + 52 * ccomps * dcomps);

            auto g_xxzz_xyy = cbuffer.data(gf_off + 53 * ccomps * dcomps);

            auto g_xxzz_xyz = cbuffer.data(gf_off + 54 * ccomps * dcomps);

            auto g_xxzz_xzz = cbuffer.data(gf_off + 55 * ccomps * dcomps);

            auto g_xxzz_yyy = cbuffer.data(gf_off + 56 * ccomps * dcomps);

            auto g_xxzz_yyz = cbuffer.data(gf_off + 57 * ccomps * dcomps);

            auto g_xxzz_yzz = cbuffer.data(gf_off + 58 * ccomps * dcomps);

            auto g_xxzz_zzz = cbuffer.data(gf_off + 59 * ccomps * dcomps);

            auto g_xyyy_xxx = cbuffer.data(gf_off + 60 * ccomps * dcomps);

            auto g_xyyy_xxy = cbuffer.data(gf_off + 61 * ccomps * dcomps);

            auto g_xyyy_xxz = cbuffer.data(gf_off + 62 * ccomps * dcomps);

            auto g_xyyy_xyy = cbuffer.data(gf_off + 63 * ccomps * dcomps);

            auto g_xyyy_xyz = cbuffer.data(gf_off + 64 * ccomps * dcomps);

            auto g_xyyy_xzz = cbuffer.data(gf_off + 65 * ccomps * dcomps);

            auto g_xyyy_yyy = cbuffer.data(gf_off + 66 * ccomps * dcomps);

            auto g_xyyy_yyz = cbuffer.data(gf_off + 67 * ccomps * dcomps);

            auto g_xyyy_yzz = cbuffer.data(gf_off + 68 * ccomps * dcomps);

            auto g_xyyy_zzz = cbuffer.data(gf_off + 69 * ccomps * dcomps);

            auto g_xyyz_xxx = cbuffer.data(gf_off + 70 * ccomps * dcomps);

            auto g_xyyz_xxy = cbuffer.data(gf_off + 71 * ccomps * dcomps);

            auto g_xyyz_xxz = cbuffer.data(gf_off + 72 * ccomps * dcomps);

            auto g_xyyz_xyy = cbuffer.data(gf_off + 73 * ccomps * dcomps);

            auto g_xyyz_xyz = cbuffer.data(gf_off + 74 * ccomps * dcomps);

            auto g_xyyz_xzz = cbuffer.data(gf_off + 75 * ccomps * dcomps);

            auto g_xyyz_yyy = cbuffer.data(gf_off + 76 * ccomps * dcomps);

            auto g_xyyz_yyz = cbuffer.data(gf_off + 77 * ccomps * dcomps);

            auto g_xyyz_yzz = cbuffer.data(gf_off + 78 * ccomps * dcomps);

            auto g_xyyz_zzz = cbuffer.data(gf_off + 79 * ccomps * dcomps);

            auto g_xyzz_xxx = cbuffer.data(gf_off + 80 * ccomps * dcomps);

            auto g_xyzz_xxy = cbuffer.data(gf_off + 81 * ccomps * dcomps);

            auto g_xyzz_xxz = cbuffer.data(gf_off + 82 * ccomps * dcomps);

            auto g_xyzz_xyy = cbuffer.data(gf_off + 83 * ccomps * dcomps);

            auto g_xyzz_xyz = cbuffer.data(gf_off + 84 * ccomps * dcomps);

            auto g_xyzz_xzz = cbuffer.data(gf_off + 85 * ccomps * dcomps);

            auto g_xyzz_yyy = cbuffer.data(gf_off + 86 * ccomps * dcomps);

            auto g_xyzz_yyz = cbuffer.data(gf_off + 87 * ccomps * dcomps);

            auto g_xyzz_yzz = cbuffer.data(gf_off + 88 * ccomps * dcomps);

            auto g_xyzz_zzz = cbuffer.data(gf_off + 89 * ccomps * dcomps);

            auto g_xzzz_xxx = cbuffer.data(gf_off + 90 * ccomps * dcomps);

            auto g_xzzz_xxy = cbuffer.data(gf_off + 91 * ccomps * dcomps);

            auto g_xzzz_xxz = cbuffer.data(gf_off + 92 * ccomps * dcomps);

            auto g_xzzz_xyy = cbuffer.data(gf_off + 93 * ccomps * dcomps);

            auto g_xzzz_xyz = cbuffer.data(gf_off + 94 * ccomps * dcomps);

            auto g_xzzz_xzz = cbuffer.data(gf_off + 95 * ccomps * dcomps);

            auto g_xzzz_yyy = cbuffer.data(gf_off + 96 * ccomps * dcomps);

            auto g_xzzz_yyz = cbuffer.data(gf_off + 97 * ccomps * dcomps);

            auto g_xzzz_yzz = cbuffer.data(gf_off + 98 * ccomps * dcomps);

            auto g_xzzz_zzz = cbuffer.data(gf_off + 99 * ccomps * dcomps);

            auto g_yyyy_xxx = cbuffer.data(gf_off + 100 * ccomps * dcomps);

            auto g_yyyy_xxy = cbuffer.data(gf_off + 101 * ccomps * dcomps);

            auto g_yyyy_xxz = cbuffer.data(gf_off + 102 * ccomps * dcomps);

            auto g_yyyy_xyy = cbuffer.data(gf_off + 103 * ccomps * dcomps);

            auto g_yyyy_xyz = cbuffer.data(gf_off + 104 * ccomps * dcomps);

            auto g_yyyy_xzz = cbuffer.data(gf_off + 105 * ccomps * dcomps);

            auto g_yyyy_yyy = cbuffer.data(gf_off + 106 * ccomps * dcomps);

            auto g_yyyy_yyz = cbuffer.data(gf_off + 107 * ccomps * dcomps);

            auto g_yyyy_yzz = cbuffer.data(gf_off + 108 * ccomps * dcomps);

            auto g_yyyy_zzz = cbuffer.data(gf_off + 109 * ccomps * dcomps);

            auto g_yyyz_xxx = cbuffer.data(gf_off + 110 * ccomps * dcomps);

            auto g_yyyz_xxy = cbuffer.data(gf_off + 111 * ccomps * dcomps);

            auto g_yyyz_xxz = cbuffer.data(gf_off + 112 * ccomps * dcomps);

            auto g_yyyz_xyy = cbuffer.data(gf_off + 113 * ccomps * dcomps);

            auto g_yyyz_xyz = cbuffer.data(gf_off + 114 * ccomps * dcomps);

            auto g_yyyz_xzz = cbuffer.data(gf_off + 115 * ccomps * dcomps);

            auto g_yyyz_yyy = cbuffer.data(gf_off + 116 * ccomps * dcomps);

            auto g_yyyz_yyz = cbuffer.data(gf_off + 117 * ccomps * dcomps);

            auto g_yyyz_yzz = cbuffer.data(gf_off + 118 * ccomps * dcomps);

            auto g_yyyz_zzz = cbuffer.data(gf_off + 119 * ccomps * dcomps);

            auto g_yyzz_xxx = cbuffer.data(gf_off + 120 * ccomps * dcomps);

            auto g_yyzz_xxy = cbuffer.data(gf_off + 121 * ccomps * dcomps);

            auto g_yyzz_xxz = cbuffer.data(gf_off + 122 * ccomps * dcomps);

            auto g_yyzz_xyy = cbuffer.data(gf_off + 123 * ccomps * dcomps);

            auto g_yyzz_xyz = cbuffer.data(gf_off + 124 * ccomps * dcomps);

            auto g_yyzz_xzz = cbuffer.data(gf_off + 125 * ccomps * dcomps);

            auto g_yyzz_yyy = cbuffer.data(gf_off + 126 * ccomps * dcomps);

            auto g_yyzz_yyz = cbuffer.data(gf_off + 127 * ccomps * dcomps);

            auto g_yyzz_yzz = cbuffer.data(gf_off + 128 * ccomps * dcomps);

            auto g_yyzz_zzz = cbuffer.data(gf_off + 129 * ccomps * dcomps);

            auto g_yzzz_xxx = cbuffer.data(gf_off + 130 * ccomps * dcomps);

            auto g_yzzz_xxy = cbuffer.data(gf_off + 131 * ccomps * dcomps);

            auto g_yzzz_xxz = cbuffer.data(gf_off + 132 * ccomps * dcomps);

            auto g_yzzz_xyy = cbuffer.data(gf_off + 133 * ccomps * dcomps);

            auto g_yzzz_xyz = cbuffer.data(gf_off + 134 * ccomps * dcomps);

            auto g_yzzz_xzz = cbuffer.data(gf_off + 135 * ccomps * dcomps);

            auto g_yzzz_yyy = cbuffer.data(gf_off + 136 * ccomps * dcomps);

            auto g_yzzz_yyz = cbuffer.data(gf_off + 137 * ccomps * dcomps);

            auto g_yzzz_yzz = cbuffer.data(gf_off + 138 * ccomps * dcomps);

            auto g_yzzz_zzz = cbuffer.data(gf_off + 139 * ccomps * dcomps);

            auto g_zzzz_xxx = cbuffer.data(gf_off + 140 * ccomps * dcomps);

            auto g_zzzz_xxy = cbuffer.data(gf_off + 141 * ccomps * dcomps);

            auto g_zzzz_xxz = cbuffer.data(gf_off + 142 * ccomps * dcomps);

            auto g_zzzz_xyy = cbuffer.data(gf_off + 143 * ccomps * dcomps);

            auto g_zzzz_xyz = cbuffer.data(gf_off + 144 * ccomps * dcomps);

            auto g_zzzz_xzz = cbuffer.data(gf_off + 145 * ccomps * dcomps);

            auto g_zzzz_yyy = cbuffer.data(gf_off + 146 * ccomps * dcomps);

            auto g_zzzz_yyz = cbuffer.data(gf_off + 147 * ccomps * dcomps);

            auto g_zzzz_yzz = cbuffer.data(gf_off + 148 * ccomps * dcomps);

            auto g_zzzz_zzz = cbuffer.data(gf_off + 149 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GFSS

            const auto gf_geom_10_off = idx_geom_10_gfxx + i * dcomps + j;

            auto g_x_0_xxxx_xxx = cbuffer.data(gf_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_xxy = cbuffer.data(gf_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_xxz = cbuffer.data(gf_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxx_xyy = cbuffer.data(gf_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxx_xyz = cbuffer.data(gf_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxx_xzz = cbuffer.data(gf_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxx_yyy = cbuffer.data(gf_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxx_yyz = cbuffer.data(gf_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxx_yzz = cbuffer.data(gf_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxx_zzz = cbuffer.data(gf_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxy_xxx = cbuffer.data(gf_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxy_xxy = cbuffer.data(gf_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxy_xxz = cbuffer.data(gf_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxy_xyy = cbuffer.data(gf_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxy_xyz = cbuffer.data(gf_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxy_xzz = cbuffer.data(gf_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxy_yyy = cbuffer.data(gf_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxy_yyz = cbuffer.data(gf_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxy_yzz = cbuffer.data(gf_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxy_zzz = cbuffer.data(gf_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxz_xxx = cbuffer.data(gf_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxz_xxy = cbuffer.data(gf_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxz_xxz = cbuffer.data(gf_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxz_xyy = cbuffer.data(gf_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxz_xyz = cbuffer.data(gf_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxz_xzz = cbuffer.data(gf_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxz_yyy = cbuffer.data(gf_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxz_yyz = cbuffer.data(gf_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxz_yzz = cbuffer.data(gf_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxz_zzz = cbuffer.data(gf_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxyy_xxx = cbuffer.data(gf_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxyy_xxy = cbuffer.data(gf_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxyy_xxz = cbuffer.data(gf_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxyy_xyy = cbuffer.data(gf_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxyy_xyz = cbuffer.data(gf_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxyy_xzz = cbuffer.data(gf_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxyy_yyy = cbuffer.data(gf_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxyy_yyz = cbuffer.data(gf_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxyy_yzz = cbuffer.data(gf_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxyy_zzz = cbuffer.data(gf_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxyz_xxx = cbuffer.data(gf_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxyz_xxy = cbuffer.data(gf_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxyz_xxz = cbuffer.data(gf_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxyz_xyy = cbuffer.data(gf_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxyz_xyz = cbuffer.data(gf_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxyz_xzz = cbuffer.data(gf_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxyz_yyy = cbuffer.data(gf_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxyz_yyz = cbuffer.data(gf_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxyz_yzz = cbuffer.data(gf_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxyz_zzz = cbuffer.data(gf_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxzz_xxx = cbuffer.data(gf_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxzz_xxy = cbuffer.data(gf_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxzz_xxz = cbuffer.data(gf_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxzz_xyy = cbuffer.data(gf_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxzz_xyz = cbuffer.data(gf_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxzz_xzz = cbuffer.data(gf_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxzz_yyy = cbuffer.data(gf_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxzz_yyz = cbuffer.data(gf_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxzz_yzz = cbuffer.data(gf_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxzz_zzz = cbuffer.data(gf_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xyyy_xxx = cbuffer.data(gf_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xyyy_xxy = cbuffer.data(gf_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xyyy_xxz = cbuffer.data(gf_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xyyy_xyy = cbuffer.data(gf_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xyyy_xyz = cbuffer.data(gf_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xyyy_xzz = cbuffer.data(gf_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xyyy_yyy = cbuffer.data(gf_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xyyy_yyz = cbuffer.data(gf_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xyyy_yzz = cbuffer.data(gf_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xyyy_zzz = cbuffer.data(gf_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xyyz_xxx = cbuffer.data(gf_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xyyz_xxy = cbuffer.data(gf_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xyyz_xxz = cbuffer.data(gf_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xyyz_xyy = cbuffer.data(gf_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xyyz_xyz = cbuffer.data(gf_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xyyz_xzz = cbuffer.data(gf_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xyyz_yyy = cbuffer.data(gf_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xyyz_yyz = cbuffer.data(gf_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xyyz_yzz = cbuffer.data(gf_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xyyz_zzz = cbuffer.data(gf_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xyzz_xxx = cbuffer.data(gf_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xyzz_xxy = cbuffer.data(gf_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xyzz_xxz = cbuffer.data(gf_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xyzz_xyy = cbuffer.data(gf_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xyzz_xyz = cbuffer.data(gf_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xyzz_xzz = cbuffer.data(gf_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xyzz_yyy = cbuffer.data(gf_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xyzz_yyz = cbuffer.data(gf_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xyzz_yzz = cbuffer.data(gf_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xyzz_zzz = cbuffer.data(gf_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xzzz_xxx = cbuffer.data(gf_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xzzz_xxy = cbuffer.data(gf_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xzzz_xxz = cbuffer.data(gf_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xzzz_xyy = cbuffer.data(gf_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xzzz_xyz = cbuffer.data(gf_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xzzz_xzz = cbuffer.data(gf_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xzzz_yyy = cbuffer.data(gf_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xzzz_yyz = cbuffer.data(gf_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xzzz_yzz = cbuffer.data(gf_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xzzz_zzz = cbuffer.data(gf_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_yyyy_xxx = cbuffer.data(gf_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_yyyy_xxy = cbuffer.data(gf_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_yyyy_xxz = cbuffer.data(gf_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_yyyy_xyy = cbuffer.data(gf_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_yyyy_xyz = cbuffer.data(gf_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_yyyy_xzz = cbuffer.data(gf_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_yyyy_yyy = cbuffer.data(gf_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_yyyy_yyz = cbuffer.data(gf_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_yyyy_yzz = cbuffer.data(gf_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_yyyy_zzz = cbuffer.data(gf_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_yyyz_xxx = cbuffer.data(gf_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_yyyz_xxy = cbuffer.data(gf_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_yyyz_xxz = cbuffer.data(gf_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_yyyz_xyy = cbuffer.data(gf_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_yyyz_xyz = cbuffer.data(gf_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_yyyz_xzz = cbuffer.data(gf_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_yyyz_yyy = cbuffer.data(gf_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_yyyz_yyz = cbuffer.data(gf_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_yyyz_yzz = cbuffer.data(gf_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_yyyz_zzz = cbuffer.data(gf_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_yyzz_xxx = cbuffer.data(gf_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_yyzz_xxy = cbuffer.data(gf_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_yyzz_xxz = cbuffer.data(gf_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_yyzz_xyy = cbuffer.data(gf_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_yyzz_xyz = cbuffer.data(gf_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_yyzz_xzz = cbuffer.data(gf_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_yyzz_yyy = cbuffer.data(gf_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_yyzz_yyz = cbuffer.data(gf_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_yyzz_yzz = cbuffer.data(gf_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_yyzz_zzz = cbuffer.data(gf_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_yzzz_xxx = cbuffer.data(gf_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_yzzz_xxy = cbuffer.data(gf_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_yzzz_xxz = cbuffer.data(gf_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_yzzz_xyy = cbuffer.data(gf_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_yzzz_xyz = cbuffer.data(gf_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_yzzz_xzz = cbuffer.data(gf_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_yzzz_yyy = cbuffer.data(gf_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_yzzz_yyz = cbuffer.data(gf_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_yzzz_yzz = cbuffer.data(gf_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_yzzz_zzz = cbuffer.data(gf_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_zzzz_xxx = cbuffer.data(gf_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_zzzz_xxy = cbuffer.data(gf_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_zzzz_xxz = cbuffer.data(gf_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_zzzz_xyy = cbuffer.data(gf_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_zzzz_xyz = cbuffer.data(gf_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_zzzz_xzz = cbuffer.data(gf_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_zzzz_yyy = cbuffer.data(gf_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_zzzz_yyz = cbuffer.data(gf_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_zzzz_yzz = cbuffer.data(gf_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_zzzz_zzz = cbuffer.data(gf_geom_10_off + 149 * ccomps * dcomps);

            auto g_y_0_xxxx_xxx = cbuffer.data(gf_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_xxxx_xxy = cbuffer.data(gf_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_xxxx_xxz = cbuffer.data(gf_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_xxxx_xyy = cbuffer.data(gf_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_xxxx_xyz = cbuffer.data(gf_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_xxxx_xzz = cbuffer.data(gf_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_xxxx_yyy = cbuffer.data(gf_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_xxxx_yyz = cbuffer.data(gf_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_xxxx_yzz = cbuffer.data(gf_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_xxxx_zzz = cbuffer.data(gf_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_xxxy_xxx = cbuffer.data(gf_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_xxxy_xxy = cbuffer.data(gf_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_xxxy_xxz = cbuffer.data(gf_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_xxxy_xyy = cbuffer.data(gf_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_xxxy_xyz = cbuffer.data(gf_geom_10_off + 164 * ccomps * dcomps);

            auto g_y_0_xxxy_xzz = cbuffer.data(gf_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_xxxy_yyy = cbuffer.data(gf_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_xxxy_yyz = cbuffer.data(gf_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_xxxy_yzz = cbuffer.data(gf_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_xxxy_zzz = cbuffer.data(gf_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_xxxz_xxx = cbuffer.data(gf_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_xxxz_xxy = cbuffer.data(gf_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_xxxz_xxz = cbuffer.data(gf_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_xxxz_xyy = cbuffer.data(gf_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_xxxz_xyz = cbuffer.data(gf_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_xxxz_xzz = cbuffer.data(gf_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_xxxz_yyy = cbuffer.data(gf_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_xxxz_yyz = cbuffer.data(gf_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_xxxz_yzz = cbuffer.data(gf_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_xxxz_zzz = cbuffer.data(gf_geom_10_off + 179 * ccomps * dcomps);

            auto g_y_0_xxyy_xxx = cbuffer.data(gf_geom_10_off + 180 * ccomps * dcomps);

            auto g_y_0_xxyy_xxy = cbuffer.data(gf_geom_10_off + 181 * ccomps * dcomps);

            auto g_y_0_xxyy_xxz = cbuffer.data(gf_geom_10_off + 182 * ccomps * dcomps);

            auto g_y_0_xxyy_xyy = cbuffer.data(gf_geom_10_off + 183 * ccomps * dcomps);

            auto g_y_0_xxyy_xyz = cbuffer.data(gf_geom_10_off + 184 * ccomps * dcomps);

            auto g_y_0_xxyy_xzz = cbuffer.data(gf_geom_10_off + 185 * ccomps * dcomps);

            auto g_y_0_xxyy_yyy = cbuffer.data(gf_geom_10_off + 186 * ccomps * dcomps);

            auto g_y_0_xxyy_yyz = cbuffer.data(gf_geom_10_off + 187 * ccomps * dcomps);

            auto g_y_0_xxyy_yzz = cbuffer.data(gf_geom_10_off + 188 * ccomps * dcomps);

            auto g_y_0_xxyy_zzz = cbuffer.data(gf_geom_10_off + 189 * ccomps * dcomps);

            auto g_y_0_xxyz_xxx = cbuffer.data(gf_geom_10_off + 190 * ccomps * dcomps);

            auto g_y_0_xxyz_xxy = cbuffer.data(gf_geom_10_off + 191 * ccomps * dcomps);

            auto g_y_0_xxyz_xxz = cbuffer.data(gf_geom_10_off + 192 * ccomps * dcomps);

            auto g_y_0_xxyz_xyy = cbuffer.data(gf_geom_10_off + 193 * ccomps * dcomps);

            auto g_y_0_xxyz_xyz = cbuffer.data(gf_geom_10_off + 194 * ccomps * dcomps);

            auto g_y_0_xxyz_xzz = cbuffer.data(gf_geom_10_off + 195 * ccomps * dcomps);

            auto g_y_0_xxyz_yyy = cbuffer.data(gf_geom_10_off + 196 * ccomps * dcomps);

            auto g_y_0_xxyz_yyz = cbuffer.data(gf_geom_10_off + 197 * ccomps * dcomps);

            auto g_y_0_xxyz_yzz = cbuffer.data(gf_geom_10_off + 198 * ccomps * dcomps);

            auto g_y_0_xxyz_zzz = cbuffer.data(gf_geom_10_off + 199 * ccomps * dcomps);

            auto g_y_0_xxzz_xxx = cbuffer.data(gf_geom_10_off + 200 * ccomps * dcomps);

            auto g_y_0_xxzz_xxy = cbuffer.data(gf_geom_10_off + 201 * ccomps * dcomps);

            auto g_y_0_xxzz_xxz = cbuffer.data(gf_geom_10_off + 202 * ccomps * dcomps);

            auto g_y_0_xxzz_xyy = cbuffer.data(gf_geom_10_off + 203 * ccomps * dcomps);

            auto g_y_0_xxzz_xyz = cbuffer.data(gf_geom_10_off + 204 * ccomps * dcomps);

            auto g_y_0_xxzz_xzz = cbuffer.data(gf_geom_10_off + 205 * ccomps * dcomps);

            auto g_y_0_xxzz_yyy = cbuffer.data(gf_geom_10_off + 206 * ccomps * dcomps);

            auto g_y_0_xxzz_yyz = cbuffer.data(gf_geom_10_off + 207 * ccomps * dcomps);

            auto g_y_0_xxzz_yzz = cbuffer.data(gf_geom_10_off + 208 * ccomps * dcomps);

            auto g_y_0_xxzz_zzz = cbuffer.data(gf_geom_10_off + 209 * ccomps * dcomps);

            auto g_y_0_xyyy_xxx = cbuffer.data(gf_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_xyyy_xxy = cbuffer.data(gf_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_xyyy_xxz = cbuffer.data(gf_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_xyyy_xyy = cbuffer.data(gf_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_xyyy_xyz = cbuffer.data(gf_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_xyyy_xzz = cbuffer.data(gf_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_xyyy_yyy = cbuffer.data(gf_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_xyyy_yyz = cbuffer.data(gf_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_xyyy_yzz = cbuffer.data(gf_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_xyyy_zzz = cbuffer.data(gf_geom_10_off + 219 * ccomps * dcomps);

            auto g_y_0_xyyz_xxx = cbuffer.data(gf_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_xyyz_xxy = cbuffer.data(gf_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_xyyz_xxz = cbuffer.data(gf_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_xyyz_xyy = cbuffer.data(gf_geom_10_off + 223 * ccomps * dcomps);

            auto g_y_0_xyyz_xyz = cbuffer.data(gf_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xyyz_xzz = cbuffer.data(gf_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xyyz_yyy = cbuffer.data(gf_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xyyz_yyz = cbuffer.data(gf_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xyyz_yzz = cbuffer.data(gf_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xyyz_zzz = cbuffer.data(gf_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xyzz_xxx = cbuffer.data(gf_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xyzz_xxy = cbuffer.data(gf_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xyzz_xxz = cbuffer.data(gf_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xyzz_xyy = cbuffer.data(gf_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xyzz_xyz = cbuffer.data(gf_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xyzz_xzz = cbuffer.data(gf_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xyzz_yyy = cbuffer.data(gf_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xyzz_yyz = cbuffer.data(gf_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xyzz_yzz = cbuffer.data(gf_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xyzz_zzz = cbuffer.data(gf_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_xzzz_xxx = cbuffer.data(gf_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xzzz_xxy = cbuffer.data(gf_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xzzz_xxz = cbuffer.data(gf_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xzzz_xyy = cbuffer.data(gf_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xzzz_xyz = cbuffer.data(gf_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_xzzz_xzz = cbuffer.data(gf_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_xzzz_yyy = cbuffer.data(gf_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_xzzz_yyz = cbuffer.data(gf_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_xzzz_yzz = cbuffer.data(gf_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_xzzz_zzz = cbuffer.data(gf_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_yyyy_xxx = cbuffer.data(gf_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_yyyy_xxy = cbuffer.data(gf_geom_10_off + 251 * ccomps * dcomps);

            auto g_y_0_yyyy_xxz = cbuffer.data(gf_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_yyyy_xyy = cbuffer.data(gf_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_yyyy_xyz = cbuffer.data(gf_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_yyyy_xzz = cbuffer.data(gf_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_yyyy_yyy = cbuffer.data(gf_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_yyyy_yyz = cbuffer.data(gf_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_yyyy_yzz = cbuffer.data(gf_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_yyyy_zzz = cbuffer.data(gf_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_yyyz_xxx = cbuffer.data(gf_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_yyyz_xxy = cbuffer.data(gf_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_yyyz_xxz = cbuffer.data(gf_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_yyyz_xyy = cbuffer.data(gf_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_yyyz_xyz = cbuffer.data(gf_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_yyyz_xzz = cbuffer.data(gf_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_yyyz_yyy = cbuffer.data(gf_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_yyyz_yyz = cbuffer.data(gf_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_yyyz_yzz = cbuffer.data(gf_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_yyyz_zzz = cbuffer.data(gf_geom_10_off + 269 * ccomps * dcomps);

            auto g_y_0_yyzz_xxx = cbuffer.data(gf_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_yyzz_xxy = cbuffer.data(gf_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_yyzz_xxz = cbuffer.data(gf_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_yyzz_xyy = cbuffer.data(gf_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_yyzz_xyz = cbuffer.data(gf_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_yyzz_xzz = cbuffer.data(gf_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_yyzz_yyy = cbuffer.data(gf_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_yyzz_yyz = cbuffer.data(gf_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_yyzz_yzz = cbuffer.data(gf_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_yyzz_zzz = cbuffer.data(gf_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_yzzz_xxx = cbuffer.data(gf_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_yzzz_xxy = cbuffer.data(gf_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_yzzz_xxz = cbuffer.data(gf_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_yzzz_xyy = cbuffer.data(gf_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_yzzz_xyz = cbuffer.data(gf_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_yzzz_xzz = cbuffer.data(gf_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_yzzz_yyy = cbuffer.data(gf_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_yzzz_yyz = cbuffer.data(gf_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_yzzz_yzz = cbuffer.data(gf_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_yzzz_zzz = cbuffer.data(gf_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_zzzz_xxx = cbuffer.data(gf_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_zzzz_xxy = cbuffer.data(gf_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_zzzz_xxz = cbuffer.data(gf_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_zzzz_xyy = cbuffer.data(gf_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_zzzz_xyz = cbuffer.data(gf_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_zzzz_xzz = cbuffer.data(gf_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_zzzz_yyy = cbuffer.data(gf_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_zzzz_yyz = cbuffer.data(gf_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_zzzz_yzz = cbuffer.data(gf_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_zzzz_zzz = cbuffer.data(gf_geom_10_off + 299 * ccomps * dcomps);

            auto g_z_0_xxxx_xxx = cbuffer.data(gf_geom_10_off + 300 * ccomps * dcomps);

            auto g_z_0_xxxx_xxy = cbuffer.data(gf_geom_10_off + 301 * ccomps * dcomps);

            auto g_z_0_xxxx_xxz = cbuffer.data(gf_geom_10_off + 302 * ccomps * dcomps);

            auto g_z_0_xxxx_xyy = cbuffer.data(gf_geom_10_off + 303 * ccomps * dcomps);

            auto g_z_0_xxxx_xyz = cbuffer.data(gf_geom_10_off + 304 * ccomps * dcomps);

            auto g_z_0_xxxx_xzz = cbuffer.data(gf_geom_10_off + 305 * ccomps * dcomps);

            auto g_z_0_xxxx_yyy = cbuffer.data(gf_geom_10_off + 306 * ccomps * dcomps);

            auto g_z_0_xxxx_yyz = cbuffer.data(gf_geom_10_off + 307 * ccomps * dcomps);

            auto g_z_0_xxxx_yzz = cbuffer.data(gf_geom_10_off + 308 * ccomps * dcomps);

            auto g_z_0_xxxx_zzz = cbuffer.data(gf_geom_10_off + 309 * ccomps * dcomps);

            auto g_z_0_xxxy_xxx = cbuffer.data(gf_geom_10_off + 310 * ccomps * dcomps);

            auto g_z_0_xxxy_xxy = cbuffer.data(gf_geom_10_off + 311 * ccomps * dcomps);

            auto g_z_0_xxxy_xxz = cbuffer.data(gf_geom_10_off + 312 * ccomps * dcomps);

            auto g_z_0_xxxy_xyy = cbuffer.data(gf_geom_10_off + 313 * ccomps * dcomps);

            auto g_z_0_xxxy_xyz = cbuffer.data(gf_geom_10_off + 314 * ccomps * dcomps);

            auto g_z_0_xxxy_xzz = cbuffer.data(gf_geom_10_off + 315 * ccomps * dcomps);

            auto g_z_0_xxxy_yyy = cbuffer.data(gf_geom_10_off + 316 * ccomps * dcomps);

            auto g_z_0_xxxy_yyz = cbuffer.data(gf_geom_10_off + 317 * ccomps * dcomps);

            auto g_z_0_xxxy_yzz = cbuffer.data(gf_geom_10_off + 318 * ccomps * dcomps);

            auto g_z_0_xxxy_zzz = cbuffer.data(gf_geom_10_off + 319 * ccomps * dcomps);

            auto g_z_0_xxxz_xxx = cbuffer.data(gf_geom_10_off + 320 * ccomps * dcomps);

            auto g_z_0_xxxz_xxy = cbuffer.data(gf_geom_10_off + 321 * ccomps * dcomps);

            auto g_z_0_xxxz_xxz = cbuffer.data(gf_geom_10_off + 322 * ccomps * dcomps);

            auto g_z_0_xxxz_xyy = cbuffer.data(gf_geom_10_off + 323 * ccomps * dcomps);

            auto g_z_0_xxxz_xyz = cbuffer.data(gf_geom_10_off + 324 * ccomps * dcomps);

            auto g_z_0_xxxz_xzz = cbuffer.data(gf_geom_10_off + 325 * ccomps * dcomps);

            auto g_z_0_xxxz_yyy = cbuffer.data(gf_geom_10_off + 326 * ccomps * dcomps);

            auto g_z_0_xxxz_yyz = cbuffer.data(gf_geom_10_off + 327 * ccomps * dcomps);

            auto g_z_0_xxxz_yzz = cbuffer.data(gf_geom_10_off + 328 * ccomps * dcomps);

            auto g_z_0_xxxz_zzz = cbuffer.data(gf_geom_10_off + 329 * ccomps * dcomps);

            auto g_z_0_xxyy_xxx = cbuffer.data(gf_geom_10_off + 330 * ccomps * dcomps);

            auto g_z_0_xxyy_xxy = cbuffer.data(gf_geom_10_off + 331 * ccomps * dcomps);

            auto g_z_0_xxyy_xxz = cbuffer.data(gf_geom_10_off + 332 * ccomps * dcomps);

            auto g_z_0_xxyy_xyy = cbuffer.data(gf_geom_10_off + 333 * ccomps * dcomps);

            auto g_z_0_xxyy_xyz = cbuffer.data(gf_geom_10_off + 334 * ccomps * dcomps);

            auto g_z_0_xxyy_xzz = cbuffer.data(gf_geom_10_off + 335 * ccomps * dcomps);

            auto g_z_0_xxyy_yyy = cbuffer.data(gf_geom_10_off + 336 * ccomps * dcomps);

            auto g_z_0_xxyy_yyz = cbuffer.data(gf_geom_10_off + 337 * ccomps * dcomps);

            auto g_z_0_xxyy_yzz = cbuffer.data(gf_geom_10_off + 338 * ccomps * dcomps);

            auto g_z_0_xxyy_zzz = cbuffer.data(gf_geom_10_off + 339 * ccomps * dcomps);

            auto g_z_0_xxyz_xxx = cbuffer.data(gf_geom_10_off + 340 * ccomps * dcomps);

            auto g_z_0_xxyz_xxy = cbuffer.data(gf_geom_10_off + 341 * ccomps * dcomps);

            auto g_z_0_xxyz_xxz = cbuffer.data(gf_geom_10_off + 342 * ccomps * dcomps);

            auto g_z_0_xxyz_xyy = cbuffer.data(gf_geom_10_off + 343 * ccomps * dcomps);

            auto g_z_0_xxyz_xyz = cbuffer.data(gf_geom_10_off + 344 * ccomps * dcomps);

            auto g_z_0_xxyz_xzz = cbuffer.data(gf_geom_10_off + 345 * ccomps * dcomps);

            auto g_z_0_xxyz_yyy = cbuffer.data(gf_geom_10_off + 346 * ccomps * dcomps);

            auto g_z_0_xxyz_yyz = cbuffer.data(gf_geom_10_off + 347 * ccomps * dcomps);

            auto g_z_0_xxyz_yzz = cbuffer.data(gf_geom_10_off + 348 * ccomps * dcomps);

            auto g_z_0_xxyz_zzz = cbuffer.data(gf_geom_10_off + 349 * ccomps * dcomps);

            auto g_z_0_xxzz_xxx = cbuffer.data(gf_geom_10_off + 350 * ccomps * dcomps);

            auto g_z_0_xxzz_xxy = cbuffer.data(gf_geom_10_off + 351 * ccomps * dcomps);

            auto g_z_0_xxzz_xxz = cbuffer.data(gf_geom_10_off + 352 * ccomps * dcomps);

            auto g_z_0_xxzz_xyy = cbuffer.data(gf_geom_10_off + 353 * ccomps * dcomps);

            auto g_z_0_xxzz_xyz = cbuffer.data(gf_geom_10_off + 354 * ccomps * dcomps);

            auto g_z_0_xxzz_xzz = cbuffer.data(gf_geom_10_off + 355 * ccomps * dcomps);

            auto g_z_0_xxzz_yyy = cbuffer.data(gf_geom_10_off + 356 * ccomps * dcomps);

            auto g_z_0_xxzz_yyz = cbuffer.data(gf_geom_10_off + 357 * ccomps * dcomps);

            auto g_z_0_xxzz_yzz = cbuffer.data(gf_geom_10_off + 358 * ccomps * dcomps);

            auto g_z_0_xxzz_zzz = cbuffer.data(gf_geom_10_off + 359 * ccomps * dcomps);

            auto g_z_0_xyyy_xxx = cbuffer.data(gf_geom_10_off + 360 * ccomps * dcomps);

            auto g_z_0_xyyy_xxy = cbuffer.data(gf_geom_10_off + 361 * ccomps * dcomps);

            auto g_z_0_xyyy_xxz = cbuffer.data(gf_geom_10_off + 362 * ccomps * dcomps);

            auto g_z_0_xyyy_xyy = cbuffer.data(gf_geom_10_off + 363 * ccomps * dcomps);

            auto g_z_0_xyyy_xyz = cbuffer.data(gf_geom_10_off + 364 * ccomps * dcomps);

            auto g_z_0_xyyy_xzz = cbuffer.data(gf_geom_10_off + 365 * ccomps * dcomps);

            auto g_z_0_xyyy_yyy = cbuffer.data(gf_geom_10_off + 366 * ccomps * dcomps);

            auto g_z_0_xyyy_yyz = cbuffer.data(gf_geom_10_off + 367 * ccomps * dcomps);

            auto g_z_0_xyyy_yzz = cbuffer.data(gf_geom_10_off + 368 * ccomps * dcomps);

            auto g_z_0_xyyy_zzz = cbuffer.data(gf_geom_10_off + 369 * ccomps * dcomps);

            auto g_z_0_xyyz_xxx = cbuffer.data(gf_geom_10_off + 370 * ccomps * dcomps);

            auto g_z_0_xyyz_xxy = cbuffer.data(gf_geom_10_off + 371 * ccomps * dcomps);

            auto g_z_0_xyyz_xxz = cbuffer.data(gf_geom_10_off + 372 * ccomps * dcomps);

            auto g_z_0_xyyz_xyy = cbuffer.data(gf_geom_10_off + 373 * ccomps * dcomps);

            auto g_z_0_xyyz_xyz = cbuffer.data(gf_geom_10_off + 374 * ccomps * dcomps);

            auto g_z_0_xyyz_xzz = cbuffer.data(gf_geom_10_off + 375 * ccomps * dcomps);

            auto g_z_0_xyyz_yyy = cbuffer.data(gf_geom_10_off + 376 * ccomps * dcomps);

            auto g_z_0_xyyz_yyz = cbuffer.data(gf_geom_10_off + 377 * ccomps * dcomps);

            auto g_z_0_xyyz_yzz = cbuffer.data(gf_geom_10_off + 378 * ccomps * dcomps);

            auto g_z_0_xyyz_zzz = cbuffer.data(gf_geom_10_off + 379 * ccomps * dcomps);

            auto g_z_0_xyzz_xxx = cbuffer.data(gf_geom_10_off + 380 * ccomps * dcomps);

            auto g_z_0_xyzz_xxy = cbuffer.data(gf_geom_10_off + 381 * ccomps * dcomps);

            auto g_z_0_xyzz_xxz = cbuffer.data(gf_geom_10_off + 382 * ccomps * dcomps);

            auto g_z_0_xyzz_xyy = cbuffer.data(gf_geom_10_off + 383 * ccomps * dcomps);

            auto g_z_0_xyzz_xyz = cbuffer.data(gf_geom_10_off + 384 * ccomps * dcomps);

            auto g_z_0_xyzz_xzz = cbuffer.data(gf_geom_10_off + 385 * ccomps * dcomps);

            auto g_z_0_xyzz_yyy = cbuffer.data(gf_geom_10_off + 386 * ccomps * dcomps);

            auto g_z_0_xyzz_yyz = cbuffer.data(gf_geom_10_off + 387 * ccomps * dcomps);

            auto g_z_0_xyzz_yzz = cbuffer.data(gf_geom_10_off + 388 * ccomps * dcomps);

            auto g_z_0_xyzz_zzz = cbuffer.data(gf_geom_10_off + 389 * ccomps * dcomps);

            auto g_z_0_xzzz_xxx = cbuffer.data(gf_geom_10_off + 390 * ccomps * dcomps);

            auto g_z_0_xzzz_xxy = cbuffer.data(gf_geom_10_off + 391 * ccomps * dcomps);

            auto g_z_0_xzzz_xxz = cbuffer.data(gf_geom_10_off + 392 * ccomps * dcomps);

            auto g_z_0_xzzz_xyy = cbuffer.data(gf_geom_10_off + 393 * ccomps * dcomps);

            auto g_z_0_xzzz_xyz = cbuffer.data(gf_geom_10_off + 394 * ccomps * dcomps);

            auto g_z_0_xzzz_xzz = cbuffer.data(gf_geom_10_off + 395 * ccomps * dcomps);

            auto g_z_0_xzzz_yyy = cbuffer.data(gf_geom_10_off + 396 * ccomps * dcomps);

            auto g_z_0_xzzz_yyz = cbuffer.data(gf_geom_10_off + 397 * ccomps * dcomps);

            auto g_z_0_xzzz_yzz = cbuffer.data(gf_geom_10_off + 398 * ccomps * dcomps);

            auto g_z_0_xzzz_zzz = cbuffer.data(gf_geom_10_off + 399 * ccomps * dcomps);

            auto g_z_0_yyyy_xxx = cbuffer.data(gf_geom_10_off + 400 * ccomps * dcomps);

            auto g_z_0_yyyy_xxy = cbuffer.data(gf_geom_10_off + 401 * ccomps * dcomps);

            auto g_z_0_yyyy_xxz = cbuffer.data(gf_geom_10_off + 402 * ccomps * dcomps);

            auto g_z_0_yyyy_xyy = cbuffer.data(gf_geom_10_off + 403 * ccomps * dcomps);

            auto g_z_0_yyyy_xyz = cbuffer.data(gf_geom_10_off + 404 * ccomps * dcomps);

            auto g_z_0_yyyy_xzz = cbuffer.data(gf_geom_10_off + 405 * ccomps * dcomps);

            auto g_z_0_yyyy_yyy = cbuffer.data(gf_geom_10_off + 406 * ccomps * dcomps);

            auto g_z_0_yyyy_yyz = cbuffer.data(gf_geom_10_off + 407 * ccomps * dcomps);

            auto g_z_0_yyyy_yzz = cbuffer.data(gf_geom_10_off + 408 * ccomps * dcomps);

            auto g_z_0_yyyy_zzz = cbuffer.data(gf_geom_10_off + 409 * ccomps * dcomps);

            auto g_z_0_yyyz_xxx = cbuffer.data(gf_geom_10_off + 410 * ccomps * dcomps);

            auto g_z_0_yyyz_xxy = cbuffer.data(gf_geom_10_off + 411 * ccomps * dcomps);

            auto g_z_0_yyyz_xxz = cbuffer.data(gf_geom_10_off + 412 * ccomps * dcomps);

            auto g_z_0_yyyz_xyy = cbuffer.data(gf_geom_10_off + 413 * ccomps * dcomps);

            auto g_z_0_yyyz_xyz = cbuffer.data(gf_geom_10_off + 414 * ccomps * dcomps);

            auto g_z_0_yyyz_xzz = cbuffer.data(gf_geom_10_off + 415 * ccomps * dcomps);

            auto g_z_0_yyyz_yyy = cbuffer.data(gf_geom_10_off + 416 * ccomps * dcomps);

            auto g_z_0_yyyz_yyz = cbuffer.data(gf_geom_10_off + 417 * ccomps * dcomps);

            auto g_z_0_yyyz_yzz = cbuffer.data(gf_geom_10_off + 418 * ccomps * dcomps);

            auto g_z_0_yyyz_zzz = cbuffer.data(gf_geom_10_off + 419 * ccomps * dcomps);

            auto g_z_0_yyzz_xxx = cbuffer.data(gf_geom_10_off + 420 * ccomps * dcomps);

            auto g_z_0_yyzz_xxy = cbuffer.data(gf_geom_10_off + 421 * ccomps * dcomps);

            auto g_z_0_yyzz_xxz = cbuffer.data(gf_geom_10_off + 422 * ccomps * dcomps);

            auto g_z_0_yyzz_xyy = cbuffer.data(gf_geom_10_off + 423 * ccomps * dcomps);

            auto g_z_0_yyzz_xyz = cbuffer.data(gf_geom_10_off + 424 * ccomps * dcomps);

            auto g_z_0_yyzz_xzz = cbuffer.data(gf_geom_10_off + 425 * ccomps * dcomps);

            auto g_z_0_yyzz_yyy = cbuffer.data(gf_geom_10_off + 426 * ccomps * dcomps);

            auto g_z_0_yyzz_yyz = cbuffer.data(gf_geom_10_off + 427 * ccomps * dcomps);

            auto g_z_0_yyzz_yzz = cbuffer.data(gf_geom_10_off + 428 * ccomps * dcomps);

            auto g_z_0_yyzz_zzz = cbuffer.data(gf_geom_10_off + 429 * ccomps * dcomps);

            auto g_z_0_yzzz_xxx = cbuffer.data(gf_geom_10_off + 430 * ccomps * dcomps);

            auto g_z_0_yzzz_xxy = cbuffer.data(gf_geom_10_off + 431 * ccomps * dcomps);

            auto g_z_0_yzzz_xxz = cbuffer.data(gf_geom_10_off + 432 * ccomps * dcomps);

            auto g_z_0_yzzz_xyy = cbuffer.data(gf_geom_10_off + 433 * ccomps * dcomps);

            auto g_z_0_yzzz_xyz = cbuffer.data(gf_geom_10_off + 434 * ccomps * dcomps);

            auto g_z_0_yzzz_xzz = cbuffer.data(gf_geom_10_off + 435 * ccomps * dcomps);

            auto g_z_0_yzzz_yyy = cbuffer.data(gf_geom_10_off + 436 * ccomps * dcomps);

            auto g_z_0_yzzz_yyz = cbuffer.data(gf_geom_10_off + 437 * ccomps * dcomps);

            auto g_z_0_yzzz_yzz = cbuffer.data(gf_geom_10_off + 438 * ccomps * dcomps);

            auto g_z_0_yzzz_zzz = cbuffer.data(gf_geom_10_off + 439 * ccomps * dcomps);

            auto g_z_0_zzzz_xxx = cbuffer.data(gf_geom_10_off + 440 * ccomps * dcomps);

            auto g_z_0_zzzz_xxy = cbuffer.data(gf_geom_10_off + 441 * ccomps * dcomps);

            auto g_z_0_zzzz_xxz = cbuffer.data(gf_geom_10_off + 442 * ccomps * dcomps);

            auto g_z_0_zzzz_xyy = cbuffer.data(gf_geom_10_off + 443 * ccomps * dcomps);

            auto g_z_0_zzzz_xyz = cbuffer.data(gf_geom_10_off + 444 * ccomps * dcomps);

            auto g_z_0_zzzz_xzz = cbuffer.data(gf_geom_10_off + 445 * ccomps * dcomps);

            auto g_z_0_zzzz_yyy = cbuffer.data(gf_geom_10_off + 446 * ccomps * dcomps);

            auto g_z_0_zzzz_yyz = cbuffer.data(gf_geom_10_off + 447 * ccomps * dcomps);

            auto g_z_0_zzzz_yzz = cbuffer.data(gf_geom_10_off + 448 * ccomps * dcomps);

            auto g_z_0_zzzz_zzz = cbuffer.data(gf_geom_10_off + 449 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GGSS

            const auto gg_geom_10_off = idx_geom_10_ggxx + i * dcomps + j;

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

            /// set up bra offset for contr_buffer_hfxx

            const auto hf_geom_10_off = idx_geom_10_hfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxx_xxx, g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxy, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xyy, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_yyy, g_x_0_xxxx_yyz, g_x_0_xxxx_yzz, g_x_0_xxxx_zzz, g_x_0_xxxxx_xxx, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_yyy, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_zzz, g_xxxx_xxx, g_xxxx_xxy, g_xxxx_xxz, g_xxxx_xyy, g_xxxx_xyz, g_xxxx_xzz, g_xxxx_yyy, g_xxxx_yyz, g_xxxx_yzz, g_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_xxx[k] = -g_xxxx_xxx[k] - g_x_0_xxxx_xxx[k] * ab_x + g_x_0_xxxx_xxxx[k];

                g_x_0_xxxxx_xxy[k] = -g_xxxx_xxy[k] - g_x_0_xxxx_xxy[k] * ab_x + g_x_0_xxxx_xxxy[k];

                g_x_0_xxxxx_xxz[k] = -g_xxxx_xxz[k] - g_x_0_xxxx_xxz[k] * ab_x + g_x_0_xxxx_xxxz[k];

                g_x_0_xxxxx_xyy[k] = -g_xxxx_xyy[k] - g_x_0_xxxx_xyy[k] * ab_x + g_x_0_xxxx_xxyy[k];

                g_x_0_xxxxx_xyz[k] = -g_xxxx_xyz[k] - g_x_0_xxxx_xyz[k] * ab_x + g_x_0_xxxx_xxyz[k];

                g_x_0_xxxxx_xzz[k] = -g_xxxx_xzz[k] - g_x_0_xxxx_xzz[k] * ab_x + g_x_0_xxxx_xxzz[k];

                g_x_0_xxxxx_yyy[k] = -g_xxxx_yyy[k] - g_x_0_xxxx_yyy[k] * ab_x + g_x_0_xxxx_xyyy[k];

                g_x_0_xxxxx_yyz[k] = -g_xxxx_yyz[k] - g_x_0_xxxx_yyz[k] * ab_x + g_x_0_xxxx_xyyz[k];

                g_x_0_xxxxx_yzz[k] = -g_xxxx_yzz[k] - g_x_0_xxxx_yzz[k] * ab_x + g_x_0_xxxx_xyzz[k];

                g_x_0_xxxxx_zzz[k] = -g_xxxx_zzz[k] - g_x_0_xxxx_zzz[k] * ab_x + g_x_0_xxxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxx_xxx, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxy, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxz, g_x_0_xxxx_xyy, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xzz, g_x_0_xxxx_yyy, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_zzz, g_x_0_xxxxy_xxx, g_x_0_xxxxy_xxy, g_x_0_xxxxy_xxz, g_x_0_xxxxy_xyy, g_x_0_xxxxy_xyz, g_x_0_xxxxy_xzz, g_x_0_xxxxy_yyy, g_x_0_xxxxy_yyz, g_x_0_xxxxy_yzz, g_x_0_xxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_xxx[k] = -g_x_0_xxxx_xxx[k] * ab_y + g_x_0_xxxx_xxxy[k];

                g_x_0_xxxxy_xxy[k] = -g_x_0_xxxx_xxy[k] * ab_y + g_x_0_xxxx_xxyy[k];

                g_x_0_xxxxy_xxz[k] = -g_x_0_xxxx_xxz[k] * ab_y + g_x_0_xxxx_xxyz[k];

                g_x_0_xxxxy_xyy[k] = -g_x_0_xxxx_xyy[k] * ab_y + g_x_0_xxxx_xyyy[k];

                g_x_0_xxxxy_xyz[k] = -g_x_0_xxxx_xyz[k] * ab_y + g_x_0_xxxx_xyyz[k];

                g_x_0_xxxxy_xzz[k] = -g_x_0_xxxx_xzz[k] * ab_y + g_x_0_xxxx_xyzz[k];

                g_x_0_xxxxy_yyy[k] = -g_x_0_xxxx_yyy[k] * ab_y + g_x_0_xxxx_yyyy[k];

                g_x_0_xxxxy_yyz[k] = -g_x_0_xxxx_yyz[k] * ab_y + g_x_0_xxxx_yyyz[k];

                g_x_0_xxxxy_yzz[k] = -g_x_0_xxxx_yzz[k] * ab_y + g_x_0_xxxx_yyzz[k];

                g_x_0_xxxxy_zzz[k] = -g_x_0_xxxx_zzz[k] * ab_y + g_x_0_xxxx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxx_xxx, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxy, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xyy, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_yyy, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_zzz, g_x_0_xxxx_zzzz, g_x_0_xxxxz_xxx, g_x_0_xxxxz_xxy, g_x_0_xxxxz_xxz, g_x_0_xxxxz_xyy, g_x_0_xxxxz_xyz, g_x_0_xxxxz_xzz, g_x_0_xxxxz_yyy, g_x_0_xxxxz_yyz, g_x_0_xxxxz_yzz, g_x_0_xxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_xxx[k] = -g_x_0_xxxx_xxx[k] * ab_z + g_x_0_xxxx_xxxz[k];

                g_x_0_xxxxz_xxy[k] = -g_x_0_xxxx_xxy[k] * ab_z + g_x_0_xxxx_xxyz[k];

                g_x_0_xxxxz_xxz[k] = -g_x_0_xxxx_xxz[k] * ab_z + g_x_0_xxxx_xxzz[k];

                g_x_0_xxxxz_xyy[k] = -g_x_0_xxxx_xyy[k] * ab_z + g_x_0_xxxx_xyyz[k];

                g_x_0_xxxxz_xyz[k] = -g_x_0_xxxx_xyz[k] * ab_z + g_x_0_xxxx_xyzz[k];

                g_x_0_xxxxz_xzz[k] = -g_x_0_xxxx_xzz[k] * ab_z + g_x_0_xxxx_xzzz[k];

                g_x_0_xxxxz_yyy[k] = -g_x_0_xxxx_yyy[k] * ab_z + g_x_0_xxxx_yyyz[k];

                g_x_0_xxxxz_yyz[k] = -g_x_0_xxxx_yyz[k] * ab_z + g_x_0_xxxx_yyzz[k];

                g_x_0_xxxxz_yzz[k] = -g_x_0_xxxx_yzz[k] * ab_z + g_x_0_xxxx_yzzz[k];

                g_x_0_xxxxz_zzz[k] = -g_x_0_xxxx_zzz[k] * ab_z + g_x_0_xxxx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxy_xxx, g_x_0_xxxy_xxxy, g_x_0_xxxy_xxy, g_x_0_xxxy_xxyy, g_x_0_xxxy_xxyz, g_x_0_xxxy_xxz, g_x_0_xxxy_xyy, g_x_0_xxxy_xyyy, g_x_0_xxxy_xyyz, g_x_0_xxxy_xyz, g_x_0_xxxy_xyzz, g_x_0_xxxy_xzz, g_x_0_xxxy_yyy, g_x_0_xxxy_yyyy, g_x_0_xxxy_yyyz, g_x_0_xxxy_yyz, g_x_0_xxxy_yyzz, g_x_0_xxxy_yzz, g_x_0_xxxy_yzzz, g_x_0_xxxy_zzz, g_x_0_xxxyy_xxx, g_x_0_xxxyy_xxy, g_x_0_xxxyy_xxz, g_x_0_xxxyy_xyy, g_x_0_xxxyy_xyz, g_x_0_xxxyy_xzz, g_x_0_xxxyy_yyy, g_x_0_xxxyy_yyz, g_x_0_xxxyy_yzz, g_x_0_xxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_xxx[k] = -g_x_0_xxxy_xxx[k] * ab_y + g_x_0_xxxy_xxxy[k];

                g_x_0_xxxyy_xxy[k] = -g_x_0_xxxy_xxy[k] * ab_y + g_x_0_xxxy_xxyy[k];

                g_x_0_xxxyy_xxz[k] = -g_x_0_xxxy_xxz[k] * ab_y + g_x_0_xxxy_xxyz[k];

                g_x_0_xxxyy_xyy[k] = -g_x_0_xxxy_xyy[k] * ab_y + g_x_0_xxxy_xyyy[k];

                g_x_0_xxxyy_xyz[k] = -g_x_0_xxxy_xyz[k] * ab_y + g_x_0_xxxy_xyyz[k];

                g_x_0_xxxyy_xzz[k] = -g_x_0_xxxy_xzz[k] * ab_y + g_x_0_xxxy_xyzz[k];

                g_x_0_xxxyy_yyy[k] = -g_x_0_xxxy_yyy[k] * ab_y + g_x_0_xxxy_yyyy[k];

                g_x_0_xxxyy_yyz[k] = -g_x_0_xxxy_yyz[k] * ab_y + g_x_0_xxxy_yyyz[k];

                g_x_0_xxxyy_yzz[k] = -g_x_0_xxxy_yzz[k] * ab_y + g_x_0_xxxy_yyzz[k];

                g_x_0_xxxyy_zzz[k] = -g_x_0_xxxy_zzz[k] * ab_y + g_x_0_xxxy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyz_xxx, g_x_0_xxxyz_xxy, g_x_0_xxxyz_xxz, g_x_0_xxxyz_xyy, g_x_0_xxxyz_xyz, g_x_0_xxxyz_xzz, g_x_0_xxxyz_yyy, g_x_0_xxxyz_yyz, g_x_0_xxxyz_yzz, g_x_0_xxxyz_zzz, g_x_0_xxxz_xxx, g_x_0_xxxz_xxxy, g_x_0_xxxz_xxy, g_x_0_xxxz_xxyy, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxz, g_x_0_xxxz_xyy, g_x_0_xxxz_xyyy, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xzz, g_x_0_xxxz_yyy, g_x_0_xxxz_yyyy, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_xxx[k] = -g_x_0_xxxz_xxx[k] * ab_y + g_x_0_xxxz_xxxy[k];

                g_x_0_xxxyz_xxy[k] = -g_x_0_xxxz_xxy[k] * ab_y + g_x_0_xxxz_xxyy[k];

                g_x_0_xxxyz_xxz[k] = -g_x_0_xxxz_xxz[k] * ab_y + g_x_0_xxxz_xxyz[k];

                g_x_0_xxxyz_xyy[k] = -g_x_0_xxxz_xyy[k] * ab_y + g_x_0_xxxz_xyyy[k];

                g_x_0_xxxyz_xyz[k] = -g_x_0_xxxz_xyz[k] * ab_y + g_x_0_xxxz_xyyz[k];

                g_x_0_xxxyz_xzz[k] = -g_x_0_xxxz_xzz[k] * ab_y + g_x_0_xxxz_xyzz[k];

                g_x_0_xxxyz_yyy[k] = -g_x_0_xxxz_yyy[k] * ab_y + g_x_0_xxxz_yyyy[k];

                g_x_0_xxxyz_yyz[k] = -g_x_0_xxxz_yyz[k] * ab_y + g_x_0_xxxz_yyyz[k];

                g_x_0_xxxyz_yzz[k] = -g_x_0_xxxz_yzz[k] * ab_y + g_x_0_xxxz_yyzz[k];

                g_x_0_xxxyz_zzz[k] = -g_x_0_xxxz_zzz[k] * ab_y + g_x_0_xxxz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxz_xxx, g_x_0_xxxz_xxxz, g_x_0_xxxz_xxy, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxz, g_x_0_xxxz_xxzz, g_x_0_xxxz_xyy, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xzz, g_x_0_xxxz_xzzz, g_x_0_xxxz_yyy, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_zzz, g_x_0_xxxz_zzzz, g_x_0_xxxzz_xxx, g_x_0_xxxzz_xxy, g_x_0_xxxzz_xxz, g_x_0_xxxzz_xyy, g_x_0_xxxzz_xyz, g_x_0_xxxzz_xzz, g_x_0_xxxzz_yyy, g_x_0_xxxzz_yyz, g_x_0_xxxzz_yzz, g_x_0_xxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_xxx[k] = -g_x_0_xxxz_xxx[k] * ab_z + g_x_0_xxxz_xxxz[k];

                g_x_0_xxxzz_xxy[k] = -g_x_0_xxxz_xxy[k] * ab_z + g_x_0_xxxz_xxyz[k];

                g_x_0_xxxzz_xxz[k] = -g_x_0_xxxz_xxz[k] * ab_z + g_x_0_xxxz_xxzz[k];

                g_x_0_xxxzz_xyy[k] = -g_x_0_xxxz_xyy[k] * ab_z + g_x_0_xxxz_xyyz[k];

                g_x_0_xxxzz_xyz[k] = -g_x_0_xxxz_xyz[k] * ab_z + g_x_0_xxxz_xyzz[k];

                g_x_0_xxxzz_xzz[k] = -g_x_0_xxxz_xzz[k] * ab_z + g_x_0_xxxz_xzzz[k];

                g_x_0_xxxzz_yyy[k] = -g_x_0_xxxz_yyy[k] * ab_z + g_x_0_xxxz_yyyz[k];

                g_x_0_xxxzz_yyz[k] = -g_x_0_xxxz_yyz[k] * ab_z + g_x_0_xxxz_yyzz[k];

                g_x_0_xxxzz_yzz[k] = -g_x_0_xxxz_yzz[k] * ab_z + g_x_0_xxxz_yzzz[k];

                g_x_0_xxxzz_zzz[k] = -g_x_0_xxxz_zzz[k] * ab_z + g_x_0_xxxz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyy_xxx, g_x_0_xxyy_xxxy, g_x_0_xxyy_xxy, g_x_0_xxyy_xxyy, g_x_0_xxyy_xxyz, g_x_0_xxyy_xxz, g_x_0_xxyy_xyy, g_x_0_xxyy_xyyy, g_x_0_xxyy_xyyz, g_x_0_xxyy_xyz, g_x_0_xxyy_xyzz, g_x_0_xxyy_xzz, g_x_0_xxyy_yyy, g_x_0_xxyy_yyyy, g_x_0_xxyy_yyyz, g_x_0_xxyy_yyz, g_x_0_xxyy_yyzz, g_x_0_xxyy_yzz, g_x_0_xxyy_yzzz, g_x_0_xxyy_zzz, g_x_0_xxyyy_xxx, g_x_0_xxyyy_xxy, g_x_0_xxyyy_xxz, g_x_0_xxyyy_xyy, g_x_0_xxyyy_xyz, g_x_0_xxyyy_xzz, g_x_0_xxyyy_yyy, g_x_0_xxyyy_yyz, g_x_0_xxyyy_yzz, g_x_0_xxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_xxx[k] = -g_x_0_xxyy_xxx[k] * ab_y + g_x_0_xxyy_xxxy[k];

                g_x_0_xxyyy_xxy[k] = -g_x_0_xxyy_xxy[k] * ab_y + g_x_0_xxyy_xxyy[k];

                g_x_0_xxyyy_xxz[k] = -g_x_0_xxyy_xxz[k] * ab_y + g_x_0_xxyy_xxyz[k];

                g_x_0_xxyyy_xyy[k] = -g_x_0_xxyy_xyy[k] * ab_y + g_x_0_xxyy_xyyy[k];

                g_x_0_xxyyy_xyz[k] = -g_x_0_xxyy_xyz[k] * ab_y + g_x_0_xxyy_xyyz[k];

                g_x_0_xxyyy_xzz[k] = -g_x_0_xxyy_xzz[k] * ab_y + g_x_0_xxyy_xyzz[k];

                g_x_0_xxyyy_yyy[k] = -g_x_0_xxyy_yyy[k] * ab_y + g_x_0_xxyy_yyyy[k];

                g_x_0_xxyyy_yyz[k] = -g_x_0_xxyy_yyz[k] * ab_y + g_x_0_xxyy_yyyz[k];

                g_x_0_xxyyy_yzz[k] = -g_x_0_xxyy_yzz[k] * ab_y + g_x_0_xxyy_yyzz[k];

                g_x_0_xxyyy_zzz[k] = -g_x_0_xxyy_zzz[k] * ab_y + g_x_0_xxyy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyz_xxx, g_x_0_xxyyz_xxy, g_x_0_xxyyz_xxz, g_x_0_xxyyz_xyy, g_x_0_xxyyz_xyz, g_x_0_xxyyz_xzz, g_x_0_xxyyz_yyy, g_x_0_xxyyz_yyz, g_x_0_xxyyz_yzz, g_x_0_xxyyz_zzz, g_x_0_xxyz_xxx, g_x_0_xxyz_xxxy, g_x_0_xxyz_xxy, g_x_0_xxyz_xxyy, g_x_0_xxyz_xxyz, g_x_0_xxyz_xxz, g_x_0_xxyz_xyy, g_x_0_xxyz_xyyy, g_x_0_xxyz_xyyz, g_x_0_xxyz_xyz, g_x_0_xxyz_xyzz, g_x_0_xxyz_xzz, g_x_0_xxyz_yyy, g_x_0_xxyz_yyyy, g_x_0_xxyz_yyyz, g_x_0_xxyz_yyz, g_x_0_xxyz_yyzz, g_x_0_xxyz_yzz, g_x_0_xxyz_yzzz, g_x_0_xxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_xxx[k] = -g_x_0_xxyz_xxx[k] * ab_y + g_x_0_xxyz_xxxy[k];

                g_x_0_xxyyz_xxy[k] = -g_x_0_xxyz_xxy[k] * ab_y + g_x_0_xxyz_xxyy[k];

                g_x_0_xxyyz_xxz[k] = -g_x_0_xxyz_xxz[k] * ab_y + g_x_0_xxyz_xxyz[k];

                g_x_0_xxyyz_xyy[k] = -g_x_0_xxyz_xyy[k] * ab_y + g_x_0_xxyz_xyyy[k];

                g_x_0_xxyyz_xyz[k] = -g_x_0_xxyz_xyz[k] * ab_y + g_x_0_xxyz_xyyz[k];

                g_x_0_xxyyz_xzz[k] = -g_x_0_xxyz_xzz[k] * ab_y + g_x_0_xxyz_xyzz[k];

                g_x_0_xxyyz_yyy[k] = -g_x_0_xxyz_yyy[k] * ab_y + g_x_0_xxyz_yyyy[k];

                g_x_0_xxyyz_yyz[k] = -g_x_0_xxyz_yyz[k] * ab_y + g_x_0_xxyz_yyyz[k];

                g_x_0_xxyyz_yzz[k] = -g_x_0_xxyz_yzz[k] * ab_y + g_x_0_xxyz_yyzz[k];

                g_x_0_xxyyz_zzz[k] = -g_x_0_xxyz_zzz[k] * ab_y + g_x_0_xxyz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyzz_xxx, g_x_0_xxyzz_xxy, g_x_0_xxyzz_xxz, g_x_0_xxyzz_xyy, g_x_0_xxyzz_xyz, g_x_0_xxyzz_xzz, g_x_0_xxyzz_yyy, g_x_0_xxyzz_yyz, g_x_0_xxyzz_yzz, g_x_0_xxyzz_zzz, g_x_0_xxzz_xxx, g_x_0_xxzz_xxxy, g_x_0_xxzz_xxy, g_x_0_xxzz_xxyy, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxz, g_x_0_xxzz_xyy, g_x_0_xxzz_xyyy, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xzz, g_x_0_xxzz_yyy, g_x_0_xxzz_yyyy, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_xxx[k] = -g_x_0_xxzz_xxx[k] * ab_y + g_x_0_xxzz_xxxy[k];

                g_x_0_xxyzz_xxy[k] = -g_x_0_xxzz_xxy[k] * ab_y + g_x_0_xxzz_xxyy[k];

                g_x_0_xxyzz_xxz[k] = -g_x_0_xxzz_xxz[k] * ab_y + g_x_0_xxzz_xxyz[k];

                g_x_0_xxyzz_xyy[k] = -g_x_0_xxzz_xyy[k] * ab_y + g_x_0_xxzz_xyyy[k];

                g_x_0_xxyzz_xyz[k] = -g_x_0_xxzz_xyz[k] * ab_y + g_x_0_xxzz_xyyz[k];

                g_x_0_xxyzz_xzz[k] = -g_x_0_xxzz_xzz[k] * ab_y + g_x_0_xxzz_xyzz[k];

                g_x_0_xxyzz_yyy[k] = -g_x_0_xxzz_yyy[k] * ab_y + g_x_0_xxzz_yyyy[k];

                g_x_0_xxyzz_yyz[k] = -g_x_0_xxzz_yyz[k] * ab_y + g_x_0_xxzz_yyyz[k];

                g_x_0_xxyzz_yzz[k] = -g_x_0_xxzz_yzz[k] * ab_y + g_x_0_xxzz_yyzz[k];

                g_x_0_xxyzz_zzz[k] = -g_x_0_xxzz_zzz[k] * ab_y + g_x_0_xxzz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxzz_xxx, g_x_0_xxzz_xxxz, g_x_0_xxzz_xxy, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxz, g_x_0_xxzz_xxzz, g_x_0_xxzz_xyy, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xzz, g_x_0_xxzz_xzzz, g_x_0_xxzz_yyy, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_zzz, g_x_0_xxzz_zzzz, g_x_0_xxzzz_xxx, g_x_0_xxzzz_xxy, g_x_0_xxzzz_xxz, g_x_0_xxzzz_xyy, g_x_0_xxzzz_xyz, g_x_0_xxzzz_xzz, g_x_0_xxzzz_yyy, g_x_0_xxzzz_yyz, g_x_0_xxzzz_yzz, g_x_0_xxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_xxx[k] = -g_x_0_xxzz_xxx[k] * ab_z + g_x_0_xxzz_xxxz[k];

                g_x_0_xxzzz_xxy[k] = -g_x_0_xxzz_xxy[k] * ab_z + g_x_0_xxzz_xxyz[k];

                g_x_0_xxzzz_xxz[k] = -g_x_0_xxzz_xxz[k] * ab_z + g_x_0_xxzz_xxzz[k];

                g_x_0_xxzzz_xyy[k] = -g_x_0_xxzz_xyy[k] * ab_z + g_x_0_xxzz_xyyz[k];

                g_x_0_xxzzz_xyz[k] = -g_x_0_xxzz_xyz[k] * ab_z + g_x_0_xxzz_xyzz[k];

                g_x_0_xxzzz_xzz[k] = -g_x_0_xxzz_xzz[k] * ab_z + g_x_0_xxzz_xzzz[k];

                g_x_0_xxzzz_yyy[k] = -g_x_0_xxzz_yyy[k] * ab_z + g_x_0_xxzz_yyyz[k];

                g_x_0_xxzzz_yyz[k] = -g_x_0_xxzz_yyz[k] * ab_z + g_x_0_xxzz_yyzz[k];

                g_x_0_xxzzz_yzz[k] = -g_x_0_xxzz_yzz[k] * ab_z + g_x_0_xxzz_yzzz[k];

                g_x_0_xxzzz_zzz[k] = -g_x_0_xxzz_zzz[k] * ab_z + g_x_0_xxzz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyy_xxx, g_x_0_xyyy_xxxy, g_x_0_xyyy_xxy, g_x_0_xyyy_xxyy, g_x_0_xyyy_xxyz, g_x_0_xyyy_xxz, g_x_0_xyyy_xyy, g_x_0_xyyy_xyyy, g_x_0_xyyy_xyyz, g_x_0_xyyy_xyz, g_x_0_xyyy_xyzz, g_x_0_xyyy_xzz, g_x_0_xyyy_yyy, g_x_0_xyyy_yyyy, g_x_0_xyyy_yyyz, g_x_0_xyyy_yyz, g_x_0_xyyy_yyzz, g_x_0_xyyy_yzz, g_x_0_xyyy_yzzz, g_x_0_xyyy_zzz, g_x_0_xyyyy_xxx, g_x_0_xyyyy_xxy, g_x_0_xyyyy_xxz, g_x_0_xyyyy_xyy, g_x_0_xyyyy_xyz, g_x_0_xyyyy_xzz, g_x_0_xyyyy_yyy, g_x_0_xyyyy_yyz, g_x_0_xyyyy_yzz, g_x_0_xyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_xxx[k] = -g_x_0_xyyy_xxx[k] * ab_y + g_x_0_xyyy_xxxy[k];

                g_x_0_xyyyy_xxy[k] = -g_x_0_xyyy_xxy[k] * ab_y + g_x_0_xyyy_xxyy[k];

                g_x_0_xyyyy_xxz[k] = -g_x_0_xyyy_xxz[k] * ab_y + g_x_0_xyyy_xxyz[k];

                g_x_0_xyyyy_xyy[k] = -g_x_0_xyyy_xyy[k] * ab_y + g_x_0_xyyy_xyyy[k];

                g_x_0_xyyyy_xyz[k] = -g_x_0_xyyy_xyz[k] * ab_y + g_x_0_xyyy_xyyz[k];

                g_x_0_xyyyy_xzz[k] = -g_x_0_xyyy_xzz[k] * ab_y + g_x_0_xyyy_xyzz[k];

                g_x_0_xyyyy_yyy[k] = -g_x_0_xyyy_yyy[k] * ab_y + g_x_0_xyyy_yyyy[k];

                g_x_0_xyyyy_yyz[k] = -g_x_0_xyyy_yyz[k] * ab_y + g_x_0_xyyy_yyyz[k];

                g_x_0_xyyyy_yzz[k] = -g_x_0_xyyy_yzz[k] * ab_y + g_x_0_xyyy_yyzz[k];

                g_x_0_xyyyy_zzz[k] = -g_x_0_xyyy_zzz[k] * ab_y + g_x_0_xyyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyz_xxx, g_x_0_xyyyz_xxy, g_x_0_xyyyz_xxz, g_x_0_xyyyz_xyy, g_x_0_xyyyz_xyz, g_x_0_xyyyz_xzz, g_x_0_xyyyz_yyy, g_x_0_xyyyz_yyz, g_x_0_xyyyz_yzz, g_x_0_xyyyz_zzz, g_x_0_xyyz_xxx, g_x_0_xyyz_xxxy, g_x_0_xyyz_xxy, g_x_0_xyyz_xxyy, g_x_0_xyyz_xxyz, g_x_0_xyyz_xxz, g_x_0_xyyz_xyy, g_x_0_xyyz_xyyy, g_x_0_xyyz_xyyz, g_x_0_xyyz_xyz, g_x_0_xyyz_xyzz, g_x_0_xyyz_xzz, g_x_0_xyyz_yyy, g_x_0_xyyz_yyyy, g_x_0_xyyz_yyyz, g_x_0_xyyz_yyz, g_x_0_xyyz_yyzz, g_x_0_xyyz_yzz, g_x_0_xyyz_yzzz, g_x_0_xyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_xxx[k] = -g_x_0_xyyz_xxx[k] * ab_y + g_x_0_xyyz_xxxy[k];

                g_x_0_xyyyz_xxy[k] = -g_x_0_xyyz_xxy[k] * ab_y + g_x_0_xyyz_xxyy[k];

                g_x_0_xyyyz_xxz[k] = -g_x_0_xyyz_xxz[k] * ab_y + g_x_0_xyyz_xxyz[k];

                g_x_0_xyyyz_xyy[k] = -g_x_0_xyyz_xyy[k] * ab_y + g_x_0_xyyz_xyyy[k];

                g_x_0_xyyyz_xyz[k] = -g_x_0_xyyz_xyz[k] * ab_y + g_x_0_xyyz_xyyz[k];

                g_x_0_xyyyz_xzz[k] = -g_x_0_xyyz_xzz[k] * ab_y + g_x_0_xyyz_xyzz[k];

                g_x_0_xyyyz_yyy[k] = -g_x_0_xyyz_yyy[k] * ab_y + g_x_0_xyyz_yyyy[k];

                g_x_0_xyyyz_yyz[k] = -g_x_0_xyyz_yyz[k] * ab_y + g_x_0_xyyz_yyyz[k];

                g_x_0_xyyyz_yzz[k] = -g_x_0_xyyz_yzz[k] * ab_y + g_x_0_xyyz_yyzz[k];

                g_x_0_xyyyz_zzz[k] = -g_x_0_xyyz_zzz[k] * ab_y + g_x_0_xyyz_yzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyzz_xxx, g_x_0_xyyzz_xxy, g_x_0_xyyzz_xxz, g_x_0_xyyzz_xyy, g_x_0_xyyzz_xyz, g_x_0_xyyzz_xzz, g_x_0_xyyzz_yyy, g_x_0_xyyzz_yyz, g_x_0_xyyzz_yzz, g_x_0_xyyzz_zzz, g_x_0_xyzz_xxx, g_x_0_xyzz_xxxy, g_x_0_xyzz_xxy, g_x_0_xyzz_xxyy, g_x_0_xyzz_xxyz, g_x_0_xyzz_xxz, g_x_0_xyzz_xyy, g_x_0_xyzz_xyyy, g_x_0_xyzz_xyyz, g_x_0_xyzz_xyz, g_x_0_xyzz_xyzz, g_x_0_xyzz_xzz, g_x_0_xyzz_yyy, g_x_0_xyzz_yyyy, g_x_0_xyzz_yyyz, g_x_0_xyzz_yyz, g_x_0_xyzz_yyzz, g_x_0_xyzz_yzz, g_x_0_xyzz_yzzz, g_x_0_xyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_xxx[k] = -g_x_0_xyzz_xxx[k] * ab_y + g_x_0_xyzz_xxxy[k];

                g_x_0_xyyzz_xxy[k] = -g_x_0_xyzz_xxy[k] * ab_y + g_x_0_xyzz_xxyy[k];

                g_x_0_xyyzz_xxz[k] = -g_x_0_xyzz_xxz[k] * ab_y + g_x_0_xyzz_xxyz[k];

                g_x_0_xyyzz_xyy[k] = -g_x_0_xyzz_xyy[k] * ab_y + g_x_0_xyzz_xyyy[k];

                g_x_0_xyyzz_xyz[k] = -g_x_0_xyzz_xyz[k] * ab_y + g_x_0_xyzz_xyyz[k];

                g_x_0_xyyzz_xzz[k] = -g_x_0_xyzz_xzz[k] * ab_y + g_x_0_xyzz_xyzz[k];

                g_x_0_xyyzz_yyy[k] = -g_x_0_xyzz_yyy[k] * ab_y + g_x_0_xyzz_yyyy[k];

                g_x_0_xyyzz_yyz[k] = -g_x_0_xyzz_yyz[k] * ab_y + g_x_0_xyzz_yyyz[k];

                g_x_0_xyyzz_yzz[k] = -g_x_0_xyzz_yzz[k] * ab_y + g_x_0_xyzz_yyzz[k];

                g_x_0_xyyzz_zzz[k] = -g_x_0_xyzz_zzz[k] * ab_y + g_x_0_xyzz_yzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzzz_xxx, g_x_0_xyzzz_xxy, g_x_0_xyzzz_xxz, g_x_0_xyzzz_xyy, g_x_0_xyzzz_xyz, g_x_0_xyzzz_xzz, g_x_0_xyzzz_yyy, g_x_0_xyzzz_yyz, g_x_0_xyzzz_yzz, g_x_0_xyzzz_zzz, g_x_0_xzzz_xxx, g_x_0_xzzz_xxxy, g_x_0_xzzz_xxy, g_x_0_xzzz_xxyy, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxz, g_x_0_xzzz_xyy, g_x_0_xzzz_xyyy, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xzz, g_x_0_xzzz_yyy, g_x_0_xzzz_yyyy, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_xxx[k] = -g_x_0_xzzz_xxx[k] * ab_y + g_x_0_xzzz_xxxy[k];

                g_x_0_xyzzz_xxy[k] = -g_x_0_xzzz_xxy[k] * ab_y + g_x_0_xzzz_xxyy[k];

                g_x_0_xyzzz_xxz[k] = -g_x_0_xzzz_xxz[k] * ab_y + g_x_0_xzzz_xxyz[k];

                g_x_0_xyzzz_xyy[k] = -g_x_0_xzzz_xyy[k] * ab_y + g_x_0_xzzz_xyyy[k];

                g_x_0_xyzzz_xyz[k] = -g_x_0_xzzz_xyz[k] * ab_y + g_x_0_xzzz_xyyz[k];

                g_x_0_xyzzz_xzz[k] = -g_x_0_xzzz_xzz[k] * ab_y + g_x_0_xzzz_xyzz[k];

                g_x_0_xyzzz_yyy[k] = -g_x_0_xzzz_yyy[k] * ab_y + g_x_0_xzzz_yyyy[k];

                g_x_0_xyzzz_yyz[k] = -g_x_0_xzzz_yyz[k] * ab_y + g_x_0_xzzz_yyyz[k];

                g_x_0_xyzzz_yzz[k] = -g_x_0_xzzz_yzz[k] * ab_y + g_x_0_xzzz_yyzz[k];

                g_x_0_xyzzz_zzz[k] = -g_x_0_xzzz_zzz[k] * ab_y + g_x_0_xzzz_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzzz_xxx, g_x_0_xzzz_xxxz, g_x_0_xzzz_xxy, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxz, g_x_0_xzzz_xxzz, g_x_0_xzzz_xyy, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xzz, g_x_0_xzzz_xzzz, g_x_0_xzzz_yyy, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_zzz, g_x_0_xzzz_zzzz, g_x_0_xzzzz_xxx, g_x_0_xzzzz_xxy, g_x_0_xzzzz_xxz, g_x_0_xzzzz_xyy, g_x_0_xzzzz_xyz, g_x_0_xzzzz_xzz, g_x_0_xzzzz_yyy, g_x_0_xzzzz_yyz, g_x_0_xzzzz_yzz, g_x_0_xzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_xxx[k] = -g_x_0_xzzz_xxx[k] * ab_z + g_x_0_xzzz_xxxz[k];

                g_x_0_xzzzz_xxy[k] = -g_x_0_xzzz_xxy[k] * ab_z + g_x_0_xzzz_xxyz[k];

                g_x_0_xzzzz_xxz[k] = -g_x_0_xzzz_xxz[k] * ab_z + g_x_0_xzzz_xxzz[k];

                g_x_0_xzzzz_xyy[k] = -g_x_0_xzzz_xyy[k] * ab_z + g_x_0_xzzz_xyyz[k];

                g_x_0_xzzzz_xyz[k] = -g_x_0_xzzz_xyz[k] * ab_z + g_x_0_xzzz_xyzz[k];

                g_x_0_xzzzz_xzz[k] = -g_x_0_xzzz_xzz[k] * ab_z + g_x_0_xzzz_xzzz[k];

                g_x_0_xzzzz_yyy[k] = -g_x_0_xzzz_yyy[k] * ab_z + g_x_0_xzzz_yyyz[k];

                g_x_0_xzzzz_yyz[k] = -g_x_0_xzzz_yyz[k] * ab_z + g_x_0_xzzz_yyzz[k];

                g_x_0_xzzzz_yzz[k] = -g_x_0_xzzz_yzz[k] * ab_z + g_x_0_xzzz_yzzz[k];

                g_x_0_xzzzz_zzz[k] = -g_x_0_xzzz_zzz[k] * ab_z + g_x_0_xzzz_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyy_xxx, g_x_0_yyyy_xxxy, g_x_0_yyyy_xxy, g_x_0_yyyy_xxyy, g_x_0_yyyy_xxyz, g_x_0_yyyy_xxz, g_x_0_yyyy_xyy, g_x_0_yyyy_xyyy, g_x_0_yyyy_xyyz, g_x_0_yyyy_xyz, g_x_0_yyyy_xyzz, g_x_0_yyyy_xzz, g_x_0_yyyy_yyy, g_x_0_yyyy_yyyy, g_x_0_yyyy_yyyz, g_x_0_yyyy_yyz, g_x_0_yyyy_yyzz, g_x_0_yyyy_yzz, g_x_0_yyyy_yzzz, g_x_0_yyyy_zzz, g_x_0_yyyyy_xxx, g_x_0_yyyyy_xxy, g_x_0_yyyyy_xxz, g_x_0_yyyyy_xyy, g_x_0_yyyyy_xyz, g_x_0_yyyyy_xzz, g_x_0_yyyyy_yyy, g_x_0_yyyyy_yyz, g_x_0_yyyyy_yzz, g_x_0_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_xxx[k] = -g_x_0_yyyy_xxx[k] * ab_y + g_x_0_yyyy_xxxy[k];

                g_x_0_yyyyy_xxy[k] = -g_x_0_yyyy_xxy[k] * ab_y + g_x_0_yyyy_xxyy[k];

                g_x_0_yyyyy_xxz[k] = -g_x_0_yyyy_xxz[k] * ab_y + g_x_0_yyyy_xxyz[k];

                g_x_0_yyyyy_xyy[k] = -g_x_0_yyyy_xyy[k] * ab_y + g_x_0_yyyy_xyyy[k];

                g_x_0_yyyyy_xyz[k] = -g_x_0_yyyy_xyz[k] * ab_y + g_x_0_yyyy_xyyz[k];

                g_x_0_yyyyy_xzz[k] = -g_x_0_yyyy_xzz[k] * ab_y + g_x_0_yyyy_xyzz[k];

                g_x_0_yyyyy_yyy[k] = -g_x_0_yyyy_yyy[k] * ab_y + g_x_0_yyyy_yyyy[k];

                g_x_0_yyyyy_yyz[k] = -g_x_0_yyyy_yyz[k] * ab_y + g_x_0_yyyy_yyyz[k];

                g_x_0_yyyyy_yzz[k] = -g_x_0_yyyy_yzz[k] * ab_y + g_x_0_yyyy_yyzz[k];

                g_x_0_yyyyy_zzz[k] = -g_x_0_yyyy_zzz[k] * ab_y + g_x_0_yyyy_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyz_xxx, g_x_0_yyyyz_xxy, g_x_0_yyyyz_xxz, g_x_0_yyyyz_xyy, g_x_0_yyyyz_xyz, g_x_0_yyyyz_xzz, g_x_0_yyyyz_yyy, g_x_0_yyyyz_yyz, g_x_0_yyyyz_yzz, g_x_0_yyyyz_zzz, g_x_0_yyyz_xxx, g_x_0_yyyz_xxxy, g_x_0_yyyz_xxy, g_x_0_yyyz_xxyy, g_x_0_yyyz_xxyz, g_x_0_yyyz_xxz, g_x_0_yyyz_xyy, g_x_0_yyyz_xyyy, g_x_0_yyyz_xyyz, g_x_0_yyyz_xyz, g_x_0_yyyz_xyzz, g_x_0_yyyz_xzz, g_x_0_yyyz_yyy, g_x_0_yyyz_yyyy, g_x_0_yyyz_yyyz, g_x_0_yyyz_yyz, g_x_0_yyyz_yyzz, g_x_0_yyyz_yzz, g_x_0_yyyz_yzzz, g_x_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_xxx[k] = -g_x_0_yyyz_xxx[k] * ab_y + g_x_0_yyyz_xxxy[k];

                g_x_0_yyyyz_xxy[k] = -g_x_0_yyyz_xxy[k] * ab_y + g_x_0_yyyz_xxyy[k];

                g_x_0_yyyyz_xxz[k] = -g_x_0_yyyz_xxz[k] * ab_y + g_x_0_yyyz_xxyz[k];

                g_x_0_yyyyz_xyy[k] = -g_x_0_yyyz_xyy[k] * ab_y + g_x_0_yyyz_xyyy[k];

                g_x_0_yyyyz_xyz[k] = -g_x_0_yyyz_xyz[k] * ab_y + g_x_0_yyyz_xyyz[k];

                g_x_0_yyyyz_xzz[k] = -g_x_0_yyyz_xzz[k] * ab_y + g_x_0_yyyz_xyzz[k];

                g_x_0_yyyyz_yyy[k] = -g_x_0_yyyz_yyy[k] * ab_y + g_x_0_yyyz_yyyy[k];

                g_x_0_yyyyz_yyz[k] = -g_x_0_yyyz_yyz[k] * ab_y + g_x_0_yyyz_yyyz[k];

                g_x_0_yyyyz_yzz[k] = -g_x_0_yyyz_yzz[k] * ab_y + g_x_0_yyyz_yyzz[k];

                g_x_0_yyyyz_zzz[k] = -g_x_0_yyyz_zzz[k] * ab_y + g_x_0_yyyz_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyzz_xxx, g_x_0_yyyzz_xxy, g_x_0_yyyzz_xxz, g_x_0_yyyzz_xyy, g_x_0_yyyzz_xyz, g_x_0_yyyzz_xzz, g_x_0_yyyzz_yyy, g_x_0_yyyzz_yyz, g_x_0_yyyzz_yzz, g_x_0_yyyzz_zzz, g_x_0_yyzz_xxx, g_x_0_yyzz_xxxy, g_x_0_yyzz_xxy, g_x_0_yyzz_xxyy, g_x_0_yyzz_xxyz, g_x_0_yyzz_xxz, g_x_0_yyzz_xyy, g_x_0_yyzz_xyyy, g_x_0_yyzz_xyyz, g_x_0_yyzz_xyz, g_x_0_yyzz_xyzz, g_x_0_yyzz_xzz, g_x_0_yyzz_yyy, g_x_0_yyzz_yyyy, g_x_0_yyzz_yyyz, g_x_0_yyzz_yyz, g_x_0_yyzz_yyzz, g_x_0_yyzz_yzz, g_x_0_yyzz_yzzz, g_x_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_xxx[k] = -g_x_0_yyzz_xxx[k] * ab_y + g_x_0_yyzz_xxxy[k];

                g_x_0_yyyzz_xxy[k] = -g_x_0_yyzz_xxy[k] * ab_y + g_x_0_yyzz_xxyy[k];

                g_x_0_yyyzz_xxz[k] = -g_x_0_yyzz_xxz[k] * ab_y + g_x_0_yyzz_xxyz[k];

                g_x_0_yyyzz_xyy[k] = -g_x_0_yyzz_xyy[k] * ab_y + g_x_0_yyzz_xyyy[k];

                g_x_0_yyyzz_xyz[k] = -g_x_0_yyzz_xyz[k] * ab_y + g_x_0_yyzz_xyyz[k];

                g_x_0_yyyzz_xzz[k] = -g_x_0_yyzz_xzz[k] * ab_y + g_x_0_yyzz_xyzz[k];

                g_x_0_yyyzz_yyy[k] = -g_x_0_yyzz_yyy[k] * ab_y + g_x_0_yyzz_yyyy[k];

                g_x_0_yyyzz_yyz[k] = -g_x_0_yyzz_yyz[k] * ab_y + g_x_0_yyzz_yyyz[k];

                g_x_0_yyyzz_yzz[k] = -g_x_0_yyzz_yzz[k] * ab_y + g_x_0_yyzz_yyzz[k];

                g_x_0_yyyzz_zzz[k] = -g_x_0_yyzz_zzz[k] * ab_y + g_x_0_yyzz_yzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzzz_xxx, g_x_0_yyzzz_xxy, g_x_0_yyzzz_xxz, g_x_0_yyzzz_xyy, g_x_0_yyzzz_xyz, g_x_0_yyzzz_xzz, g_x_0_yyzzz_yyy, g_x_0_yyzzz_yyz, g_x_0_yyzzz_yzz, g_x_0_yyzzz_zzz, g_x_0_yzzz_xxx, g_x_0_yzzz_xxxy, g_x_0_yzzz_xxy, g_x_0_yzzz_xxyy, g_x_0_yzzz_xxyz, g_x_0_yzzz_xxz, g_x_0_yzzz_xyy, g_x_0_yzzz_xyyy, g_x_0_yzzz_xyyz, g_x_0_yzzz_xyz, g_x_0_yzzz_xyzz, g_x_0_yzzz_xzz, g_x_0_yzzz_yyy, g_x_0_yzzz_yyyy, g_x_0_yzzz_yyyz, g_x_0_yzzz_yyz, g_x_0_yzzz_yyzz, g_x_0_yzzz_yzz, g_x_0_yzzz_yzzz, g_x_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_xxx[k] = -g_x_0_yzzz_xxx[k] * ab_y + g_x_0_yzzz_xxxy[k];

                g_x_0_yyzzz_xxy[k] = -g_x_0_yzzz_xxy[k] * ab_y + g_x_0_yzzz_xxyy[k];

                g_x_0_yyzzz_xxz[k] = -g_x_0_yzzz_xxz[k] * ab_y + g_x_0_yzzz_xxyz[k];

                g_x_0_yyzzz_xyy[k] = -g_x_0_yzzz_xyy[k] * ab_y + g_x_0_yzzz_xyyy[k];

                g_x_0_yyzzz_xyz[k] = -g_x_0_yzzz_xyz[k] * ab_y + g_x_0_yzzz_xyyz[k];

                g_x_0_yyzzz_xzz[k] = -g_x_0_yzzz_xzz[k] * ab_y + g_x_0_yzzz_xyzz[k];

                g_x_0_yyzzz_yyy[k] = -g_x_0_yzzz_yyy[k] * ab_y + g_x_0_yzzz_yyyy[k];

                g_x_0_yyzzz_yyz[k] = -g_x_0_yzzz_yyz[k] * ab_y + g_x_0_yzzz_yyyz[k];

                g_x_0_yyzzz_yzz[k] = -g_x_0_yzzz_yzz[k] * ab_y + g_x_0_yzzz_yyzz[k];

                g_x_0_yyzzz_zzz[k] = -g_x_0_yzzz_zzz[k] * ab_y + g_x_0_yzzz_yzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzzz_xxx, g_x_0_yzzzz_xxy, g_x_0_yzzzz_xxz, g_x_0_yzzzz_xyy, g_x_0_yzzzz_xyz, g_x_0_yzzzz_xzz, g_x_0_yzzzz_yyy, g_x_0_yzzzz_yyz, g_x_0_yzzzz_yzz, g_x_0_yzzzz_zzz, g_x_0_zzzz_xxx, g_x_0_zzzz_xxxy, g_x_0_zzzz_xxy, g_x_0_zzzz_xxyy, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxz, g_x_0_zzzz_xyy, g_x_0_zzzz_xyyy, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xzz, g_x_0_zzzz_yyy, g_x_0_zzzz_yyyy, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_xxx[k] = -g_x_0_zzzz_xxx[k] * ab_y + g_x_0_zzzz_xxxy[k];

                g_x_0_yzzzz_xxy[k] = -g_x_0_zzzz_xxy[k] * ab_y + g_x_0_zzzz_xxyy[k];

                g_x_0_yzzzz_xxz[k] = -g_x_0_zzzz_xxz[k] * ab_y + g_x_0_zzzz_xxyz[k];

                g_x_0_yzzzz_xyy[k] = -g_x_0_zzzz_xyy[k] * ab_y + g_x_0_zzzz_xyyy[k];

                g_x_0_yzzzz_xyz[k] = -g_x_0_zzzz_xyz[k] * ab_y + g_x_0_zzzz_xyyz[k];

                g_x_0_yzzzz_xzz[k] = -g_x_0_zzzz_xzz[k] * ab_y + g_x_0_zzzz_xyzz[k];

                g_x_0_yzzzz_yyy[k] = -g_x_0_zzzz_yyy[k] * ab_y + g_x_0_zzzz_yyyy[k];

                g_x_0_yzzzz_yyz[k] = -g_x_0_zzzz_yyz[k] * ab_y + g_x_0_zzzz_yyyz[k];

                g_x_0_yzzzz_yzz[k] = -g_x_0_zzzz_yzz[k] * ab_y + g_x_0_zzzz_yyzz[k];

                g_x_0_yzzzz_zzz[k] = -g_x_0_zzzz_zzz[k] * ab_y + g_x_0_zzzz_yzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzz_xxx, g_x_0_zzzz_xxxz, g_x_0_zzzz_xxy, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxz, g_x_0_zzzz_xxzz, g_x_0_zzzz_xyy, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xzz, g_x_0_zzzz_xzzz, g_x_0_zzzz_yyy, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_zzz, g_x_0_zzzz_zzzz, g_x_0_zzzzz_xxx, g_x_0_zzzzz_xxy, g_x_0_zzzzz_xxz, g_x_0_zzzzz_xyy, g_x_0_zzzzz_xyz, g_x_0_zzzzz_xzz, g_x_0_zzzzz_yyy, g_x_0_zzzzz_yyz, g_x_0_zzzzz_yzz, g_x_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_xxx[k] = -g_x_0_zzzz_xxx[k] * ab_z + g_x_0_zzzz_xxxz[k];

                g_x_0_zzzzz_xxy[k] = -g_x_0_zzzz_xxy[k] * ab_z + g_x_0_zzzz_xxyz[k];

                g_x_0_zzzzz_xxz[k] = -g_x_0_zzzz_xxz[k] * ab_z + g_x_0_zzzz_xxzz[k];

                g_x_0_zzzzz_xyy[k] = -g_x_0_zzzz_xyy[k] * ab_z + g_x_0_zzzz_xyyz[k];

                g_x_0_zzzzz_xyz[k] = -g_x_0_zzzz_xyz[k] * ab_z + g_x_0_zzzz_xyzz[k];

                g_x_0_zzzzz_xzz[k] = -g_x_0_zzzz_xzz[k] * ab_z + g_x_0_zzzz_xzzz[k];

                g_x_0_zzzzz_yyy[k] = -g_x_0_zzzz_yyy[k] * ab_z + g_x_0_zzzz_yyyz[k];

                g_x_0_zzzzz_yyz[k] = -g_x_0_zzzz_yyz[k] * ab_z + g_x_0_zzzz_yyzz[k];

                g_x_0_zzzzz_yzz[k] = -g_x_0_zzzz_yzz[k] * ab_z + g_x_0_zzzz_yzzz[k];

                g_x_0_zzzzz_zzz[k] = -g_x_0_zzzz_zzz[k] * ab_z + g_x_0_zzzz_zzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxx_xxx, g_y_0_xxxx_xxxx, g_y_0_xxxx_xxxy, g_y_0_xxxx_xxxz, g_y_0_xxxx_xxy, g_y_0_xxxx_xxyy, g_y_0_xxxx_xxyz, g_y_0_xxxx_xxz, g_y_0_xxxx_xxzz, g_y_0_xxxx_xyy, g_y_0_xxxx_xyyy, g_y_0_xxxx_xyyz, g_y_0_xxxx_xyz, g_y_0_xxxx_xyzz, g_y_0_xxxx_xzz, g_y_0_xxxx_xzzz, g_y_0_xxxx_yyy, g_y_0_xxxx_yyz, g_y_0_xxxx_yzz, g_y_0_xxxx_zzz, g_y_0_xxxxx_xxx, g_y_0_xxxxx_xxy, g_y_0_xxxxx_xxz, g_y_0_xxxxx_xyy, g_y_0_xxxxx_xyz, g_y_0_xxxxx_xzz, g_y_0_xxxxx_yyy, g_y_0_xxxxx_yyz, g_y_0_xxxxx_yzz, g_y_0_xxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_xxx[k] = -g_y_0_xxxx_xxx[k] * ab_x + g_y_0_xxxx_xxxx[k];

                g_y_0_xxxxx_xxy[k] = -g_y_0_xxxx_xxy[k] * ab_x + g_y_0_xxxx_xxxy[k];

                g_y_0_xxxxx_xxz[k] = -g_y_0_xxxx_xxz[k] * ab_x + g_y_0_xxxx_xxxz[k];

                g_y_0_xxxxx_xyy[k] = -g_y_0_xxxx_xyy[k] * ab_x + g_y_0_xxxx_xxyy[k];

                g_y_0_xxxxx_xyz[k] = -g_y_0_xxxx_xyz[k] * ab_x + g_y_0_xxxx_xxyz[k];

                g_y_0_xxxxx_xzz[k] = -g_y_0_xxxx_xzz[k] * ab_x + g_y_0_xxxx_xxzz[k];

                g_y_0_xxxxx_yyy[k] = -g_y_0_xxxx_yyy[k] * ab_x + g_y_0_xxxx_xyyy[k];

                g_y_0_xxxxx_yyz[k] = -g_y_0_xxxx_yyz[k] * ab_x + g_y_0_xxxx_xyyz[k];

                g_y_0_xxxxx_yzz[k] = -g_y_0_xxxx_yzz[k] * ab_x + g_y_0_xxxx_xyzz[k];

                g_y_0_xxxxx_zzz[k] = -g_y_0_xxxx_zzz[k] * ab_x + g_y_0_xxxx_xzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 223 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxy_xxx, g_y_0_xxxxy_xxy, g_y_0_xxxxy_xxz, g_y_0_xxxxy_xyy, g_y_0_xxxxy_xyz, g_y_0_xxxxy_xzz, g_y_0_xxxxy_yyy, g_y_0_xxxxy_yyz, g_y_0_xxxxy_yzz, g_y_0_xxxxy_zzz, g_y_0_xxxy_xxx, g_y_0_xxxy_xxxx, g_y_0_xxxy_xxxy, g_y_0_xxxy_xxxz, g_y_0_xxxy_xxy, g_y_0_xxxy_xxyy, g_y_0_xxxy_xxyz, g_y_0_xxxy_xxz, g_y_0_xxxy_xxzz, g_y_0_xxxy_xyy, g_y_0_xxxy_xyyy, g_y_0_xxxy_xyyz, g_y_0_xxxy_xyz, g_y_0_xxxy_xyzz, g_y_0_xxxy_xzz, g_y_0_xxxy_xzzz, g_y_0_xxxy_yyy, g_y_0_xxxy_yyz, g_y_0_xxxy_yzz, g_y_0_xxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_xxx[k] = -g_y_0_xxxy_xxx[k] * ab_x + g_y_0_xxxy_xxxx[k];

                g_y_0_xxxxy_xxy[k] = -g_y_0_xxxy_xxy[k] * ab_x + g_y_0_xxxy_xxxy[k];

                g_y_0_xxxxy_xxz[k] = -g_y_0_xxxy_xxz[k] * ab_x + g_y_0_xxxy_xxxz[k];

                g_y_0_xxxxy_xyy[k] = -g_y_0_xxxy_xyy[k] * ab_x + g_y_0_xxxy_xxyy[k];

                g_y_0_xxxxy_xyz[k] = -g_y_0_xxxy_xyz[k] * ab_x + g_y_0_xxxy_xxyz[k];

                g_y_0_xxxxy_xzz[k] = -g_y_0_xxxy_xzz[k] * ab_x + g_y_0_xxxy_xxzz[k];

                g_y_0_xxxxy_yyy[k] = -g_y_0_xxxy_yyy[k] * ab_x + g_y_0_xxxy_xyyy[k];

                g_y_0_xxxxy_yyz[k] = -g_y_0_xxxy_yyz[k] * ab_x + g_y_0_xxxy_xyyz[k];

                g_y_0_xxxxy_yzz[k] = -g_y_0_xxxy_yzz[k] * ab_x + g_y_0_xxxy_xyzz[k];

                g_y_0_xxxxy_zzz[k] = -g_y_0_xxxy_zzz[k] * ab_x + g_y_0_xxxy_xzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxz_xxx, g_y_0_xxxxz_xxy, g_y_0_xxxxz_xxz, g_y_0_xxxxz_xyy, g_y_0_xxxxz_xyz, g_y_0_xxxxz_xzz, g_y_0_xxxxz_yyy, g_y_0_xxxxz_yyz, g_y_0_xxxxz_yzz, g_y_0_xxxxz_zzz, g_y_0_xxxz_xxx, g_y_0_xxxz_xxxx, g_y_0_xxxz_xxxy, g_y_0_xxxz_xxxz, g_y_0_xxxz_xxy, g_y_0_xxxz_xxyy, g_y_0_xxxz_xxyz, g_y_0_xxxz_xxz, g_y_0_xxxz_xxzz, g_y_0_xxxz_xyy, g_y_0_xxxz_xyyy, g_y_0_xxxz_xyyz, g_y_0_xxxz_xyz, g_y_0_xxxz_xyzz, g_y_0_xxxz_xzz, g_y_0_xxxz_xzzz, g_y_0_xxxz_yyy, g_y_0_xxxz_yyz, g_y_0_xxxz_yzz, g_y_0_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_xxx[k] = -g_y_0_xxxz_xxx[k] * ab_x + g_y_0_xxxz_xxxx[k];

                g_y_0_xxxxz_xxy[k] = -g_y_0_xxxz_xxy[k] * ab_x + g_y_0_xxxz_xxxy[k];

                g_y_0_xxxxz_xxz[k] = -g_y_0_xxxz_xxz[k] * ab_x + g_y_0_xxxz_xxxz[k];

                g_y_0_xxxxz_xyy[k] = -g_y_0_xxxz_xyy[k] * ab_x + g_y_0_xxxz_xxyy[k];

                g_y_0_xxxxz_xyz[k] = -g_y_0_xxxz_xyz[k] * ab_x + g_y_0_xxxz_xxyz[k];

                g_y_0_xxxxz_xzz[k] = -g_y_0_xxxz_xzz[k] * ab_x + g_y_0_xxxz_xxzz[k];

                g_y_0_xxxxz_yyy[k] = -g_y_0_xxxz_yyy[k] * ab_x + g_y_0_xxxz_xyyy[k];

                g_y_0_xxxxz_yyz[k] = -g_y_0_xxxz_yyz[k] * ab_x + g_y_0_xxxz_xyyz[k];

                g_y_0_xxxxz_yzz[k] = -g_y_0_xxxz_yzz[k] * ab_x + g_y_0_xxxz_xyzz[k];

                g_y_0_xxxxz_zzz[k] = -g_y_0_xxxz_zzz[k] * ab_x + g_y_0_xxxz_xzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyy_xxx, g_y_0_xxxyy_xxy, g_y_0_xxxyy_xxz, g_y_0_xxxyy_xyy, g_y_0_xxxyy_xyz, g_y_0_xxxyy_xzz, g_y_0_xxxyy_yyy, g_y_0_xxxyy_yyz, g_y_0_xxxyy_yzz, g_y_0_xxxyy_zzz, g_y_0_xxyy_xxx, g_y_0_xxyy_xxxx, g_y_0_xxyy_xxxy, g_y_0_xxyy_xxxz, g_y_0_xxyy_xxy, g_y_0_xxyy_xxyy, g_y_0_xxyy_xxyz, g_y_0_xxyy_xxz, g_y_0_xxyy_xxzz, g_y_0_xxyy_xyy, g_y_0_xxyy_xyyy, g_y_0_xxyy_xyyz, g_y_0_xxyy_xyz, g_y_0_xxyy_xyzz, g_y_0_xxyy_xzz, g_y_0_xxyy_xzzz, g_y_0_xxyy_yyy, g_y_0_xxyy_yyz, g_y_0_xxyy_yzz, g_y_0_xxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_xxx[k] = -g_y_0_xxyy_xxx[k] * ab_x + g_y_0_xxyy_xxxx[k];

                g_y_0_xxxyy_xxy[k] = -g_y_0_xxyy_xxy[k] * ab_x + g_y_0_xxyy_xxxy[k];

                g_y_0_xxxyy_xxz[k] = -g_y_0_xxyy_xxz[k] * ab_x + g_y_0_xxyy_xxxz[k];

                g_y_0_xxxyy_xyy[k] = -g_y_0_xxyy_xyy[k] * ab_x + g_y_0_xxyy_xxyy[k];

                g_y_0_xxxyy_xyz[k] = -g_y_0_xxyy_xyz[k] * ab_x + g_y_0_xxyy_xxyz[k];

                g_y_0_xxxyy_xzz[k] = -g_y_0_xxyy_xzz[k] * ab_x + g_y_0_xxyy_xxzz[k];

                g_y_0_xxxyy_yyy[k] = -g_y_0_xxyy_yyy[k] * ab_x + g_y_0_xxyy_xyyy[k];

                g_y_0_xxxyy_yyz[k] = -g_y_0_xxyy_yyz[k] * ab_x + g_y_0_xxyy_xyyz[k];

                g_y_0_xxxyy_yzz[k] = -g_y_0_xxyy_yzz[k] * ab_x + g_y_0_xxyy_xyzz[k];

                g_y_0_xxxyy_zzz[k] = -g_y_0_xxyy_zzz[k] * ab_x + g_y_0_xxyy_xzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 251 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 259 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyz_xxx, g_y_0_xxxyz_xxy, g_y_0_xxxyz_xxz, g_y_0_xxxyz_xyy, g_y_0_xxxyz_xyz, g_y_0_xxxyz_xzz, g_y_0_xxxyz_yyy, g_y_0_xxxyz_yyz, g_y_0_xxxyz_yzz, g_y_0_xxxyz_zzz, g_y_0_xxyz_xxx, g_y_0_xxyz_xxxx, g_y_0_xxyz_xxxy, g_y_0_xxyz_xxxz, g_y_0_xxyz_xxy, g_y_0_xxyz_xxyy, g_y_0_xxyz_xxyz, g_y_0_xxyz_xxz, g_y_0_xxyz_xxzz, g_y_0_xxyz_xyy, g_y_0_xxyz_xyyy, g_y_0_xxyz_xyyz, g_y_0_xxyz_xyz, g_y_0_xxyz_xyzz, g_y_0_xxyz_xzz, g_y_0_xxyz_xzzz, g_y_0_xxyz_yyy, g_y_0_xxyz_yyz, g_y_0_xxyz_yzz, g_y_0_xxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_xxx[k] = -g_y_0_xxyz_xxx[k] * ab_x + g_y_0_xxyz_xxxx[k];

                g_y_0_xxxyz_xxy[k] = -g_y_0_xxyz_xxy[k] * ab_x + g_y_0_xxyz_xxxy[k];

                g_y_0_xxxyz_xxz[k] = -g_y_0_xxyz_xxz[k] * ab_x + g_y_0_xxyz_xxxz[k];

                g_y_0_xxxyz_xyy[k] = -g_y_0_xxyz_xyy[k] * ab_x + g_y_0_xxyz_xxyy[k];

                g_y_0_xxxyz_xyz[k] = -g_y_0_xxyz_xyz[k] * ab_x + g_y_0_xxyz_xxyz[k];

                g_y_0_xxxyz_xzz[k] = -g_y_0_xxyz_xzz[k] * ab_x + g_y_0_xxyz_xxzz[k];

                g_y_0_xxxyz_yyy[k] = -g_y_0_xxyz_yyy[k] * ab_x + g_y_0_xxyz_xyyy[k];

                g_y_0_xxxyz_yyz[k] = -g_y_0_xxyz_yyz[k] * ab_x + g_y_0_xxyz_xyyz[k];

                g_y_0_xxxyz_yzz[k] = -g_y_0_xxyz_yzz[k] * ab_x + g_y_0_xxyz_xyzz[k];

                g_y_0_xxxyz_zzz[k] = -g_y_0_xxyz_zzz[k] * ab_x + g_y_0_xxyz_xzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxzz_xxx, g_y_0_xxxzz_xxy, g_y_0_xxxzz_xxz, g_y_0_xxxzz_xyy, g_y_0_xxxzz_xyz, g_y_0_xxxzz_xzz, g_y_0_xxxzz_yyy, g_y_0_xxxzz_yyz, g_y_0_xxxzz_yzz, g_y_0_xxxzz_zzz, g_y_0_xxzz_xxx, g_y_0_xxzz_xxxx, g_y_0_xxzz_xxxy, g_y_0_xxzz_xxxz, g_y_0_xxzz_xxy, g_y_0_xxzz_xxyy, g_y_0_xxzz_xxyz, g_y_0_xxzz_xxz, g_y_0_xxzz_xxzz, g_y_0_xxzz_xyy, g_y_0_xxzz_xyyy, g_y_0_xxzz_xyyz, g_y_0_xxzz_xyz, g_y_0_xxzz_xyzz, g_y_0_xxzz_xzz, g_y_0_xxzz_xzzz, g_y_0_xxzz_yyy, g_y_0_xxzz_yyz, g_y_0_xxzz_yzz, g_y_0_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_xxx[k] = -g_y_0_xxzz_xxx[k] * ab_x + g_y_0_xxzz_xxxx[k];

                g_y_0_xxxzz_xxy[k] = -g_y_0_xxzz_xxy[k] * ab_x + g_y_0_xxzz_xxxy[k];

                g_y_0_xxxzz_xxz[k] = -g_y_0_xxzz_xxz[k] * ab_x + g_y_0_xxzz_xxxz[k];

                g_y_0_xxxzz_xyy[k] = -g_y_0_xxzz_xyy[k] * ab_x + g_y_0_xxzz_xxyy[k];

                g_y_0_xxxzz_xyz[k] = -g_y_0_xxzz_xyz[k] * ab_x + g_y_0_xxzz_xxyz[k];

                g_y_0_xxxzz_xzz[k] = -g_y_0_xxzz_xzz[k] * ab_x + g_y_0_xxzz_xxzz[k];

                g_y_0_xxxzz_yyy[k] = -g_y_0_xxzz_yyy[k] * ab_x + g_y_0_xxzz_xyyy[k];

                g_y_0_xxxzz_yyz[k] = -g_y_0_xxzz_yyz[k] * ab_x + g_y_0_xxzz_xyyz[k];

                g_y_0_xxxzz_yzz[k] = -g_y_0_xxzz_yzz[k] * ab_x + g_y_0_xxzz_xyzz[k];

                g_y_0_xxxzz_zzz[k] = -g_y_0_xxzz_zzz[k] * ab_x + g_y_0_xxzz_xzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyy_xxx, g_y_0_xxyyy_xxy, g_y_0_xxyyy_xxz, g_y_0_xxyyy_xyy, g_y_0_xxyyy_xyz, g_y_0_xxyyy_xzz, g_y_0_xxyyy_yyy, g_y_0_xxyyy_yyz, g_y_0_xxyyy_yzz, g_y_0_xxyyy_zzz, g_y_0_xyyy_xxx, g_y_0_xyyy_xxxx, g_y_0_xyyy_xxxy, g_y_0_xyyy_xxxz, g_y_0_xyyy_xxy, g_y_0_xyyy_xxyy, g_y_0_xyyy_xxyz, g_y_0_xyyy_xxz, g_y_0_xyyy_xxzz, g_y_0_xyyy_xyy, g_y_0_xyyy_xyyy, g_y_0_xyyy_xyyz, g_y_0_xyyy_xyz, g_y_0_xyyy_xyzz, g_y_0_xyyy_xzz, g_y_0_xyyy_xzzz, g_y_0_xyyy_yyy, g_y_0_xyyy_yyz, g_y_0_xyyy_yzz, g_y_0_xyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_xxx[k] = -g_y_0_xyyy_xxx[k] * ab_x + g_y_0_xyyy_xxxx[k];

                g_y_0_xxyyy_xxy[k] = -g_y_0_xyyy_xxy[k] * ab_x + g_y_0_xyyy_xxxy[k];

                g_y_0_xxyyy_xxz[k] = -g_y_0_xyyy_xxz[k] * ab_x + g_y_0_xyyy_xxxz[k];

                g_y_0_xxyyy_xyy[k] = -g_y_0_xyyy_xyy[k] * ab_x + g_y_0_xyyy_xxyy[k];

                g_y_0_xxyyy_xyz[k] = -g_y_0_xyyy_xyz[k] * ab_x + g_y_0_xyyy_xxyz[k];

                g_y_0_xxyyy_xzz[k] = -g_y_0_xyyy_xzz[k] * ab_x + g_y_0_xyyy_xxzz[k];

                g_y_0_xxyyy_yyy[k] = -g_y_0_xyyy_yyy[k] * ab_x + g_y_0_xyyy_xyyy[k];

                g_y_0_xxyyy_yyz[k] = -g_y_0_xyyy_yyz[k] * ab_x + g_y_0_xyyy_xyyz[k];

                g_y_0_xxyyy_yzz[k] = -g_y_0_xyyy_yzz[k] * ab_x + g_y_0_xyyy_xyzz[k];

                g_y_0_xxyyy_zzz[k] = -g_y_0_xyyy_zzz[k] * ab_x + g_y_0_xyyy_xzzz[k];
            }

            /// Set up 280-290 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 289 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyz_xxx, g_y_0_xxyyz_xxy, g_y_0_xxyyz_xxz, g_y_0_xxyyz_xyy, g_y_0_xxyyz_xyz, g_y_0_xxyyz_xzz, g_y_0_xxyyz_yyy, g_y_0_xxyyz_yyz, g_y_0_xxyyz_yzz, g_y_0_xxyyz_zzz, g_y_0_xyyz_xxx, g_y_0_xyyz_xxxx, g_y_0_xyyz_xxxy, g_y_0_xyyz_xxxz, g_y_0_xyyz_xxy, g_y_0_xyyz_xxyy, g_y_0_xyyz_xxyz, g_y_0_xyyz_xxz, g_y_0_xyyz_xxzz, g_y_0_xyyz_xyy, g_y_0_xyyz_xyyy, g_y_0_xyyz_xyyz, g_y_0_xyyz_xyz, g_y_0_xyyz_xyzz, g_y_0_xyyz_xzz, g_y_0_xyyz_xzzz, g_y_0_xyyz_yyy, g_y_0_xyyz_yyz, g_y_0_xyyz_yzz, g_y_0_xyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_xxx[k] = -g_y_0_xyyz_xxx[k] * ab_x + g_y_0_xyyz_xxxx[k];

                g_y_0_xxyyz_xxy[k] = -g_y_0_xyyz_xxy[k] * ab_x + g_y_0_xyyz_xxxy[k];

                g_y_0_xxyyz_xxz[k] = -g_y_0_xyyz_xxz[k] * ab_x + g_y_0_xyyz_xxxz[k];

                g_y_0_xxyyz_xyy[k] = -g_y_0_xyyz_xyy[k] * ab_x + g_y_0_xyyz_xxyy[k];

                g_y_0_xxyyz_xyz[k] = -g_y_0_xyyz_xyz[k] * ab_x + g_y_0_xyyz_xxyz[k];

                g_y_0_xxyyz_xzz[k] = -g_y_0_xyyz_xzz[k] * ab_x + g_y_0_xyyz_xxzz[k];

                g_y_0_xxyyz_yyy[k] = -g_y_0_xyyz_yyy[k] * ab_x + g_y_0_xyyz_xyyy[k];

                g_y_0_xxyyz_yyz[k] = -g_y_0_xyyz_yyz[k] * ab_x + g_y_0_xyyz_xyyz[k];

                g_y_0_xxyyz_yzz[k] = -g_y_0_xyyz_yzz[k] * ab_x + g_y_0_xyyz_xyzz[k];

                g_y_0_xxyyz_zzz[k] = -g_y_0_xyyz_zzz[k] * ab_x + g_y_0_xyyz_xzzz[k];
            }

            /// Set up 290-300 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyzz_xxx, g_y_0_xxyzz_xxy, g_y_0_xxyzz_xxz, g_y_0_xxyzz_xyy, g_y_0_xxyzz_xyz, g_y_0_xxyzz_xzz, g_y_0_xxyzz_yyy, g_y_0_xxyzz_yyz, g_y_0_xxyzz_yzz, g_y_0_xxyzz_zzz, g_y_0_xyzz_xxx, g_y_0_xyzz_xxxx, g_y_0_xyzz_xxxy, g_y_0_xyzz_xxxz, g_y_0_xyzz_xxy, g_y_0_xyzz_xxyy, g_y_0_xyzz_xxyz, g_y_0_xyzz_xxz, g_y_0_xyzz_xxzz, g_y_0_xyzz_xyy, g_y_0_xyzz_xyyy, g_y_0_xyzz_xyyz, g_y_0_xyzz_xyz, g_y_0_xyzz_xyzz, g_y_0_xyzz_xzz, g_y_0_xyzz_xzzz, g_y_0_xyzz_yyy, g_y_0_xyzz_yyz, g_y_0_xyzz_yzz, g_y_0_xyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_xxx[k] = -g_y_0_xyzz_xxx[k] * ab_x + g_y_0_xyzz_xxxx[k];

                g_y_0_xxyzz_xxy[k] = -g_y_0_xyzz_xxy[k] * ab_x + g_y_0_xyzz_xxxy[k];

                g_y_0_xxyzz_xxz[k] = -g_y_0_xyzz_xxz[k] * ab_x + g_y_0_xyzz_xxxz[k];

                g_y_0_xxyzz_xyy[k] = -g_y_0_xyzz_xyy[k] * ab_x + g_y_0_xyzz_xxyy[k];

                g_y_0_xxyzz_xyz[k] = -g_y_0_xyzz_xyz[k] * ab_x + g_y_0_xyzz_xxyz[k];

                g_y_0_xxyzz_xzz[k] = -g_y_0_xyzz_xzz[k] * ab_x + g_y_0_xyzz_xxzz[k];

                g_y_0_xxyzz_yyy[k] = -g_y_0_xyzz_yyy[k] * ab_x + g_y_0_xyzz_xyyy[k];

                g_y_0_xxyzz_yyz[k] = -g_y_0_xyzz_yyz[k] * ab_x + g_y_0_xyzz_xyyz[k];

                g_y_0_xxyzz_yzz[k] = -g_y_0_xyzz_yzz[k] * ab_x + g_y_0_xyzz_xyzz[k];

                g_y_0_xxyzz_zzz[k] = -g_y_0_xyzz_zzz[k] * ab_x + g_y_0_xyzz_xzzz[k];
            }

            /// Set up 300-310 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 309 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzzz_xxx, g_y_0_xxzzz_xxy, g_y_0_xxzzz_xxz, g_y_0_xxzzz_xyy, g_y_0_xxzzz_xyz, g_y_0_xxzzz_xzz, g_y_0_xxzzz_yyy, g_y_0_xxzzz_yyz, g_y_0_xxzzz_yzz, g_y_0_xxzzz_zzz, g_y_0_xzzz_xxx, g_y_0_xzzz_xxxx, g_y_0_xzzz_xxxy, g_y_0_xzzz_xxxz, g_y_0_xzzz_xxy, g_y_0_xzzz_xxyy, g_y_0_xzzz_xxyz, g_y_0_xzzz_xxz, g_y_0_xzzz_xxzz, g_y_0_xzzz_xyy, g_y_0_xzzz_xyyy, g_y_0_xzzz_xyyz, g_y_0_xzzz_xyz, g_y_0_xzzz_xyzz, g_y_0_xzzz_xzz, g_y_0_xzzz_xzzz, g_y_0_xzzz_yyy, g_y_0_xzzz_yyz, g_y_0_xzzz_yzz, g_y_0_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_xxx[k] = -g_y_0_xzzz_xxx[k] * ab_x + g_y_0_xzzz_xxxx[k];

                g_y_0_xxzzz_xxy[k] = -g_y_0_xzzz_xxy[k] * ab_x + g_y_0_xzzz_xxxy[k];

                g_y_0_xxzzz_xxz[k] = -g_y_0_xzzz_xxz[k] * ab_x + g_y_0_xzzz_xxxz[k];

                g_y_0_xxzzz_xyy[k] = -g_y_0_xzzz_xyy[k] * ab_x + g_y_0_xzzz_xxyy[k];

                g_y_0_xxzzz_xyz[k] = -g_y_0_xzzz_xyz[k] * ab_x + g_y_0_xzzz_xxyz[k];

                g_y_0_xxzzz_xzz[k] = -g_y_0_xzzz_xzz[k] * ab_x + g_y_0_xzzz_xxzz[k];

                g_y_0_xxzzz_yyy[k] = -g_y_0_xzzz_yyy[k] * ab_x + g_y_0_xzzz_xyyy[k];

                g_y_0_xxzzz_yyz[k] = -g_y_0_xzzz_yyz[k] * ab_x + g_y_0_xzzz_xyyz[k];

                g_y_0_xxzzz_yzz[k] = -g_y_0_xzzz_yzz[k] * ab_x + g_y_0_xzzz_xyzz[k];

                g_y_0_xxzzz_zzz[k] = -g_y_0_xzzz_zzz[k] * ab_x + g_y_0_xzzz_xzzz[k];
            }

            /// Set up 310-320 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 319 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyy_xxx, g_y_0_xyyyy_xxy, g_y_0_xyyyy_xxz, g_y_0_xyyyy_xyy, g_y_0_xyyyy_xyz, g_y_0_xyyyy_xzz, g_y_0_xyyyy_yyy, g_y_0_xyyyy_yyz, g_y_0_xyyyy_yzz, g_y_0_xyyyy_zzz, g_y_0_yyyy_xxx, g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxy, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyy, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyy, g_y_0_yyyy_yyz, g_y_0_yyyy_yzz, g_y_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_xxx[k] = -g_y_0_yyyy_xxx[k] * ab_x + g_y_0_yyyy_xxxx[k];

                g_y_0_xyyyy_xxy[k] = -g_y_0_yyyy_xxy[k] * ab_x + g_y_0_yyyy_xxxy[k];

                g_y_0_xyyyy_xxz[k] = -g_y_0_yyyy_xxz[k] * ab_x + g_y_0_yyyy_xxxz[k];

                g_y_0_xyyyy_xyy[k] = -g_y_0_yyyy_xyy[k] * ab_x + g_y_0_yyyy_xxyy[k];

                g_y_0_xyyyy_xyz[k] = -g_y_0_yyyy_xyz[k] * ab_x + g_y_0_yyyy_xxyz[k];

                g_y_0_xyyyy_xzz[k] = -g_y_0_yyyy_xzz[k] * ab_x + g_y_0_yyyy_xxzz[k];

                g_y_0_xyyyy_yyy[k] = -g_y_0_yyyy_yyy[k] * ab_x + g_y_0_yyyy_xyyy[k];

                g_y_0_xyyyy_yyz[k] = -g_y_0_yyyy_yyz[k] * ab_x + g_y_0_yyyy_xyyz[k];

                g_y_0_xyyyy_yzz[k] = -g_y_0_yyyy_yzz[k] * ab_x + g_y_0_yyyy_xyzz[k];

                g_y_0_xyyyy_zzz[k] = -g_y_0_yyyy_zzz[k] * ab_x + g_y_0_yyyy_xzzz[k];
            }

            /// Set up 320-330 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyz_xxx, g_y_0_xyyyz_xxy, g_y_0_xyyyz_xxz, g_y_0_xyyyz_xyy, g_y_0_xyyyz_xyz, g_y_0_xyyyz_xzz, g_y_0_xyyyz_yyy, g_y_0_xyyyz_yyz, g_y_0_xyyyz_yzz, g_y_0_xyyyz_zzz, g_y_0_yyyz_xxx, g_y_0_yyyz_xxxx, g_y_0_yyyz_xxxy, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxy, g_y_0_yyyz_xxyy, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xyy, g_y_0_yyyz_xyyy, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_yyy, g_y_0_yyyz_yyz, g_y_0_yyyz_yzz, g_y_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_xxx[k] = -g_y_0_yyyz_xxx[k] * ab_x + g_y_0_yyyz_xxxx[k];

                g_y_0_xyyyz_xxy[k] = -g_y_0_yyyz_xxy[k] * ab_x + g_y_0_yyyz_xxxy[k];

                g_y_0_xyyyz_xxz[k] = -g_y_0_yyyz_xxz[k] * ab_x + g_y_0_yyyz_xxxz[k];

                g_y_0_xyyyz_xyy[k] = -g_y_0_yyyz_xyy[k] * ab_x + g_y_0_yyyz_xxyy[k];

                g_y_0_xyyyz_xyz[k] = -g_y_0_yyyz_xyz[k] * ab_x + g_y_0_yyyz_xxyz[k];

                g_y_0_xyyyz_xzz[k] = -g_y_0_yyyz_xzz[k] * ab_x + g_y_0_yyyz_xxzz[k];

                g_y_0_xyyyz_yyy[k] = -g_y_0_yyyz_yyy[k] * ab_x + g_y_0_yyyz_xyyy[k];

                g_y_0_xyyyz_yyz[k] = -g_y_0_yyyz_yyz[k] * ab_x + g_y_0_yyyz_xyyz[k];

                g_y_0_xyyyz_yzz[k] = -g_y_0_yyyz_yzz[k] * ab_x + g_y_0_yyyz_xyzz[k];

                g_y_0_xyyyz_zzz[k] = -g_y_0_yyyz_zzz[k] * ab_x + g_y_0_yyyz_xzzz[k];
            }

            /// Set up 330-340 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 339 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyzz_xxx, g_y_0_xyyzz_xxy, g_y_0_xyyzz_xxz, g_y_0_xyyzz_xyy, g_y_0_xyyzz_xyz, g_y_0_xyyzz_xzz, g_y_0_xyyzz_yyy, g_y_0_xyyzz_yyz, g_y_0_xyyzz_yzz, g_y_0_xyyzz_zzz, g_y_0_yyzz_xxx, g_y_0_yyzz_xxxx, g_y_0_yyzz_xxxy, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxy, g_y_0_yyzz_xxyy, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xyy, g_y_0_yyzz_xyyy, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_yyy, g_y_0_yyzz_yyz, g_y_0_yyzz_yzz, g_y_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_xxx[k] = -g_y_0_yyzz_xxx[k] * ab_x + g_y_0_yyzz_xxxx[k];

                g_y_0_xyyzz_xxy[k] = -g_y_0_yyzz_xxy[k] * ab_x + g_y_0_yyzz_xxxy[k];

                g_y_0_xyyzz_xxz[k] = -g_y_0_yyzz_xxz[k] * ab_x + g_y_0_yyzz_xxxz[k];

                g_y_0_xyyzz_xyy[k] = -g_y_0_yyzz_xyy[k] * ab_x + g_y_0_yyzz_xxyy[k];

                g_y_0_xyyzz_xyz[k] = -g_y_0_yyzz_xyz[k] * ab_x + g_y_0_yyzz_xxyz[k];

                g_y_0_xyyzz_xzz[k] = -g_y_0_yyzz_xzz[k] * ab_x + g_y_0_yyzz_xxzz[k];

                g_y_0_xyyzz_yyy[k] = -g_y_0_yyzz_yyy[k] * ab_x + g_y_0_yyzz_xyyy[k];

                g_y_0_xyyzz_yyz[k] = -g_y_0_yyzz_yyz[k] * ab_x + g_y_0_yyzz_xyyz[k];

                g_y_0_xyyzz_yzz[k] = -g_y_0_yyzz_yzz[k] * ab_x + g_y_0_yyzz_xyzz[k];

                g_y_0_xyyzz_zzz[k] = -g_y_0_yyzz_zzz[k] * ab_x + g_y_0_yyzz_xzzz[k];
            }

            /// Set up 340-350 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 349 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzzz_xxx, g_y_0_xyzzz_xxy, g_y_0_xyzzz_xxz, g_y_0_xyzzz_xyy, g_y_0_xyzzz_xyz, g_y_0_xyzzz_xzz, g_y_0_xyzzz_yyy, g_y_0_xyzzz_yyz, g_y_0_xyzzz_yzz, g_y_0_xyzzz_zzz, g_y_0_yzzz_xxx, g_y_0_yzzz_xxxx, g_y_0_yzzz_xxxy, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxy, g_y_0_yzzz_xxyy, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xyy, g_y_0_yzzz_xyyy, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_yyy, g_y_0_yzzz_yyz, g_y_0_yzzz_yzz, g_y_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_xxx[k] = -g_y_0_yzzz_xxx[k] * ab_x + g_y_0_yzzz_xxxx[k];

                g_y_0_xyzzz_xxy[k] = -g_y_0_yzzz_xxy[k] * ab_x + g_y_0_yzzz_xxxy[k];

                g_y_0_xyzzz_xxz[k] = -g_y_0_yzzz_xxz[k] * ab_x + g_y_0_yzzz_xxxz[k];

                g_y_0_xyzzz_xyy[k] = -g_y_0_yzzz_xyy[k] * ab_x + g_y_0_yzzz_xxyy[k];

                g_y_0_xyzzz_xyz[k] = -g_y_0_yzzz_xyz[k] * ab_x + g_y_0_yzzz_xxyz[k];

                g_y_0_xyzzz_xzz[k] = -g_y_0_yzzz_xzz[k] * ab_x + g_y_0_yzzz_xxzz[k];

                g_y_0_xyzzz_yyy[k] = -g_y_0_yzzz_yyy[k] * ab_x + g_y_0_yzzz_xyyy[k];

                g_y_0_xyzzz_yyz[k] = -g_y_0_yzzz_yyz[k] * ab_x + g_y_0_yzzz_xyyz[k];

                g_y_0_xyzzz_yzz[k] = -g_y_0_yzzz_yzz[k] * ab_x + g_y_0_yzzz_xyzz[k];

                g_y_0_xyzzz_zzz[k] = -g_y_0_yzzz_zzz[k] * ab_x + g_y_0_yzzz_xzzz[k];
            }

            /// Set up 350-360 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzzz_xxx, g_y_0_xzzzz_xxy, g_y_0_xzzzz_xxz, g_y_0_xzzzz_xyy, g_y_0_xzzzz_xyz, g_y_0_xzzzz_xzz, g_y_0_xzzzz_yyy, g_y_0_xzzzz_yyz, g_y_0_xzzzz_yzz, g_y_0_xzzzz_zzz, g_y_0_zzzz_xxx, g_y_0_zzzz_xxxx, g_y_0_zzzz_xxxy, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxy, g_y_0_zzzz_xxyy, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xyy, g_y_0_zzzz_xyyy, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_yyy, g_y_0_zzzz_yyz, g_y_0_zzzz_yzz, g_y_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_xxx[k] = -g_y_0_zzzz_xxx[k] * ab_x + g_y_0_zzzz_xxxx[k];

                g_y_0_xzzzz_xxy[k] = -g_y_0_zzzz_xxy[k] * ab_x + g_y_0_zzzz_xxxy[k];

                g_y_0_xzzzz_xxz[k] = -g_y_0_zzzz_xxz[k] * ab_x + g_y_0_zzzz_xxxz[k];

                g_y_0_xzzzz_xyy[k] = -g_y_0_zzzz_xyy[k] * ab_x + g_y_0_zzzz_xxyy[k];

                g_y_0_xzzzz_xyz[k] = -g_y_0_zzzz_xyz[k] * ab_x + g_y_0_zzzz_xxyz[k];

                g_y_0_xzzzz_xzz[k] = -g_y_0_zzzz_xzz[k] * ab_x + g_y_0_zzzz_xxzz[k];

                g_y_0_xzzzz_yyy[k] = -g_y_0_zzzz_yyy[k] * ab_x + g_y_0_zzzz_xyyy[k];

                g_y_0_xzzzz_yyz[k] = -g_y_0_zzzz_yyz[k] * ab_x + g_y_0_zzzz_xyyz[k];

                g_y_0_xzzzz_yzz[k] = -g_y_0_zzzz_yzz[k] * ab_x + g_y_0_zzzz_xyzz[k];

                g_y_0_xzzzz_zzz[k] = -g_y_0_zzzz_zzz[k] * ab_x + g_y_0_zzzz_xzzz[k];
            }

            /// Set up 360-370 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 369 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyy_xxx, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxy, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxz, g_y_0_yyyy_xyy, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzz, g_y_0_yyyy_yyy, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzz, g_y_0_yyyyy_xxx, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_yyy, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_zzz, g_yyyy_xxx, g_yyyy_xxy, g_yyyy_xxz, g_yyyy_xyy, g_yyyy_xyz, g_yyyy_xzz, g_yyyy_yyy, g_yyyy_yyz, g_yyyy_yzz, g_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_xxx[k] = -g_yyyy_xxx[k] - g_y_0_yyyy_xxx[k] * ab_y + g_y_0_yyyy_xxxy[k];

                g_y_0_yyyyy_xxy[k] = -g_yyyy_xxy[k] - g_y_0_yyyy_xxy[k] * ab_y + g_y_0_yyyy_xxyy[k];

                g_y_0_yyyyy_xxz[k] = -g_yyyy_xxz[k] - g_y_0_yyyy_xxz[k] * ab_y + g_y_0_yyyy_xxyz[k];

                g_y_0_yyyyy_xyy[k] = -g_yyyy_xyy[k] - g_y_0_yyyy_xyy[k] * ab_y + g_y_0_yyyy_xyyy[k];

                g_y_0_yyyyy_xyz[k] = -g_yyyy_xyz[k] - g_y_0_yyyy_xyz[k] * ab_y + g_y_0_yyyy_xyyz[k];

                g_y_0_yyyyy_xzz[k] = -g_yyyy_xzz[k] - g_y_0_yyyy_xzz[k] * ab_y + g_y_0_yyyy_xyzz[k];

                g_y_0_yyyyy_yyy[k] = -g_yyyy_yyy[k] - g_y_0_yyyy_yyy[k] * ab_y + g_y_0_yyyy_yyyy[k];

                g_y_0_yyyyy_yyz[k] = -g_yyyy_yyz[k] - g_y_0_yyyy_yyz[k] * ab_y + g_y_0_yyyy_yyyz[k];

                g_y_0_yyyyy_yzz[k] = -g_yyyy_yzz[k] - g_y_0_yyyy_yzz[k] * ab_y + g_y_0_yyyy_yyzz[k];

                g_y_0_yyyyy_zzz[k] = -g_yyyy_zzz[k] - g_y_0_yyyy_zzz[k] * ab_y + g_y_0_yyyy_yzzz[k];
            }

            /// Set up 370-380 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 379 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyy_xxx, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzz, g_y_0_yyyy_zzzz, g_y_0_yyyyz_xxx, g_y_0_yyyyz_xxy, g_y_0_yyyyz_xxz, g_y_0_yyyyz_xyy, g_y_0_yyyyz_xyz, g_y_0_yyyyz_xzz, g_y_0_yyyyz_yyy, g_y_0_yyyyz_yyz, g_y_0_yyyyz_yzz, g_y_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_xxx[k] = -g_y_0_yyyy_xxx[k] * ab_z + g_y_0_yyyy_xxxz[k];

                g_y_0_yyyyz_xxy[k] = -g_y_0_yyyy_xxy[k] * ab_z + g_y_0_yyyy_xxyz[k];

                g_y_0_yyyyz_xxz[k] = -g_y_0_yyyy_xxz[k] * ab_z + g_y_0_yyyy_xxzz[k];

                g_y_0_yyyyz_xyy[k] = -g_y_0_yyyy_xyy[k] * ab_z + g_y_0_yyyy_xyyz[k];

                g_y_0_yyyyz_xyz[k] = -g_y_0_yyyy_xyz[k] * ab_z + g_y_0_yyyy_xyzz[k];

                g_y_0_yyyyz_xzz[k] = -g_y_0_yyyy_xzz[k] * ab_z + g_y_0_yyyy_xzzz[k];

                g_y_0_yyyyz_yyy[k] = -g_y_0_yyyy_yyy[k] * ab_z + g_y_0_yyyy_yyyz[k];

                g_y_0_yyyyz_yyz[k] = -g_y_0_yyyy_yyz[k] * ab_z + g_y_0_yyyy_yyzz[k];

                g_y_0_yyyyz_yzz[k] = -g_y_0_yyyy_yzz[k] * ab_z + g_y_0_yyyy_yzzz[k];

                g_y_0_yyyyz_zzz[k] = -g_y_0_yyyy_zzz[k] * ab_z + g_y_0_yyyy_zzzz[k];
            }

            /// Set up 380-390 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyz_xxx, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxy, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xyy, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_yyy, g_y_0_yyyz_yyyz, g_y_0_yyyz_yyz, g_y_0_yyyz_yyzz, g_y_0_yyyz_yzz, g_y_0_yyyz_yzzz, g_y_0_yyyz_zzz, g_y_0_yyyz_zzzz, g_y_0_yyyzz_xxx, g_y_0_yyyzz_xxy, g_y_0_yyyzz_xxz, g_y_0_yyyzz_xyy, g_y_0_yyyzz_xyz, g_y_0_yyyzz_xzz, g_y_0_yyyzz_yyy, g_y_0_yyyzz_yyz, g_y_0_yyyzz_yzz, g_y_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_xxx[k] = -g_y_0_yyyz_xxx[k] * ab_z + g_y_0_yyyz_xxxz[k];

                g_y_0_yyyzz_xxy[k] = -g_y_0_yyyz_xxy[k] * ab_z + g_y_0_yyyz_xxyz[k];

                g_y_0_yyyzz_xxz[k] = -g_y_0_yyyz_xxz[k] * ab_z + g_y_0_yyyz_xxzz[k];

                g_y_0_yyyzz_xyy[k] = -g_y_0_yyyz_xyy[k] * ab_z + g_y_0_yyyz_xyyz[k];

                g_y_0_yyyzz_xyz[k] = -g_y_0_yyyz_xyz[k] * ab_z + g_y_0_yyyz_xyzz[k];

                g_y_0_yyyzz_xzz[k] = -g_y_0_yyyz_xzz[k] * ab_z + g_y_0_yyyz_xzzz[k];

                g_y_0_yyyzz_yyy[k] = -g_y_0_yyyz_yyy[k] * ab_z + g_y_0_yyyz_yyyz[k];

                g_y_0_yyyzz_yyz[k] = -g_y_0_yyyz_yyz[k] * ab_z + g_y_0_yyyz_yyzz[k];

                g_y_0_yyyzz_yzz[k] = -g_y_0_yyyz_yzz[k] * ab_z + g_y_0_yyyz_yzzz[k];

                g_y_0_yyyzz_zzz[k] = -g_y_0_yyyz_zzz[k] * ab_z + g_y_0_yyyz_zzzz[k];
            }

            /// Set up 390-400 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 399 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyzz_xxx, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxy, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xyy, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_yyy, g_y_0_yyzz_yyyz, g_y_0_yyzz_yyz, g_y_0_yyzz_yyzz, g_y_0_yyzz_yzz, g_y_0_yyzz_yzzz, g_y_0_yyzz_zzz, g_y_0_yyzz_zzzz, g_y_0_yyzzz_xxx, g_y_0_yyzzz_xxy, g_y_0_yyzzz_xxz, g_y_0_yyzzz_xyy, g_y_0_yyzzz_xyz, g_y_0_yyzzz_xzz, g_y_0_yyzzz_yyy, g_y_0_yyzzz_yyz, g_y_0_yyzzz_yzz, g_y_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_xxx[k] = -g_y_0_yyzz_xxx[k] * ab_z + g_y_0_yyzz_xxxz[k];

                g_y_0_yyzzz_xxy[k] = -g_y_0_yyzz_xxy[k] * ab_z + g_y_0_yyzz_xxyz[k];

                g_y_0_yyzzz_xxz[k] = -g_y_0_yyzz_xxz[k] * ab_z + g_y_0_yyzz_xxzz[k];

                g_y_0_yyzzz_xyy[k] = -g_y_0_yyzz_xyy[k] * ab_z + g_y_0_yyzz_xyyz[k];

                g_y_0_yyzzz_xyz[k] = -g_y_0_yyzz_xyz[k] * ab_z + g_y_0_yyzz_xyzz[k];

                g_y_0_yyzzz_xzz[k] = -g_y_0_yyzz_xzz[k] * ab_z + g_y_0_yyzz_xzzz[k];

                g_y_0_yyzzz_yyy[k] = -g_y_0_yyzz_yyy[k] * ab_z + g_y_0_yyzz_yyyz[k];

                g_y_0_yyzzz_yyz[k] = -g_y_0_yyzz_yyz[k] * ab_z + g_y_0_yyzz_yyzz[k];

                g_y_0_yyzzz_yzz[k] = -g_y_0_yyzz_yzz[k] * ab_z + g_y_0_yyzz_yzzz[k];

                g_y_0_yyzzz_zzz[k] = -g_y_0_yyzz_zzz[k] * ab_z + g_y_0_yyzz_zzzz[k];
            }

            /// Set up 400-410 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 409 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzzz_xxx, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxy, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xyy, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_yyy, g_y_0_yzzz_yyyz, g_y_0_yzzz_yyz, g_y_0_yzzz_yyzz, g_y_0_yzzz_yzz, g_y_0_yzzz_yzzz, g_y_0_yzzz_zzz, g_y_0_yzzz_zzzz, g_y_0_yzzzz_xxx, g_y_0_yzzzz_xxy, g_y_0_yzzzz_xxz, g_y_0_yzzzz_xyy, g_y_0_yzzzz_xyz, g_y_0_yzzzz_xzz, g_y_0_yzzzz_yyy, g_y_0_yzzzz_yyz, g_y_0_yzzzz_yzz, g_y_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_xxx[k] = -g_y_0_yzzz_xxx[k] * ab_z + g_y_0_yzzz_xxxz[k];

                g_y_0_yzzzz_xxy[k] = -g_y_0_yzzz_xxy[k] * ab_z + g_y_0_yzzz_xxyz[k];

                g_y_0_yzzzz_xxz[k] = -g_y_0_yzzz_xxz[k] * ab_z + g_y_0_yzzz_xxzz[k];

                g_y_0_yzzzz_xyy[k] = -g_y_0_yzzz_xyy[k] * ab_z + g_y_0_yzzz_xyyz[k];

                g_y_0_yzzzz_xyz[k] = -g_y_0_yzzz_xyz[k] * ab_z + g_y_0_yzzz_xyzz[k];

                g_y_0_yzzzz_xzz[k] = -g_y_0_yzzz_xzz[k] * ab_z + g_y_0_yzzz_xzzz[k];

                g_y_0_yzzzz_yyy[k] = -g_y_0_yzzz_yyy[k] * ab_z + g_y_0_yzzz_yyyz[k];

                g_y_0_yzzzz_yyz[k] = -g_y_0_yzzz_yyz[k] * ab_z + g_y_0_yzzz_yyzz[k];

                g_y_0_yzzzz_yzz[k] = -g_y_0_yzzz_yzz[k] * ab_z + g_y_0_yzzz_yzzz[k];

                g_y_0_yzzzz_zzz[k] = -g_y_0_yzzz_zzz[k] * ab_z + g_y_0_yzzz_zzzz[k];
            }

            /// Set up 410-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzz_xxx, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxy, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xyy, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_yyy, g_y_0_zzzz_yyyz, g_y_0_zzzz_yyz, g_y_0_zzzz_yyzz, g_y_0_zzzz_yzz, g_y_0_zzzz_yzzz, g_y_0_zzzz_zzz, g_y_0_zzzz_zzzz, g_y_0_zzzzz_xxx, g_y_0_zzzzz_xxy, g_y_0_zzzzz_xxz, g_y_0_zzzzz_xyy, g_y_0_zzzzz_xyz, g_y_0_zzzzz_xzz, g_y_0_zzzzz_yyy, g_y_0_zzzzz_yyz, g_y_0_zzzzz_yzz, g_y_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_xxx[k] = -g_y_0_zzzz_xxx[k] * ab_z + g_y_0_zzzz_xxxz[k];

                g_y_0_zzzzz_xxy[k] = -g_y_0_zzzz_xxy[k] * ab_z + g_y_0_zzzz_xxyz[k];

                g_y_0_zzzzz_xxz[k] = -g_y_0_zzzz_xxz[k] * ab_z + g_y_0_zzzz_xxzz[k];

                g_y_0_zzzzz_xyy[k] = -g_y_0_zzzz_xyy[k] * ab_z + g_y_0_zzzz_xyyz[k];

                g_y_0_zzzzz_xyz[k] = -g_y_0_zzzz_xyz[k] * ab_z + g_y_0_zzzz_xyzz[k];

                g_y_0_zzzzz_xzz[k] = -g_y_0_zzzz_xzz[k] * ab_z + g_y_0_zzzz_xzzz[k];

                g_y_0_zzzzz_yyy[k] = -g_y_0_zzzz_yyy[k] * ab_z + g_y_0_zzzz_yyyz[k];

                g_y_0_zzzzz_yyz[k] = -g_y_0_zzzz_yyz[k] * ab_z + g_y_0_zzzz_yyzz[k];

                g_y_0_zzzzz_yzz[k] = -g_y_0_zzzz_yzz[k] * ab_z + g_y_0_zzzz_yzzz[k];

                g_y_0_zzzzz_zzz[k] = -g_y_0_zzzz_zzz[k] * ab_z + g_y_0_zzzz_zzzz[k];
            }

            /// Set up 420-430 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 420 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 421 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 422 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 423 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 424 * ccomps * dcomps);

            auto g_z_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 425 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 426 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 427 * ccomps * dcomps);

            auto g_z_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 428 * ccomps * dcomps);

            auto g_z_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 429 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxx_xxx, g_z_0_xxxx_xxxx, g_z_0_xxxx_xxxy, g_z_0_xxxx_xxxz, g_z_0_xxxx_xxy, g_z_0_xxxx_xxyy, g_z_0_xxxx_xxyz, g_z_0_xxxx_xxz, g_z_0_xxxx_xxzz, g_z_0_xxxx_xyy, g_z_0_xxxx_xyyy, g_z_0_xxxx_xyyz, g_z_0_xxxx_xyz, g_z_0_xxxx_xyzz, g_z_0_xxxx_xzz, g_z_0_xxxx_xzzz, g_z_0_xxxx_yyy, g_z_0_xxxx_yyz, g_z_0_xxxx_yzz, g_z_0_xxxx_zzz, g_z_0_xxxxx_xxx, g_z_0_xxxxx_xxy, g_z_0_xxxxx_xxz, g_z_0_xxxxx_xyy, g_z_0_xxxxx_xyz, g_z_0_xxxxx_xzz, g_z_0_xxxxx_yyy, g_z_0_xxxxx_yyz, g_z_0_xxxxx_yzz, g_z_0_xxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_xxx[k] = -g_z_0_xxxx_xxx[k] * ab_x + g_z_0_xxxx_xxxx[k];

                g_z_0_xxxxx_xxy[k] = -g_z_0_xxxx_xxy[k] * ab_x + g_z_0_xxxx_xxxy[k];

                g_z_0_xxxxx_xxz[k] = -g_z_0_xxxx_xxz[k] * ab_x + g_z_0_xxxx_xxxz[k];

                g_z_0_xxxxx_xyy[k] = -g_z_0_xxxx_xyy[k] * ab_x + g_z_0_xxxx_xxyy[k];

                g_z_0_xxxxx_xyz[k] = -g_z_0_xxxx_xyz[k] * ab_x + g_z_0_xxxx_xxyz[k];

                g_z_0_xxxxx_xzz[k] = -g_z_0_xxxx_xzz[k] * ab_x + g_z_0_xxxx_xxzz[k];

                g_z_0_xxxxx_yyy[k] = -g_z_0_xxxx_yyy[k] * ab_x + g_z_0_xxxx_xyyy[k];

                g_z_0_xxxxx_yyz[k] = -g_z_0_xxxx_yyz[k] * ab_x + g_z_0_xxxx_xyyz[k];

                g_z_0_xxxxx_yzz[k] = -g_z_0_xxxx_yzz[k] * ab_x + g_z_0_xxxx_xyzz[k];

                g_z_0_xxxxx_zzz[k] = -g_z_0_xxxx_zzz[k] * ab_x + g_z_0_xxxx_xzzz[k];
            }

            /// Set up 430-440 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 430 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 431 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 432 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 433 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 434 * ccomps * dcomps);

            auto g_z_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 435 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 436 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 437 * ccomps * dcomps);

            auto g_z_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 438 * ccomps * dcomps);

            auto g_z_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 439 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxy_xxx, g_z_0_xxxxy_xxy, g_z_0_xxxxy_xxz, g_z_0_xxxxy_xyy, g_z_0_xxxxy_xyz, g_z_0_xxxxy_xzz, g_z_0_xxxxy_yyy, g_z_0_xxxxy_yyz, g_z_0_xxxxy_yzz, g_z_0_xxxxy_zzz, g_z_0_xxxy_xxx, g_z_0_xxxy_xxxx, g_z_0_xxxy_xxxy, g_z_0_xxxy_xxxz, g_z_0_xxxy_xxy, g_z_0_xxxy_xxyy, g_z_0_xxxy_xxyz, g_z_0_xxxy_xxz, g_z_0_xxxy_xxzz, g_z_0_xxxy_xyy, g_z_0_xxxy_xyyy, g_z_0_xxxy_xyyz, g_z_0_xxxy_xyz, g_z_0_xxxy_xyzz, g_z_0_xxxy_xzz, g_z_0_xxxy_xzzz, g_z_0_xxxy_yyy, g_z_0_xxxy_yyz, g_z_0_xxxy_yzz, g_z_0_xxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_xxx[k] = -g_z_0_xxxy_xxx[k] * ab_x + g_z_0_xxxy_xxxx[k];

                g_z_0_xxxxy_xxy[k] = -g_z_0_xxxy_xxy[k] * ab_x + g_z_0_xxxy_xxxy[k];

                g_z_0_xxxxy_xxz[k] = -g_z_0_xxxy_xxz[k] * ab_x + g_z_0_xxxy_xxxz[k];

                g_z_0_xxxxy_xyy[k] = -g_z_0_xxxy_xyy[k] * ab_x + g_z_0_xxxy_xxyy[k];

                g_z_0_xxxxy_xyz[k] = -g_z_0_xxxy_xyz[k] * ab_x + g_z_0_xxxy_xxyz[k];

                g_z_0_xxxxy_xzz[k] = -g_z_0_xxxy_xzz[k] * ab_x + g_z_0_xxxy_xxzz[k];

                g_z_0_xxxxy_yyy[k] = -g_z_0_xxxy_yyy[k] * ab_x + g_z_0_xxxy_xyyy[k];

                g_z_0_xxxxy_yyz[k] = -g_z_0_xxxy_yyz[k] * ab_x + g_z_0_xxxy_xyyz[k];

                g_z_0_xxxxy_yzz[k] = -g_z_0_xxxy_yzz[k] * ab_x + g_z_0_xxxy_xyzz[k];

                g_z_0_xxxxy_zzz[k] = -g_z_0_xxxy_zzz[k] * ab_x + g_z_0_xxxy_xzzz[k];
            }

            /// Set up 440-450 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 440 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 441 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 442 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 443 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 444 * ccomps * dcomps);

            auto g_z_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 445 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 446 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 447 * ccomps * dcomps);

            auto g_z_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 448 * ccomps * dcomps);

            auto g_z_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxz_xxx, g_z_0_xxxxz_xxy, g_z_0_xxxxz_xxz, g_z_0_xxxxz_xyy, g_z_0_xxxxz_xyz, g_z_0_xxxxz_xzz, g_z_0_xxxxz_yyy, g_z_0_xxxxz_yyz, g_z_0_xxxxz_yzz, g_z_0_xxxxz_zzz, g_z_0_xxxz_xxx, g_z_0_xxxz_xxxx, g_z_0_xxxz_xxxy, g_z_0_xxxz_xxxz, g_z_0_xxxz_xxy, g_z_0_xxxz_xxyy, g_z_0_xxxz_xxyz, g_z_0_xxxz_xxz, g_z_0_xxxz_xxzz, g_z_0_xxxz_xyy, g_z_0_xxxz_xyyy, g_z_0_xxxz_xyyz, g_z_0_xxxz_xyz, g_z_0_xxxz_xyzz, g_z_0_xxxz_xzz, g_z_0_xxxz_xzzz, g_z_0_xxxz_yyy, g_z_0_xxxz_yyz, g_z_0_xxxz_yzz, g_z_0_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_xxx[k] = -g_z_0_xxxz_xxx[k] * ab_x + g_z_0_xxxz_xxxx[k];

                g_z_0_xxxxz_xxy[k] = -g_z_0_xxxz_xxy[k] * ab_x + g_z_0_xxxz_xxxy[k];

                g_z_0_xxxxz_xxz[k] = -g_z_0_xxxz_xxz[k] * ab_x + g_z_0_xxxz_xxxz[k];

                g_z_0_xxxxz_xyy[k] = -g_z_0_xxxz_xyy[k] * ab_x + g_z_0_xxxz_xxyy[k];

                g_z_0_xxxxz_xyz[k] = -g_z_0_xxxz_xyz[k] * ab_x + g_z_0_xxxz_xxyz[k];

                g_z_0_xxxxz_xzz[k] = -g_z_0_xxxz_xzz[k] * ab_x + g_z_0_xxxz_xxzz[k];

                g_z_0_xxxxz_yyy[k] = -g_z_0_xxxz_yyy[k] * ab_x + g_z_0_xxxz_xyyy[k];

                g_z_0_xxxxz_yyz[k] = -g_z_0_xxxz_yyz[k] * ab_x + g_z_0_xxxz_xyyz[k];

                g_z_0_xxxxz_yzz[k] = -g_z_0_xxxz_yzz[k] * ab_x + g_z_0_xxxz_xyzz[k];

                g_z_0_xxxxz_zzz[k] = -g_z_0_xxxz_zzz[k] * ab_x + g_z_0_xxxz_xzzz[k];
            }

            /// Set up 450-460 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 450 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 451 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 452 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 453 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 454 * ccomps * dcomps);

            auto g_z_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 455 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 456 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 457 * ccomps * dcomps);

            auto g_z_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 458 * ccomps * dcomps);

            auto g_z_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 459 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyy_xxx, g_z_0_xxxyy_xxy, g_z_0_xxxyy_xxz, g_z_0_xxxyy_xyy, g_z_0_xxxyy_xyz, g_z_0_xxxyy_xzz, g_z_0_xxxyy_yyy, g_z_0_xxxyy_yyz, g_z_0_xxxyy_yzz, g_z_0_xxxyy_zzz, g_z_0_xxyy_xxx, g_z_0_xxyy_xxxx, g_z_0_xxyy_xxxy, g_z_0_xxyy_xxxz, g_z_0_xxyy_xxy, g_z_0_xxyy_xxyy, g_z_0_xxyy_xxyz, g_z_0_xxyy_xxz, g_z_0_xxyy_xxzz, g_z_0_xxyy_xyy, g_z_0_xxyy_xyyy, g_z_0_xxyy_xyyz, g_z_0_xxyy_xyz, g_z_0_xxyy_xyzz, g_z_0_xxyy_xzz, g_z_0_xxyy_xzzz, g_z_0_xxyy_yyy, g_z_0_xxyy_yyz, g_z_0_xxyy_yzz, g_z_0_xxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_xxx[k] = -g_z_0_xxyy_xxx[k] * ab_x + g_z_0_xxyy_xxxx[k];

                g_z_0_xxxyy_xxy[k] = -g_z_0_xxyy_xxy[k] * ab_x + g_z_0_xxyy_xxxy[k];

                g_z_0_xxxyy_xxz[k] = -g_z_0_xxyy_xxz[k] * ab_x + g_z_0_xxyy_xxxz[k];

                g_z_0_xxxyy_xyy[k] = -g_z_0_xxyy_xyy[k] * ab_x + g_z_0_xxyy_xxyy[k];

                g_z_0_xxxyy_xyz[k] = -g_z_0_xxyy_xyz[k] * ab_x + g_z_0_xxyy_xxyz[k];

                g_z_0_xxxyy_xzz[k] = -g_z_0_xxyy_xzz[k] * ab_x + g_z_0_xxyy_xxzz[k];

                g_z_0_xxxyy_yyy[k] = -g_z_0_xxyy_yyy[k] * ab_x + g_z_0_xxyy_xyyy[k];

                g_z_0_xxxyy_yyz[k] = -g_z_0_xxyy_yyz[k] * ab_x + g_z_0_xxyy_xyyz[k];

                g_z_0_xxxyy_yzz[k] = -g_z_0_xxyy_yzz[k] * ab_x + g_z_0_xxyy_xyzz[k];

                g_z_0_xxxyy_zzz[k] = -g_z_0_xxyy_zzz[k] * ab_x + g_z_0_xxyy_xzzz[k];
            }

            /// Set up 460-470 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 460 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 461 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 462 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 463 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 464 * ccomps * dcomps);

            auto g_z_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 465 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 466 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 467 * ccomps * dcomps);

            auto g_z_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 468 * ccomps * dcomps);

            auto g_z_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 469 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyz_xxx, g_z_0_xxxyz_xxy, g_z_0_xxxyz_xxz, g_z_0_xxxyz_xyy, g_z_0_xxxyz_xyz, g_z_0_xxxyz_xzz, g_z_0_xxxyz_yyy, g_z_0_xxxyz_yyz, g_z_0_xxxyz_yzz, g_z_0_xxxyz_zzz, g_z_0_xxyz_xxx, g_z_0_xxyz_xxxx, g_z_0_xxyz_xxxy, g_z_0_xxyz_xxxz, g_z_0_xxyz_xxy, g_z_0_xxyz_xxyy, g_z_0_xxyz_xxyz, g_z_0_xxyz_xxz, g_z_0_xxyz_xxzz, g_z_0_xxyz_xyy, g_z_0_xxyz_xyyy, g_z_0_xxyz_xyyz, g_z_0_xxyz_xyz, g_z_0_xxyz_xyzz, g_z_0_xxyz_xzz, g_z_0_xxyz_xzzz, g_z_0_xxyz_yyy, g_z_0_xxyz_yyz, g_z_0_xxyz_yzz, g_z_0_xxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_xxx[k] = -g_z_0_xxyz_xxx[k] * ab_x + g_z_0_xxyz_xxxx[k];

                g_z_0_xxxyz_xxy[k] = -g_z_0_xxyz_xxy[k] * ab_x + g_z_0_xxyz_xxxy[k];

                g_z_0_xxxyz_xxz[k] = -g_z_0_xxyz_xxz[k] * ab_x + g_z_0_xxyz_xxxz[k];

                g_z_0_xxxyz_xyy[k] = -g_z_0_xxyz_xyy[k] * ab_x + g_z_0_xxyz_xxyy[k];

                g_z_0_xxxyz_xyz[k] = -g_z_0_xxyz_xyz[k] * ab_x + g_z_0_xxyz_xxyz[k];

                g_z_0_xxxyz_xzz[k] = -g_z_0_xxyz_xzz[k] * ab_x + g_z_0_xxyz_xxzz[k];

                g_z_0_xxxyz_yyy[k] = -g_z_0_xxyz_yyy[k] * ab_x + g_z_0_xxyz_xyyy[k];

                g_z_0_xxxyz_yyz[k] = -g_z_0_xxyz_yyz[k] * ab_x + g_z_0_xxyz_xyyz[k];

                g_z_0_xxxyz_yzz[k] = -g_z_0_xxyz_yzz[k] * ab_x + g_z_0_xxyz_xyzz[k];

                g_z_0_xxxyz_zzz[k] = -g_z_0_xxyz_zzz[k] * ab_x + g_z_0_xxyz_xzzz[k];
            }

            /// Set up 470-480 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 470 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 471 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 472 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 473 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 474 * ccomps * dcomps);

            auto g_z_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 475 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 476 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 477 * ccomps * dcomps);

            auto g_z_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 478 * ccomps * dcomps);

            auto g_z_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxzz_xxx, g_z_0_xxxzz_xxy, g_z_0_xxxzz_xxz, g_z_0_xxxzz_xyy, g_z_0_xxxzz_xyz, g_z_0_xxxzz_xzz, g_z_0_xxxzz_yyy, g_z_0_xxxzz_yyz, g_z_0_xxxzz_yzz, g_z_0_xxxzz_zzz, g_z_0_xxzz_xxx, g_z_0_xxzz_xxxx, g_z_0_xxzz_xxxy, g_z_0_xxzz_xxxz, g_z_0_xxzz_xxy, g_z_0_xxzz_xxyy, g_z_0_xxzz_xxyz, g_z_0_xxzz_xxz, g_z_0_xxzz_xxzz, g_z_0_xxzz_xyy, g_z_0_xxzz_xyyy, g_z_0_xxzz_xyyz, g_z_0_xxzz_xyz, g_z_0_xxzz_xyzz, g_z_0_xxzz_xzz, g_z_0_xxzz_xzzz, g_z_0_xxzz_yyy, g_z_0_xxzz_yyz, g_z_0_xxzz_yzz, g_z_0_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_xxx[k] = -g_z_0_xxzz_xxx[k] * ab_x + g_z_0_xxzz_xxxx[k];

                g_z_0_xxxzz_xxy[k] = -g_z_0_xxzz_xxy[k] * ab_x + g_z_0_xxzz_xxxy[k];

                g_z_0_xxxzz_xxz[k] = -g_z_0_xxzz_xxz[k] * ab_x + g_z_0_xxzz_xxxz[k];

                g_z_0_xxxzz_xyy[k] = -g_z_0_xxzz_xyy[k] * ab_x + g_z_0_xxzz_xxyy[k];

                g_z_0_xxxzz_xyz[k] = -g_z_0_xxzz_xyz[k] * ab_x + g_z_0_xxzz_xxyz[k];

                g_z_0_xxxzz_xzz[k] = -g_z_0_xxzz_xzz[k] * ab_x + g_z_0_xxzz_xxzz[k];

                g_z_0_xxxzz_yyy[k] = -g_z_0_xxzz_yyy[k] * ab_x + g_z_0_xxzz_xyyy[k];

                g_z_0_xxxzz_yyz[k] = -g_z_0_xxzz_yyz[k] * ab_x + g_z_0_xxzz_xyyz[k];

                g_z_0_xxxzz_yzz[k] = -g_z_0_xxzz_yzz[k] * ab_x + g_z_0_xxzz_xyzz[k];

                g_z_0_xxxzz_zzz[k] = -g_z_0_xxzz_zzz[k] * ab_x + g_z_0_xxzz_xzzz[k];
            }

            /// Set up 480-490 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 480 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 481 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 482 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 483 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 484 * ccomps * dcomps);

            auto g_z_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 485 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 486 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 487 * ccomps * dcomps);

            auto g_z_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 488 * ccomps * dcomps);

            auto g_z_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 489 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyy_xxx, g_z_0_xxyyy_xxy, g_z_0_xxyyy_xxz, g_z_0_xxyyy_xyy, g_z_0_xxyyy_xyz, g_z_0_xxyyy_xzz, g_z_0_xxyyy_yyy, g_z_0_xxyyy_yyz, g_z_0_xxyyy_yzz, g_z_0_xxyyy_zzz, g_z_0_xyyy_xxx, g_z_0_xyyy_xxxx, g_z_0_xyyy_xxxy, g_z_0_xyyy_xxxz, g_z_0_xyyy_xxy, g_z_0_xyyy_xxyy, g_z_0_xyyy_xxyz, g_z_0_xyyy_xxz, g_z_0_xyyy_xxzz, g_z_0_xyyy_xyy, g_z_0_xyyy_xyyy, g_z_0_xyyy_xyyz, g_z_0_xyyy_xyz, g_z_0_xyyy_xyzz, g_z_0_xyyy_xzz, g_z_0_xyyy_xzzz, g_z_0_xyyy_yyy, g_z_0_xyyy_yyz, g_z_0_xyyy_yzz, g_z_0_xyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_xxx[k] = -g_z_0_xyyy_xxx[k] * ab_x + g_z_0_xyyy_xxxx[k];

                g_z_0_xxyyy_xxy[k] = -g_z_0_xyyy_xxy[k] * ab_x + g_z_0_xyyy_xxxy[k];

                g_z_0_xxyyy_xxz[k] = -g_z_0_xyyy_xxz[k] * ab_x + g_z_0_xyyy_xxxz[k];

                g_z_0_xxyyy_xyy[k] = -g_z_0_xyyy_xyy[k] * ab_x + g_z_0_xyyy_xxyy[k];

                g_z_0_xxyyy_xyz[k] = -g_z_0_xyyy_xyz[k] * ab_x + g_z_0_xyyy_xxyz[k];

                g_z_0_xxyyy_xzz[k] = -g_z_0_xyyy_xzz[k] * ab_x + g_z_0_xyyy_xxzz[k];

                g_z_0_xxyyy_yyy[k] = -g_z_0_xyyy_yyy[k] * ab_x + g_z_0_xyyy_xyyy[k];

                g_z_0_xxyyy_yyz[k] = -g_z_0_xyyy_yyz[k] * ab_x + g_z_0_xyyy_xyyz[k];

                g_z_0_xxyyy_yzz[k] = -g_z_0_xyyy_yzz[k] * ab_x + g_z_0_xyyy_xyzz[k];

                g_z_0_xxyyy_zzz[k] = -g_z_0_xyyy_zzz[k] * ab_x + g_z_0_xyyy_xzzz[k];
            }

            /// Set up 490-500 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 490 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 491 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 492 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 493 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 494 * ccomps * dcomps);

            auto g_z_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 495 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 496 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 497 * ccomps * dcomps);

            auto g_z_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 498 * ccomps * dcomps);

            auto g_z_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 499 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyz_xxx, g_z_0_xxyyz_xxy, g_z_0_xxyyz_xxz, g_z_0_xxyyz_xyy, g_z_0_xxyyz_xyz, g_z_0_xxyyz_xzz, g_z_0_xxyyz_yyy, g_z_0_xxyyz_yyz, g_z_0_xxyyz_yzz, g_z_0_xxyyz_zzz, g_z_0_xyyz_xxx, g_z_0_xyyz_xxxx, g_z_0_xyyz_xxxy, g_z_0_xyyz_xxxz, g_z_0_xyyz_xxy, g_z_0_xyyz_xxyy, g_z_0_xyyz_xxyz, g_z_0_xyyz_xxz, g_z_0_xyyz_xxzz, g_z_0_xyyz_xyy, g_z_0_xyyz_xyyy, g_z_0_xyyz_xyyz, g_z_0_xyyz_xyz, g_z_0_xyyz_xyzz, g_z_0_xyyz_xzz, g_z_0_xyyz_xzzz, g_z_0_xyyz_yyy, g_z_0_xyyz_yyz, g_z_0_xyyz_yzz, g_z_0_xyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_xxx[k] = -g_z_0_xyyz_xxx[k] * ab_x + g_z_0_xyyz_xxxx[k];

                g_z_0_xxyyz_xxy[k] = -g_z_0_xyyz_xxy[k] * ab_x + g_z_0_xyyz_xxxy[k];

                g_z_0_xxyyz_xxz[k] = -g_z_0_xyyz_xxz[k] * ab_x + g_z_0_xyyz_xxxz[k];

                g_z_0_xxyyz_xyy[k] = -g_z_0_xyyz_xyy[k] * ab_x + g_z_0_xyyz_xxyy[k];

                g_z_0_xxyyz_xyz[k] = -g_z_0_xyyz_xyz[k] * ab_x + g_z_0_xyyz_xxyz[k];

                g_z_0_xxyyz_xzz[k] = -g_z_0_xyyz_xzz[k] * ab_x + g_z_0_xyyz_xxzz[k];

                g_z_0_xxyyz_yyy[k] = -g_z_0_xyyz_yyy[k] * ab_x + g_z_0_xyyz_xyyy[k];

                g_z_0_xxyyz_yyz[k] = -g_z_0_xyyz_yyz[k] * ab_x + g_z_0_xyyz_xyyz[k];

                g_z_0_xxyyz_yzz[k] = -g_z_0_xyyz_yzz[k] * ab_x + g_z_0_xyyz_xyzz[k];

                g_z_0_xxyyz_zzz[k] = -g_z_0_xyyz_zzz[k] * ab_x + g_z_0_xyyz_xzzz[k];
            }

            /// Set up 500-510 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 500 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 501 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 502 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 503 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 504 * ccomps * dcomps);

            auto g_z_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 505 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 506 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 507 * ccomps * dcomps);

            auto g_z_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 508 * ccomps * dcomps);

            auto g_z_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyzz_xxx, g_z_0_xxyzz_xxy, g_z_0_xxyzz_xxz, g_z_0_xxyzz_xyy, g_z_0_xxyzz_xyz, g_z_0_xxyzz_xzz, g_z_0_xxyzz_yyy, g_z_0_xxyzz_yyz, g_z_0_xxyzz_yzz, g_z_0_xxyzz_zzz, g_z_0_xyzz_xxx, g_z_0_xyzz_xxxx, g_z_0_xyzz_xxxy, g_z_0_xyzz_xxxz, g_z_0_xyzz_xxy, g_z_0_xyzz_xxyy, g_z_0_xyzz_xxyz, g_z_0_xyzz_xxz, g_z_0_xyzz_xxzz, g_z_0_xyzz_xyy, g_z_0_xyzz_xyyy, g_z_0_xyzz_xyyz, g_z_0_xyzz_xyz, g_z_0_xyzz_xyzz, g_z_0_xyzz_xzz, g_z_0_xyzz_xzzz, g_z_0_xyzz_yyy, g_z_0_xyzz_yyz, g_z_0_xyzz_yzz, g_z_0_xyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_xxx[k] = -g_z_0_xyzz_xxx[k] * ab_x + g_z_0_xyzz_xxxx[k];

                g_z_0_xxyzz_xxy[k] = -g_z_0_xyzz_xxy[k] * ab_x + g_z_0_xyzz_xxxy[k];

                g_z_0_xxyzz_xxz[k] = -g_z_0_xyzz_xxz[k] * ab_x + g_z_0_xyzz_xxxz[k];

                g_z_0_xxyzz_xyy[k] = -g_z_0_xyzz_xyy[k] * ab_x + g_z_0_xyzz_xxyy[k];

                g_z_0_xxyzz_xyz[k] = -g_z_0_xyzz_xyz[k] * ab_x + g_z_0_xyzz_xxyz[k];

                g_z_0_xxyzz_xzz[k] = -g_z_0_xyzz_xzz[k] * ab_x + g_z_0_xyzz_xxzz[k];

                g_z_0_xxyzz_yyy[k] = -g_z_0_xyzz_yyy[k] * ab_x + g_z_0_xyzz_xyyy[k];

                g_z_0_xxyzz_yyz[k] = -g_z_0_xyzz_yyz[k] * ab_x + g_z_0_xyzz_xyyz[k];

                g_z_0_xxyzz_yzz[k] = -g_z_0_xyzz_yzz[k] * ab_x + g_z_0_xyzz_xyzz[k];

                g_z_0_xxyzz_zzz[k] = -g_z_0_xyzz_zzz[k] * ab_x + g_z_0_xyzz_xzzz[k];
            }

            /// Set up 510-520 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 510 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 511 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 512 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 513 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 514 * ccomps * dcomps);

            auto g_z_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 515 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 516 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 517 * ccomps * dcomps);

            auto g_z_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 518 * ccomps * dcomps);

            auto g_z_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 519 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzzz_xxx, g_z_0_xxzzz_xxy, g_z_0_xxzzz_xxz, g_z_0_xxzzz_xyy, g_z_0_xxzzz_xyz, g_z_0_xxzzz_xzz, g_z_0_xxzzz_yyy, g_z_0_xxzzz_yyz, g_z_0_xxzzz_yzz, g_z_0_xxzzz_zzz, g_z_0_xzzz_xxx, g_z_0_xzzz_xxxx, g_z_0_xzzz_xxxy, g_z_0_xzzz_xxxz, g_z_0_xzzz_xxy, g_z_0_xzzz_xxyy, g_z_0_xzzz_xxyz, g_z_0_xzzz_xxz, g_z_0_xzzz_xxzz, g_z_0_xzzz_xyy, g_z_0_xzzz_xyyy, g_z_0_xzzz_xyyz, g_z_0_xzzz_xyz, g_z_0_xzzz_xyzz, g_z_0_xzzz_xzz, g_z_0_xzzz_xzzz, g_z_0_xzzz_yyy, g_z_0_xzzz_yyz, g_z_0_xzzz_yzz, g_z_0_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_xxx[k] = -g_z_0_xzzz_xxx[k] * ab_x + g_z_0_xzzz_xxxx[k];

                g_z_0_xxzzz_xxy[k] = -g_z_0_xzzz_xxy[k] * ab_x + g_z_0_xzzz_xxxy[k];

                g_z_0_xxzzz_xxz[k] = -g_z_0_xzzz_xxz[k] * ab_x + g_z_0_xzzz_xxxz[k];

                g_z_0_xxzzz_xyy[k] = -g_z_0_xzzz_xyy[k] * ab_x + g_z_0_xzzz_xxyy[k];

                g_z_0_xxzzz_xyz[k] = -g_z_0_xzzz_xyz[k] * ab_x + g_z_0_xzzz_xxyz[k];

                g_z_0_xxzzz_xzz[k] = -g_z_0_xzzz_xzz[k] * ab_x + g_z_0_xzzz_xxzz[k];

                g_z_0_xxzzz_yyy[k] = -g_z_0_xzzz_yyy[k] * ab_x + g_z_0_xzzz_xyyy[k];

                g_z_0_xxzzz_yyz[k] = -g_z_0_xzzz_yyz[k] * ab_x + g_z_0_xzzz_xyyz[k];

                g_z_0_xxzzz_yzz[k] = -g_z_0_xzzz_yzz[k] * ab_x + g_z_0_xzzz_xyzz[k];

                g_z_0_xxzzz_zzz[k] = -g_z_0_xzzz_zzz[k] * ab_x + g_z_0_xzzz_xzzz[k];
            }

            /// Set up 520-530 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 520 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 521 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 522 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 523 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 524 * ccomps * dcomps);

            auto g_z_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 525 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 526 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 527 * ccomps * dcomps);

            auto g_z_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 528 * ccomps * dcomps);

            auto g_z_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 529 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyy_xxx, g_z_0_xyyyy_xxy, g_z_0_xyyyy_xxz, g_z_0_xyyyy_xyy, g_z_0_xyyyy_xyz, g_z_0_xyyyy_xzz, g_z_0_xyyyy_yyy, g_z_0_xyyyy_yyz, g_z_0_xyyyy_yzz, g_z_0_xyyyy_zzz, g_z_0_yyyy_xxx, g_z_0_yyyy_xxxx, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxxz, g_z_0_yyyy_xxy, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxz, g_z_0_yyyy_xxzz, g_z_0_yyyy_xyy, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xzz, g_z_0_yyyy_xzzz, g_z_0_yyyy_yyy, g_z_0_yyyy_yyz, g_z_0_yyyy_yzz, g_z_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_xxx[k] = -g_z_0_yyyy_xxx[k] * ab_x + g_z_0_yyyy_xxxx[k];

                g_z_0_xyyyy_xxy[k] = -g_z_0_yyyy_xxy[k] * ab_x + g_z_0_yyyy_xxxy[k];

                g_z_0_xyyyy_xxz[k] = -g_z_0_yyyy_xxz[k] * ab_x + g_z_0_yyyy_xxxz[k];

                g_z_0_xyyyy_xyy[k] = -g_z_0_yyyy_xyy[k] * ab_x + g_z_0_yyyy_xxyy[k];

                g_z_0_xyyyy_xyz[k] = -g_z_0_yyyy_xyz[k] * ab_x + g_z_0_yyyy_xxyz[k];

                g_z_0_xyyyy_xzz[k] = -g_z_0_yyyy_xzz[k] * ab_x + g_z_0_yyyy_xxzz[k];

                g_z_0_xyyyy_yyy[k] = -g_z_0_yyyy_yyy[k] * ab_x + g_z_0_yyyy_xyyy[k];

                g_z_0_xyyyy_yyz[k] = -g_z_0_yyyy_yyz[k] * ab_x + g_z_0_yyyy_xyyz[k];

                g_z_0_xyyyy_yzz[k] = -g_z_0_yyyy_yzz[k] * ab_x + g_z_0_yyyy_xyzz[k];

                g_z_0_xyyyy_zzz[k] = -g_z_0_yyyy_zzz[k] * ab_x + g_z_0_yyyy_xzzz[k];
            }

            /// Set up 530-540 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 530 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 531 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 532 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 533 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 534 * ccomps * dcomps);

            auto g_z_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 535 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 536 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 537 * ccomps * dcomps);

            auto g_z_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 538 * ccomps * dcomps);

            auto g_z_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyz_xxx, g_z_0_xyyyz_xxy, g_z_0_xyyyz_xxz, g_z_0_xyyyz_xyy, g_z_0_xyyyz_xyz, g_z_0_xyyyz_xzz, g_z_0_xyyyz_yyy, g_z_0_xyyyz_yyz, g_z_0_xyyyz_yzz, g_z_0_xyyyz_zzz, g_z_0_yyyz_xxx, g_z_0_yyyz_xxxx, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxxz, g_z_0_yyyz_xxy, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxz, g_z_0_yyyz_xxzz, g_z_0_yyyz_xyy, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xzz, g_z_0_yyyz_xzzz, g_z_0_yyyz_yyy, g_z_0_yyyz_yyz, g_z_0_yyyz_yzz, g_z_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_xxx[k] = -g_z_0_yyyz_xxx[k] * ab_x + g_z_0_yyyz_xxxx[k];

                g_z_0_xyyyz_xxy[k] = -g_z_0_yyyz_xxy[k] * ab_x + g_z_0_yyyz_xxxy[k];

                g_z_0_xyyyz_xxz[k] = -g_z_0_yyyz_xxz[k] * ab_x + g_z_0_yyyz_xxxz[k];

                g_z_0_xyyyz_xyy[k] = -g_z_0_yyyz_xyy[k] * ab_x + g_z_0_yyyz_xxyy[k];

                g_z_0_xyyyz_xyz[k] = -g_z_0_yyyz_xyz[k] * ab_x + g_z_0_yyyz_xxyz[k];

                g_z_0_xyyyz_xzz[k] = -g_z_0_yyyz_xzz[k] * ab_x + g_z_0_yyyz_xxzz[k];

                g_z_0_xyyyz_yyy[k] = -g_z_0_yyyz_yyy[k] * ab_x + g_z_0_yyyz_xyyy[k];

                g_z_0_xyyyz_yyz[k] = -g_z_0_yyyz_yyz[k] * ab_x + g_z_0_yyyz_xyyz[k];

                g_z_0_xyyyz_yzz[k] = -g_z_0_yyyz_yzz[k] * ab_x + g_z_0_yyyz_xyzz[k];

                g_z_0_xyyyz_zzz[k] = -g_z_0_yyyz_zzz[k] * ab_x + g_z_0_yyyz_xzzz[k];
            }

            /// Set up 540-550 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 540 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 541 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 542 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 543 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 544 * ccomps * dcomps);

            auto g_z_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 545 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 546 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 547 * ccomps * dcomps);

            auto g_z_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 548 * ccomps * dcomps);

            auto g_z_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 549 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyzz_xxx, g_z_0_xyyzz_xxy, g_z_0_xyyzz_xxz, g_z_0_xyyzz_xyy, g_z_0_xyyzz_xyz, g_z_0_xyyzz_xzz, g_z_0_xyyzz_yyy, g_z_0_xyyzz_yyz, g_z_0_xyyzz_yzz, g_z_0_xyyzz_zzz, g_z_0_yyzz_xxx, g_z_0_yyzz_xxxx, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxxz, g_z_0_yyzz_xxy, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxz, g_z_0_yyzz_xxzz, g_z_0_yyzz_xyy, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xzz, g_z_0_yyzz_xzzz, g_z_0_yyzz_yyy, g_z_0_yyzz_yyz, g_z_0_yyzz_yzz, g_z_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_xxx[k] = -g_z_0_yyzz_xxx[k] * ab_x + g_z_0_yyzz_xxxx[k];

                g_z_0_xyyzz_xxy[k] = -g_z_0_yyzz_xxy[k] * ab_x + g_z_0_yyzz_xxxy[k];

                g_z_0_xyyzz_xxz[k] = -g_z_0_yyzz_xxz[k] * ab_x + g_z_0_yyzz_xxxz[k];

                g_z_0_xyyzz_xyy[k] = -g_z_0_yyzz_xyy[k] * ab_x + g_z_0_yyzz_xxyy[k];

                g_z_0_xyyzz_xyz[k] = -g_z_0_yyzz_xyz[k] * ab_x + g_z_0_yyzz_xxyz[k];

                g_z_0_xyyzz_xzz[k] = -g_z_0_yyzz_xzz[k] * ab_x + g_z_0_yyzz_xxzz[k];

                g_z_0_xyyzz_yyy[k] = -g_z_0_yyzz_yyy[k] * ab_x + g_z_0_yyzz_xyyy[k];

                g_z_0_xyyzz_yyz[k] = -g_z_0_yyzz_yyz[k] * ab_x + g_z_0_yyzz_xyyz[k];

                g_z_0_xyyzz_yzz[k] = -g_z_0_yyzz_yzz[k] * ab_x + g_z_0_yyzz_xyzz[k];

                g_z_0_xyyzz_zzz[k] = -g_z_0_yyzz_zzz[k] * ab_x + g_z_0_yyzz_xzzz[k];
            }

            /// Set up 550-560 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 550 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 551 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 552 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 553 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 554 * ccomps * dcomps);

            auto g_z_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 555 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 556 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 557 * ccomps * dcomps);

            auto g_z_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 558 * ccomps * dcomps);

            auto g_z_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzzz_xxx, g_z_0_xyzzz_xxy, g_z_0_xyzzz_xxz, g_z_0_xyzzz_xyy, g_z_0_xyzzz_xyz, g_z_0_xyzzz_xzz, g_z_0_xyzzz_yyy, g_z_0_xyzzz_yyz, g_z_0_xyzzz_yzz, g_z_0_xyzzz_zzz, g_z_0_yzzz_xxx, g_z_0_yzzz_xxxx, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxxz, g_z_0_yzzz_xxy, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxz, g_z_0_yzzz_xxzz, g_z_0_yzzz_xyy, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xzz, g_z_0_yzzz_xzzz, g_z_0_yzzz_yyy, g_z_0_yzzz_yyz, g_z_0_yzzz_yzz, g_z_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_xxx[k] = -g_z_0_yzzz_xxx[k] * ab_x + g_z_0_yzzz_xxxx[k];

                g_z_0_xyzzz_xxy[k] = -g_z_0_yzzz_xxy[k] * ab_x + g_z_0_yzzz_xxxy[k];

                g_z_0_xyzzz_xxz[k] = -g_z_0_yzzz_xxz[k] * ab_x + g_z_0_yzzz_xxxz[k];

                g_z_0_xyzzz_xyy[k] = -g_z_0_yzzz_xyy[k] * ab_x + g_z_0_yzzz_xxyy[k];

                g_z_0_xyzzz_xyz[k] = -g_z_0_yzzz_xyz[k] * ab_x + g_z_0_yzzz_xxyz[k];

                g_z_0_xyzzz_xzz[k] = -g_z_0_yzzz_xzz[k] * ab_x + g_z_0_yzzz_xxzz[k];

                g_z_0_xyzzz_yyy[k] = -g_z_0_yzzz_yyy[k] * ab_x + g_z_0_yzzz_xyyy[k];

                g_z_0_xyzzz_yyz[k] = -g_z_0_yzzz_yyz[k] * ab_x + g_z_0_yzzz_xyyz[k];

                g_z_0_xyzzz_yzz[k] = -g_z_0_yzzz_yzz[k] * ab_x + g_z_0_yzzz_xyzz[k];

                g_z_0_xyzzz_zzz[k] = -g_z_0_yzzz_zzz[k] * ab_x + g_z_0_yzzz_xzzz[k];
            }

            /// Set up 560-570 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 566 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzzz_xxx, g_z_0_xzzzz_xxy, g_z_0_xzzzz_xxz, g_z_0_xzzzz_xyy, g_z_0_xzzzz_xyz, g_z_0_xzzzz_xzz, g_z_0_xzzzz_yyy, g_z_0_xzzzz_yyz, g_z_0_xzzzz_yzz, g_z_0_xzzzz_zzz, g_z_0_zzzz_xxx, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxy, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyy, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyy, g_z_0_zzzz_yyz, g_z_0_zzzz_yzz, g_z_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_xxx[k] = -g_z_0_zzzz_xxx[k] * ab_x + g_z_0_zzzz_xxxx[k];

                g_z_0_xzzzz_xxy[k] = -g_z_0_zzzz_xxy[k] * ab_x + g_z_0_zzzz_xxxy[k];

                g_z_0_xzzzz_xxz[k] = -g_z_0_zzzz_xxz[k] * ab_x + g_z_0_zzzz_xxxz[k];

                g_z_0_xzzzz_xyy[k] = -g_z_0_zzzz_xyy[k] * ab_x + g_z_0_zzzz_xxyy[k];

                g_z_0_xzzzz_xyz[k] = -g_z_0_zzzz_xyz[k] * ab_x + g_z_0_zzzz_xxyz[k];

                g_z_0_xzzzz_xzz[k] = -g_z_0_zzzz_xzz[k] * ab_x + g_z_0_zzzz_xxzz[k];

                g_z_0_xzzzz_yyy[k] = -g_z_0_zzzz_yyy[k] * ab_x + g_z_0_zzzz_xyyy[k];

                g_z_0_xzzzz_yyz[k] = -g_z_0_zzzz_yyz[k] * ab_x + g_z_0_zzzz_xyyz[k];

                g_z_0_xzzzz_yzz[k] = -g_z_0_zzzz_yzz[k] * ab_x + g_z_0_zzzz_xyzz[k];

                g_z_0_xzzzz_zzz[k] = -g_z_0_zzzz_zzz[k] * ab_x + g_z_0_zzzz_xzzz[k];
            }

            /// Set up 570-580 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 575 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 579 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyy_xxx, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxy, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxz, g_z_0_yyyy_xyy, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xzz, g_z_0_yyyy_yyy, g_z_0_yyyy_yyyy, g_z_0_yyyy_yyyz, g_z_0_yyyy_yyz, g_z_0_yyyy_yyzz, g_z_0_yyyy_yzz, g_z_0_yyyy_yzzz, g_z_0_yyyy_zzz, g_z_0_yyyyy_xxx, g_z_0_yyyyy_xxy, g_z_0_yyyyy_xxz, g_z_0_yyyyy_xyy, g_z_0_yyyyy_xyz, g_z_0_yyyyy_xzz, g_z_0_yyyyy_yyy, g_z_0_yyyyy_yyz, g_z_0_yyyyy_yzz, g_z_0_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_xxx[k] = -g_z_0_yyyy_xxx[k] * ab_y + g_z_0_yyyy_xxxy[k];

                g_z_0_yyyyy_xxy[k] = -g_z_0_yyyy_xxy[k] * ab_y + g_z_0_yyyy_xxyy[k];

                g_z_0_yyyyy_xxz[k] = -g_z_0_yyyy_xxz[k] * ab_y + g_z_0_yyyy_xxyz[k];

                g_z_0_yyyyy_xyy[k] = -g_z_0_yyyy_xyy[k] * ab_y + g_z_0_yyyy_xyyy[k];

                g_z_0_yyyyy_xyz[k] = -g_z_0_yyyy_xyz[k] * ab_y + g_z_0_yyyy_xyyz[k];

                g_z_0_yyyyy_xzz[k] = -g_z_0_yyyy_xzz[k] * ab_y + g_z_0_yyyy_xyzz[k];

                g_z_0_yyyyy_yyy[k] = -g_z_0_yyyy_yyy[k] * ab_y + g_z_0_yyyy_yyyy[k];

                g_z_0_yyyyy_yyz[k] = -g_z_0_yyyy_yyz[k] * ab_y + g_z_0_yyyy_yyyz[k];

                g_z_0_yyyyy_yzz[k] = -g_z_0_yyyy_yzz[k] * ab_y + g_z_0_yyyy_yyzz[k];

                g_z_0_yyyyy_zzz[k] = -g_z_0_yyyy_zzz[k] * ab_y + g_z_0_yyyy_yzzz[k];
            }

            /// Set up 580-590 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 584 * ccomps * dcomps);

            auto g_z_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 587 * ccomps * dcomps);

            auto g_z_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 589 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyz_xxx, g_z_0_yyyyz_xxy, g_z_0_yyyyz_xxz, g_z_0_yyyyz_xyy, g_z_0_yyyyz_xyz, g_z_0_yyyyz_xzz, g_z_0_yyyyz_yyy, g_z_0_yyyyz_yyz, g_z_0_yyyyz_yzz, g_z_0_yyyyz_zzz, g_z_0_yyyz_xxx, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxy, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxz, g_z_0_yyyz_xyy, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xzz, g_z_0_yyyz_yyy, g_z_0_yyyz_yyyy, g_z_0_yyyz_yyyz, g_z_0_yyyz_yyz, g_z_0_yyyz_yyzz, g_z_0_yyyz_yzz, g_z_0_yyyz_yzzz, g_z_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_xxx[k] = -g_z_0_yyyz_xxx[k] * ab_y + g_z_0_yyyz_xxxy[k];

                g_z_0_yyyyz_xxy[k] = -g_z_0_yyyz_xxy[k] * ab_y + g_z_0_yyyz_xxyy[k];

                g_z_0_yyyyz_xxz[k] = -g_z_0_yyyz_xxz[k] * ab_y + g_z_0_yyyz_xxyz[k];

                g_z_0_yyyyz_xyy[k] = -g_z_0_yyyz_xyy[k] * ab_y + g_z_0_yyyz_xyyy[k];

                g_z_0_yyyyz_xyz[k] = -g_z_0_yyyz_xyz[k] * ab_y + g_z_0_yyyz_xyyz[k];

                g_z_0_yyyyz_xzz[k] = -g_z_0_yyyz_xzz[k] * ab_y + g_z_0_yyyz_xyzz[k];

                g_z_0_yyyyz_yyy[k] = -g_z_0_yyyz_yyy[k] * ab_y + g_z_0_yyyz_yyyy[k];

                g_z_0_yyyyz_yyz[k] = -g_z_0_yyyz_yyz[k] * ab_y + g_z_0_yyyz_yyyz[k];

                g_z_0_yyyyz_yzz[k] = -g_z_0_yyyz_yzz[k] * ab_y + g_z_0_yyyz_yyzz[k];

                g_z_0_yyyyz_zzz[k] = -g_z_0_yyyz_zzz[k] * ab_y + g_z_0_yyyz_yzzz[k];
            }

            /// Set up 590-600 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyzz_xxx, g_z_0_yyyzz_xxy, g_z_0_yyyzz_xxz, g_z_0_yyyzz_xyy, g_z_0_yyyzz_xyz, g_z_0_yyyzz_xzz, g_z_0_yyyzz_yyy, g_z_0_yyyzz_yyz, g_z_0_yyyzz_yzz, g_z_0_yyyzz_zzz, g_z_0_yyzz_xxx, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxy, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxz, g_z_0_yyzz_xyy, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xzz, g_z_0_yyzz_yyy, g_z_0_yyzz_yyyy, g_z_0_yyzz_yyyz, g_z_0_yyzz_yyz, g_z_0_yyzz_yyzz, g_z_0_yyzz_yzz, g_z_0_yyzz_yzzz, g_z_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_xxx[k] = -g_z_0_yyzz_xxx[k] * ab_y + g_z_0_yyzz_xxxy[k];

                g_z_0_yyyzz_xxy[k] = -g_z_0_yyzz_xxy[k] * ab_y + g_z_0_yyzz_xxyy[k];

                g_z_0_yyyzz_xxz[k] = -g_z_0_yyzz_xxz[k] * ab_y + g_z_0_yyzz_xxyz[k];

                g_z_0_yyyzz_xyy[k] = -g_z_0_yyzz_xyy[k] * ab_y + g_z_0_yyzz_xyyy[k];

                g_z_0_yyyzz_xyz[k] = -g_z_0_yyzz_xyz[k] * ab_y + g_z_0_yyzz_xyyz[k];

                g_z_0_yyyzz_xzz[k] = -g_z_0_yyzz_xzz[k] * ab_y + g_z_0_yyzz_xyzz[k];

                g_z_0_yyyzz_yyy[k] = -g_z_0_yyzz_yyy[k] * ab_y + g_z_0_yyzz_yyyy[k];

                g_z_0_yyyzz_yyz[k] = -g_z_0_yyzz_yyz[k] * ab_y + g_z_0_yyzz_yyyz[k];

                g_z_0_yyyzz_yzz[k] = -g_z_0_yyzz_yzz[k] * ab_y + g_z_0_yyzz_yyzz[k];

                g_z_0_yyyzz_zzz[k] = -g_z_0_yyzz_zzz[k] * ab_y + g_z_0_yyzz_yzzz[k];
            }

            /// Set up 600-610 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 608 * ccomps * dcomps);

            auto g_z_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 609 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzzz_xxx, g_z_0_yyzzz_xxy, g_z_0_yyzzz_xxz, g_z_0_yyzzz_xyy, g_z_0_yyzzz_xyz, g_z_0_yyzzz_xzz, g_z_0_yyzzz_yyy, g_z_0_yyzzz_yyz, g_z_0_yyzzz_yzz, g_z_0_yyzzz_zzz, g_z_0_yzzz_xxx, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxy, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxz, g_z_0_yzzz_xyy, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xzz, g_z_0_yzzz_yyy, g_z_0_yzzz_yyyy, g_z_0_yzzz_yyyz, g_z_0_yzzz_yyz, g_z_0_yzzz_yyzz, g_z_0_yzzz_yzz, g_z_0_yzzz_yzzz, g_z_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_xxx[k] = -g_z_0_yzzz_xxx[k] * ab_y + g_z_0_yzzz_xxxy[k];

                g_z_0_yyzzz_xxy[k] = -g_z_0_yzzz_xxy[k] * ab_y + g_z_0_yzzz_xxyy[k];

                g_z_0_yyzzz_xxz[k] = -g_z_0_yzzz_xxz[k] * ab_y + g_z_0_yzzz_xxyz[k];

                g_z_0_yyzzz_xyy[k] = -g_z_0_yzzz_xyy[k] * ab_y + g_z_0_yzzz_xyyy[k];

                g_z_0_yyzzz_xyz[k] = -g_z_0_yzzz_xyz[k] * ab_y + g_z_0_yzzz_xyyz[k];

                g_z_0_yyzzz_xzz[k] = -g_z_0_yzzz_xzz[k] * ab_y + g_z_0_yzzz_xyzz[k];

                g_z_0_yyzzz_yyy[k] = -g_z_0_yzzz_yyy[k] * ab_y + g_z_0_yzzz_yyyy[k];

                g_z_0_yyzzz_yyz[k] = -g_z_0_yzzz_yyz[k] * ab_y + g_z_0_yzzz_yyyz[k];

                g_z_0_yyzzz_yzz[k] = -g_z_0_yzzz_yzz[k] * ab_y + g_z_0_yzzz_yyzz[k];

                g_z_0_yyzzz_zzz[k] = -g_z_0_yzzz_zzz[k] * ab_y + g_z_0_yzzz_yzzz[k];
            }

            /// Set up 610-620 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 611 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 614 * ccomps * dcomps);

            auto g_z_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 615 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 619 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzzz_xxx, g_z_0_yzzzz_xxy, g_z_0_yzzzz_xxz, g_z_0_yzzzz_xyy, g_z_0_yzzzz_xyz, g_z_0_yzzzz_xzz, g_z_0_yzzzz_yyy, g_z_0_yzzzz_yyz, g_z_0_yzzzz_yzz, g_z_0_yzzzz_zzz, g_z_0_zzzz_xxx, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxy, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxz, g_z_0_zzzz_xyy, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzz, g_z_0_zzzz_yyy, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_xxx[k] = -g_z_0_zzzz_xxx[k] * ab_y + g_z_0_zzzz_xxxy[k];

                g_z_0_yzzzz_xxy[k] = -g_z_0_zzzz_xxy[k] * ab_y + g_z_0_zzzz_xxyy[k];

                g_z_0_yzzzz_xxz[k] = -g_z_0_zzzz_xxz[k] * ab_y + g_z_0_zzzz_xxyz[k];

                g_z_0_yzzzz_xyy[k] = -g_z_0_zzzz_xyy[k] * ab_y + g_z_0_zzzz_xyyy[k];

                g_z_0_yzzzz_xyz[k] = -g_z_0_zzzz_xyz[k] * ab_y + g_z_0_zzzz_xyyz[k];

                g_z_0_yzzzz_xzz[k] = -g_z_0_zzzz_xzz[k] * ab_y + g_z_0_zzzz_xyzz[k];

                g_z_0_yzzzz_yyy[k] = -g_z_0_zzzz_yyy[k] * ab_y + g_z_0_zzzz_yyyy[k];

                g_z_0_yzzzz_yyz[k] = -g_z_0_zzzz_yyz[k] * ab_y + g_z_0_zzzz_yyyz[k];

                g_z_0_yzzzz_yzz[k] = -g_z_0_zzzz_yzz[k] * ab_y + g_z_0_zzzz_yyzz[k];

                g_z_0_yzzzz_zzz[k] = -g_z_0_zzzz_zzz[k] * ab_y + g_z_0_zzzz_yzzz[k];
            }

            /// Set up 620-630 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzz_xxx, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzz, g_z_0_zzzz_zzzz, g_z_0_zzzzz_xxx, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_yyy, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_zzz, g_zzzz_xxx, g_zzzz_xxy, g_zzzz_xxz, g_zzzz_xyy, g_zzzz_xyz, g_zzzz_xzz, g_zzzz_yyy, g_zzzz_yyz, g_zzzz_yzz, g_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_xxx[k] = -g_zzzz_xxx[k] - g_z_0_zzzz_xxx[k] * ab_z + g_z_0_zzzz_xxxz[k];

                g_z_0_zzzzz_xxy[k] = -g_zzzz_xxy[k] - g_z_0_zzzz_xxy[k] * ab_z + g_z_0_zzzz_xxyz[k];

                g_z_0_zzzzz_xxz[k] = -g_zzzz_xxz[k] - g_z_0_zzzz_xxz[k] * ab_z + g_z_0_zzzz_xxzz[k];

                g_z_0_zzzzz_xyy[k] = -g_zzzz_xyy[k] - g_z_0_zzzz_xyy[k] * ab_z + g_z_0_zzzz_xyyz[k];

                g_z_0_zzzzz_xyz[k] = -g_zzzz_xyz[k] - g_z_0_zzzz_xyz[k] * ab_z + g_z_0_zzzz_xyzz[k];

                g_z_0_zzzzz_xzz[k] = -g_zzzz_xzz[k] - g_z_0_zzzz_xzz[k] * ab_z + g_z_0_zzzz_xzzz[k];

                g_z_0_zzzzz_yyy[k] = -g_zzzz_yyy[k] - g_z_0_zzzz_yyy[k] * ab_z + g_z_0_zzzz_yyyz[k];

                g_z_0_zzzzz_yyz[k] = -g_zzzz_yyz[k] - g_z_0_zzzz_yyz[k] * ab_z + g_z_0_zzzz_yyzz[k];

                g_z_0_zzzzz_yzz[k] = -g_zzzz_yzz[k] - g_z_0_zzzz_yzz[k] * ab_z + g_z_0_zzzz_yzzz[k];

                g_z_0_zzzzz_zzz[k] = -g_zzzz_zzz[k] - g_z_0_zzzz_zzz[k] * ab_z + g_z_0_zzzz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

