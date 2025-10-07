#include "ElectronRepulsionGeom2000ContrRecDFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_dfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_dfxx,
                                            const size_t idx_geom_10_pfxx,
                                            const size_t idx_geom_20_pfxx,
                                            const size_t idx_geom_20_pgxx,
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
            /// Set up components of auxilary buffer : PFSS

            const auto pf_geom_10_off = idx_geom_10_pfxx + i * dcomps + j;

            auto g_x_0_x_xxx = cbuffer.data(pf_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxy = cbuffer.data(pf_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxz = cbuffer.data(pf_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xyy = cbuffer.data(pf_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xyz = cbuffer.data(pf_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xzz = cbuffer.data(pf_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_yyy = cbuffer.data(pf_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_yyz = cbuffer.data(pf_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_yzz = cbuffer.data(pf_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_zzz = cbuffer.data(pf_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_y_xxx = cbuffer.data(pf_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_y_xxy = cbuffer.data(pf_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_y_xxz = cbuffer.data(pf_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_y_xyy = cbuffer.data(pf_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_y_xyz = cbuffer.data(pf_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_y_xzz = cbuffer.data(pf_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_y_yyy = cbuffer.data(pf_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_y_yyz = cbuffer.data(pf_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_y_yzz = cbuffer.data(pf_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_y_zzz = cbuffer.data(pf_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_z_xxx = cbuffer.data(pf_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_z_xxy = cbuffer.data(pf_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_z_xxz = cbuffer.data(pf_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_z_xyy = cbuffer.data(pf_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_z_xyz = cbuffer.data(pf_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_z_xzz = cbuffer.data(pf_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_z_yyy = cbuffer.data(pf_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_z_yyz = cbuffer.data(pf_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_z_yzz = cbuffer.data(pf_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_z_zzz = cbuffer.data(pf_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_x_xxx = cbuffer.data(pf_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_x_xxy = cbuffer.data(pf_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_x_xxz = cbuffer.data(pf_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_x_xyy = cbuffer.data(pf_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_x_xyz = cbuffer.data(pf_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_x_xzz = cbuffer.data(pf_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_x_yyy = cbuffer.data(pf_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_x_yyz = cbuffer.data(pf_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_x_yzz = cbuffer.data(pf_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_x_zzz = cbuffer.data(pf_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_y_xxx = cbuffer.data(pf_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_y_xxy = cbuffer.data(pf_geom_10_off + 41 * ccomps * dcomps);

            auto g_y_0_y_xxz = cbuffer.data(pf_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_y_xyy = cbuffer.data(pf_geom_10_off + 43 * ccomps * dcomps);

            auto g_y_0_y_xyz = cbuffer.data(pf_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_y_xzz = cbuffer.data(pf_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_y_yyy = cbuffer.data(pf_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_y_yyz = cbuffer.data(pf_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_y_yzz = cbuffer.data(pf_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_y_zzz = cbuffer.data(pf_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_z_xxx = cbuffer.data(pf_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_z_xxy = cbuffer.data(pf_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_z_xxz = cbuffer.data(pf_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_z_xyy = cbuffer.data(pf_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_z_xyz = cbuffer.data(pf_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_z_xzz = cbuffer.data(pf_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_z_yyy = cbuffer.data(pf_geom_10_off + 56 * ccomps * dcomps);

            auto g_y_0_z_yyz = cbuffer.data(pf_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_z_yzz = cbuffer.data(pf_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_z_zzz = cbuffer.data(pf_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_x_xxx = cbuffer.data(pf_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_x_xxy = cbuffer.data(pf_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_x_xxz = cbuffer.data(pf_geom_10_off + 62 * ccomps * dcomps);

            auto g_z_0_x_xyy = cbuffer.data(pf_geom_10_off + 63 * ccomps * dcomps);

            auto g_z_0_x_xyz = cbuffer.data(pf_geom_10_off + 64 * ccomps * dcomps);

            auto g_z_0_x_xzz = cbuffer.data(pf_geom_10_off + 65 * ccomps * dcomps);

            auto g_z_0_x_yyy = cbuffer.data(pf_geom_10_off + 66 * ccomps * dcomps);

            auto g_z_0_x_yyz = cbuffer.data(pf_geom_10_off + 67 * ccomps * dcomps);

            auto g_z_0_x_yzz = cbuffer.data(pf_geom_10_off + 68 * ccomps * dcomps);

            auto g_z_0_x_zzz = cbuffer.data(pf_geom_10_off + 69 * ccomps * dcomps);

            auto g_z_0_y_xxx = cbuffer.data(pf_geom_10_off + 70 * ccomps * dcomps);

            auto g_z_0_y_xxy = cbuffer.data(pf_geom_10_off + 71 * ccomps * dcomps);

            auto g_z_0_y_xxz = cbuffer.data(pf_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_y_xyy = cbuffer.data(pf_geom_10_off + 73 * ccomps * dcomps);

            auto g_z_0_y_xyz = cbuffer.data(pf_geom_10_off + 74 * ccomps * dcomps);

            auto g_z_0_y_xzz = cbuffer.data(pf_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_y_yyy = cbuffer.data(pf_geom_10_off + 76 * ccomps * dcomps);

            auto g_z_0_y_yyz = cbuffer.data(pf_geom_10_off + 77 * ccomps * dcomps);

            auto g_z_0_y_yzz = cbuffer.data(pf_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_y_zzz = cbuffer.data(pf_geom_10_off + 79 * ccomps * dcomps);

            auto g_z_0_z_xxx = cbuffer.data(pf_geom_10_off + 80 * ccomps * dcomps);

            auto g_z_0_z_xxy = cbuffer.data(pf_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_z_xxz = cbuffer.data(pf_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_z_xyy = cbuffer.data(pf_geom_10_off + 83 * ccomps * dcomps);

            auto g_z_0_z_xyz = cbuffer.data(pf_geom_10_off + 84 * ccomps * dcomps);

            auto g_z_0_z_xzz = cbuffer.data(pf_geom_10_off + 85 * ccomps * dcomps);

            auto g_z_0_z_yyy = cbuffer.data(pf_geom_10_off + 86 * ccomps * dcomps);

            auto g_z_0_z_yyz = cbuffer.data(pf_geom_10_off + 87 * ccomps * dcomps);

            auto g_z_0_z_yzz = cbuffer.data(pf_geom_10_off + 88 * ccomps * dcomps);

            auto g_z_0_z_zzz = cbuffer.data(pf_geom_10_off + 89 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PFSS

            const auto pf_geom_20_off = idx_geom_20_pfxx + i * dcomps + j;

            auto g_xx_0_x_xxx = cbuffer.data(pf_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_xxy = cbuffer.data(pf_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_xxz = cbuffer.data(pf_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_x_xyy = cbuffer.data(pf_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_x_xyz = cbuffer.data(pf_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_x_xzz = cbuffer.data(pf_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_x_yyy = cbuffer.data(pf_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_x_yyz = cbuffer.data(pf_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_x_yzz = cbuffer.data(pf_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_x_zzz = cbuffer.data(pf_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_y_xxx = cbuffer.data(pf_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_y_xxy = cbuffer.data(pf_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_y_xxz = cbuffer.data(pf_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_y_xyy = cbuffer.data(pf_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_y_xyz = cbuffer.data(pf_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_y_xzz = cbuffer.data(pf_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_y_yyy = cbuffer.data(pf_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_y_yyz = cbuffer.data(pf_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_y_yzz = cbuffer.data(pf_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_y_zzz = cbuffer.data(pf_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_z_xxx = cbuffer.data(pf_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_z_xxy = cbuffer.data(pf_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_z_xxz = cbuffer.data(pf_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_z_xyy = cbuffer.data(pf_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_z_xyz = cbuffer.data(pf_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_z_xzz = cbuffer.data(pf_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_z_yyy = cbuffer.data(pf_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_z_yyz = cbuffer.data(pf_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_z_yzz = cbuffer.data(pf_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_z_zzz = cbuffer.data(pf_geom_20_off + 29 * ccomps * dcomps);

            auto g_xy_0_x_xxx = cbuffer.data(pf_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_x_xxy = cbuffer.data(pf_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_x_xxz = cbuffer.data(pf_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_x_xyy = cbuffer.data(pf_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_x_xyz = cbuffer.data(pf_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_x_xzz = cbuffer.data(pf_geom_20_off + 35 * ccomps * dcomps);

            auto g_xy_0_x_yyy = cbuffer.data(pf_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_x_yyz = cbuffer.data(pf_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_x_yzz = cbuffer.data(pf_geom_20_off + 38 * ccomps * dcomps);

            auto g_xy_0_x_zzz = cbuffer.data(pf_geom_20_off + 39 * ccomps * dcomps);

            auto g_xy_0_y_xxx = cbuffer.data(pf_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_y_xxy = cbuffer.data(pf_geom_20_off + 41 * ccomps * dcomps);

            auto g_xy_0_y_xxz = cbuffer.data(pf_geom_20_off + 42 * ccomps * dcomps);

            auto g_xy_0_y_xyy = cbuffer.data(pf_geom_20_off + 43 * ccomps * dcomps);

            auto g_xy_0_y_xyz = cbuffer.data(pf_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_y_xzz = cbuffer.data(pf_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_y_yyy = cbuffer.data(pf_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_y_yyz = cbuffer.data(pf_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_y_yzz = cbuffer.data(pf_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_y_zzz = cbuffer.data(pf_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_z_xxx = cbuffer.data(pf_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_z_xxy = cbuffer.data(pf_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_z_xxz = cbuffer.data(pf_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_z_xyy = cbuffer.data(pf_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_z_xyz = cbuffer.data(pf_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_z_xzz = cbuffer.data(pf_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_z_yyy = cbuffer.data(pf_geom_20_off + 56 * ccomps * dcomps);

            auto g_xy_0_z_yyz = cbuffer.data(pf_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_z_yzz = cbuffer.data(pf_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_z_zzz = cbuffer.data(pf_geom_20_off + 59 * ccomps * dcomps);

            auto g_xz_0_x_xxx = cbuffer.data(pf_geom_20_off + 60 * ccomps * dcomps);

            auto g_xz_0_x_xxy = cbuffer.data(pf_geom_20_off + 61 * ccomps * dcomps);

            auto g_xz_0_x_xxz = cbuffer.data(pf_geom_20_off + 62 * ccomps * dcomps);

            auto g_xz_0_x_xyy = cbuffer.data(pf_geom_20_off + 63 * ccomps * dcomps);

            auto g_xz_0_x_xyz = cbuffer.data(pf_geom_20_off + 64 * ccomps * dcomps);

            auto g_xz_0_x_xzz = cbuffer.data(pf_geom_20_off + 65 * ccomps * dcomps);

            auto g_xz_0_x_yyy = cbuffer.data(pf_geom_20_off + 66 * ccomps * dcomps);

            auto g_xz_0_x_yyz = cbuffer.data(pf_geom_20_off + 67 * ccomps * dcomps);

            auto g_xz_0_x_yzz = cbuffer.data(pf_geom_20_off + 68 * ccomps * dcomps);

            auto g_xz_0_x_zzz = cbuffer.data(pf_geom_20_off + 69 * ccomps * dcomps);

            auto g_xz_0_y_xxx = cbuffer.data(pf_geom_20_off + 70 * ccomps * dcomps);

            auto g_xz_0_y_xxy = cbuffer.data(pf_geom_20_off + 71 * ccomps * dcomps);

            auto g_xz_0_y_xxz = cbuffer.data(pf_geom_20_off + 72 * ccomps * dcomps);

            auto g_xz_0_y_xyy = cbuffer.data(pf_geom_20_off + 73 * ccomps * dcomps);

            auto g_xz_0_y_xyz = cbuffer.data(pf_geom_20_off + 74 * ccomps * dcomps);

            auto g_xz_0_y_xzz = cbuffer.data(pf_geom_20_off + 75 * ccomps * dcomps);

            auto g_xz_0_y_yyy = cbuffer.data(pf_geom_20_off + 76 * ccomps * dcomps);

            auto g_xz_0_y_yyz = cbuffer.data(pf_geom_20_off + 77 * ccomps * dcomps);

            auto g_xz_0_y_yzz = cbuffer.data(pf_geom_20_off + 78 * ccomps * dcomps);

            auto g_xz_0_y_zzz = cbuffer.data(pf_geom_20_off + 79 * ccomps * dcomps);

            auto g_xz_0_z_xxx = cbuffer.data(pf_geom_20_off + 80 * ccomps * dcomps);

            auto g_xz_0_z_xxy = cbuffer.data(pf_geom_20_off + 81 * ccomps * dcomps);

            auto g_xz_0_z_xxz = cbuffer.data(pf_geom_20_off + 82 * ccomps * dcomps);

            auto g_xz_0_z_xyy = cbuffer.data(pf_geom_20_off + 83 * ccomps * dcomps);

            auto g_xz_0_z_xyz = cbuffer.data(pf_geom_20_off + 84 * ccomps * dcomps);

            auto g_xz_0_z_xzz = cbuffer.data(pf_geom_20_off + 85 * ccomps * dcomps);

            auto g_xz_0_z_yyy = cbuffer.data(pf_geom_20_off + 86 * ccomps * dcomps);

            auto g_xz_0_z_yyz = cbuffer.data(pf_geom_20_off + 87 * ccomps * dcomps);

            auto g_xz_0_z_yzz = cbuffer.data(pf_geom_20_off + 88 * ccomps * dcomps);

            auto g_xz_0_z_zzz = cbuffer.data(pf_geom_20_off + 89 * ccomps * dcomps);

            auto g_yy_0_x_xxx = cbuffer.data(pf_geom_20_off + 90 * ccomps * dcomps);

            auto g_yy_0_x_xxy = cbuffer.data(pf_geom_20_off + 91 * ccomps * dcomps);

            auto g_yy_0_x_xxz = cbuffer.data(pf_geom_20_off + 92 * ccomps * dcomps);

            auto g_yy_0_x_xyy = cbuffer.data(pf_geom_20_off + 93 * ccomps * dcomps);

            auto g_yy_0_x_xyz = cbuffer.data(pf_geom_20_off + 94 * ccomps * dcomps);

            auto g_yy_0_x_xzz = cbuffer.data(pf_geom_20_off + 95 * ccomps * dcomps);

            auto g_yy_0_x_yyy = cbuffer.data(pf_geom_20_off + 96 * ccomps * dcomps);

            auto g_yy_0_x_yyz = cbuffer.data(pf_geom_20_off + 97 * ccomps * dcomps);

            auto g_yy_0_x_yzz = cbuffer.data(pf_geom_20_off + 98 * ccomps * dcomps);

            auto g_yy_0_x_zzz = cbuffer.data(pf_geom_20_off + 99 * ccomps * dcomps);

            auto g_yy_0_y_xxx = cbuffer.data(pf_geom_20_off + 100 * ccomps * dcomps);

            auto g_yy_0_y_xxy = cbuffer.data(pf_geom_20_off + 101 * ccomps * dcomps);

            auto g_yy_0_y_xxz = cbuffer.data(pf_geom_20_off + 102 * ccomps * dcomps);

            auto g_yy_0_y_xyy = cbuffer.data(pf_geom_20_off + 103 * ccomps * dcomps);

            auto g_yy_0_y_xyz = cbuffer.data(pf_geom_20_off + 104 * ccomps * dcomps);

            auto g_yy_0_y_xzz = cbuffer.data(pf_geom_20_off + 105 * ccomps * dcomps);

            auto g_yy_0_y_yyy = cbuffer.data(pf_geom_20_off + 106 * ccomps * dcomps);

            auto g_yy_0_y_yyz = cbuffer.data(pf_geom_20_off + 107 * ccomps * dcomps);

            auto g_yy_0_y_yzz = cbuffer.data(pf_geom_20_off + 108 * ccomps * dcomps);

            auto g_yy_0_y_zzz = cbuffer.data(pf_geom_20_off + 109 * ccomps * dcomps);

            auto g_yy_0_z_xxx = cbuffer.data(pf_geom_20_off + 110 * ccomps * dcomps);

            auto g_yy_0_z_xxy = cbuffer.data(pf_geom_20_off + 111 * ccomps * dcomps);

            auto g_yy_0_z_xxz = cbuffer.data(pf_geom_20_off + 112 * ccomps * dcomps);

            auto g_yy_0_z_xyy = cbuffer.data(pf_geom_20_off + 113 * ccomps * dcomps);

            auto g_yy_0_z_xyz = cbuffer.data(pf_geom_20_off + 114 * ccomps * dcomps);

            auto g_yy_0_z_xzz = cbuffer.data(pf_geom_20_off + 115 * ccomps * dcomps);

            auto g_yy_0_z_yyy = cbuffer.data(pf_geom_20_off + 116 * ccomps * dcomps);

            auto g_yy_0_z_yyz = cbuffer.data(pf_geom_20_off + 117 * ccomps * dcomps);

            auto g_yy_0_z_yzz = cbuffer.data(pf_geom_20_off + 118 * ccomps * dcomps);

            auto g_yy_0_z_zzz = cbuffer.data(pf_geom_20_off + 119 * ccomps * dcomps);

            auto g_yz_0_x_xxx = cbuffer.data(pf_geom_20_off + 120 * ccomps * dcomps);

            auto g_yz_0_x_xxy = cbuffer.data(pf_geom_20_off + 121 * ccomps * dcomps);

            auto g_yz_0_x_xxz = cbuffer.data(pf_geom_20_off + 122 * ccomps * dcomps);

            auto g_yz_0_x_xyy = cbuffer.data(pf_geom_20_off + 123 * ccomps * dcomps);

            auto g_yz_0_x_xyz = cbuffer.data(pf_geom_20_off + 124 * ccomps * dcomps);

            auto g_yz_0_x_xzz = cbuffer.data(pf_geom_20_off + 125 * ccomps * dcomps);

            auto g_yz_0_x_yyy = cbuffer.data(pf_geom_20_off + 126 * ccomps * dcomps);

            auto g_yz_0_x_yyz = cbuffer.data(pf_geom_20_off + 127 * ccomps * dcomps);

            auto g_yz_0_x_yzz = cbuffer.data(pf_geom_20_off + 128 * ccomps * dcomps);

            auto g_yz_0_x_zzz = cbuffer.data(pf_geom_20_off + 129 * ccomps * dcomps);

            auto g_yz_0_y_xxx = cbuffer.data(pf_geom_20_off + 130 * ccomps * dcomps);

            auto g_yz_0_y_xxy = cbuffer.data(pf_geom_20_off + 131 * ccomps * dcomps);

            auto g_yz_0_y_xxz = cbuffer.data(pf_geom_20_off + 132 * ccomps * dcomps);

            auto g_yz_0_y_xyy = cbuffer.data(pf_geom_20_off + 133 * ccomps * dcomps);

            auto g_yz_0_y_xyz = cbuffer.data(pf_geom_20_off + 134 * ccomps * dcomps);

            auto g_yz_0_y_xzz = cbuffer.data(pf_geom_20_off + 135 * ccomps * dcomps);

            auto g_yz_0_y_yyy = cbuffer.data(pf_geom_20_off + 136 * ccomps * dcomps);

            auto g_yz_0_y_yyz = cbuffer.data(pf_geom_20_off + 137 * ccomps * dcomps);

            auto g_yz_0_y_yzz = cbuffer.data(pf_geom_20_off + 138 * ccomps * dcomps);

            auto g_yz_0_y_zzz = cbuffer.data(pf_geom_20_off + 139 * ccomps * dcomps);

            auto g_yz_0_z_xxx = cbuffer.data(pf_geom_20_off + 140 * ccomps * dcomps);

            auto g_yz_0_z_xxy = cbuffer.data(pf_geom_20_off + 141 * ccomps * dcomps);

            auto g_yz_0_z_xxz = cbuffer.data(pf_geom_20_off + 142 * ccomps * dcomps);

            auto g_yz_0_z_xyy = cbuffer.data(pf_geom_20_off + 143 * ccomps * dcomps);

            auto g_yz_0_z_xyz = cbuffer.data(pf_geom_20_off + 144 * ccomps * dcomps);

            auto g_yz_0_z_xzz = cbuffer.data(pf_geom_20_off + 145 * ccomps * dcomps);

            auto g_yz_0_z_yyy = cbuffer.data(pf_geom_20_off + 146 * ccomps * dcomps);

            auto g_yz_0_z_yyz = cbuffer.data(pf_geom_20_off + 147 * ccomps * dcomps);

            auto g_yz_0_z_yzz = cbuffer.data(pf_geom_20_off + 148 * ccomps * dcomps);

            auto g_yz_0_z_zzz = cbuffer.data(pf_geom_20_off + 149 * ccomps * dcomps);

            auto g_zz_0_x_xxx = cbuffer.data(pf_geom_20_off + 150 * ccomps * dcomps);

            auto g_zz_0_x_xxy = cbuffer.data(pf_geom_20_off + 151 * ccomps * dcomps);

            auto g_zz_0_x_xxz = cbuffer.data(pf_geom_20_off + 152 * ccomps * dcomps);

            auto g_zz_0_x_xyy = cbuffer.data(pf_geom_20_off + 153 * ccomps * dcomps);

            auto g_zz_0_x_xyz = cbuffer.data(pf_geom_20_off + 154 * ccomps * dcomps);

            auto g_zz_0_x_xzz = cbuffer.data(pf_geom_20_off + 155 * ccomps * dcomps);

            auto g_zz_0_x_yyy = cbuffer.data(pf_geom_20_off + 156 * ccomps * dcomps);

            auto g_zz_0_x_yyz = cbuffer.data(pf_geom_20_off + 157 * ccomps * dcomps);

            auto g_zz_0_x_yzz = cbuffer.data(pf_geom_20_off + 158 * ccomps * dcomps);

            auto g_zz_0_x_zzz = cbuffer.data(pf_geom_20_off + 159 * ccomps * dcomps);

            auto g_zz_0_y_xxx = cbuffer.data(pf_geom_20_off + 160 * ccomps * dcomps);

            auto g_zz_0_y_xxy = cbuffer.data(pf_geom_20_off + 161 * ccomps * dcomps);

            auto g_zz_0_y_xxz = cbuffer.data(pf_geom_20_off + 162 * ccomps * dcomps);

            auto g_zz_0_y_xyy = cbuffer.data(pf_geom_20_off + 163 * ccomps * dcomps);

            auto g_zz_0_y_xyz = cbuffer.data(pf_geom_20_off + 164 * ccomps * dcomps);

            auto g_zz_0_y_xzz = cbuffer.data(pf_geom_20_off + 165 * ccomps * dcomps);

            auto g_zz_0_y_yyy = cbuffer.data(pf_geom_20_off + 166 * ccomps * dcomps);

            auto g_zz_0_y_yyz = cbuffer.data(pf_geom_20_off + 167 * ccomps * dcomps);

            auto g_zz_0_y_yzz = cbuffer.data(pf_geom_20_off + 168 * ccomps * dcomps);

            auto g_zz_0_y_zzz = cbuffer.data(pf_geom_20_off + 169 * ccomps * dcomps);

            auto g_zz_0_z_xxx = cbuffer.data(pf_geom_20_off + 170 * ccomps * dcomps);

            auto g_zz_0_z_xxy = cbuffer.data(pf_geom_20_off + 171 * ccomps * dcomps);

            auto g_zz_0_z_xxz = cbuffer.data(pf_geom_20_off + 172 * ccomps * dcomps);

            auto g_zz_0_z_xyy = cbuffer.data(pf_geom_20_off + 173 * ccomps * dcomps);

            auto g_zz_0_z_xyz = cbuffer.data(pf_geom_20_off + 174 * ccomps * dcomps);

            auto g_zz_0_z_xzz = cbuffer.data(pf_geom_20_off + 175 * ccomps * dcomps);

            auto g_zz_0_z_yyy = cbuffer.data(pf_geom_20_off + 176 * ccomps * dcomps);

            auto g_zz_0_z_yyz = cbuffer.data(pf_geom_20_off + 177 * ccomps * dcomps);

            auto g_zz_0_z_yzz = cbuffer.data(pf_geom_20_off + 178 * ccomps * dcomps);

            auto g_zz_0_z_zzz = cbuffer.data(pf_geom_20_off + 179 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PGSS

            const auto pg_geom_20_off = idx_geom_20_pgxx + i * dcomps + j;

            auto g_xx_0_x_xxxx = cbuffer.data(pg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_xxxy = cbuffer.data(pg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_xxxz = cbuffer.data(pg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_x_xxyy = cbuffer.data(pg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_x_xxyz = cbuffer.data(pg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_x_xxzz = cbuffer.data(pg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_x_xyyy = cbuffer.data(pg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_x_xyyz = cbuffer.data(pg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_x_xyzz = cbuffer.data(pg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_x_xzzz = cbuffer.data(pg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_x_yyyy = cbuffer.data(pg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_x_yyyz = cbuffer.data(pg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_x_yyzz = cbuffer.data(pg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_x_yzzz = cbuffer.data(pg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_x_zzzz = cbuffer.data(pg_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_y_xxxx = cbuffer.data(pg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_y_xxxy = cbuffer.data(pg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_y_xxxz = cbuffer.data(pg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_y_xxyy = cbuffer.data(pg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_y_xxyz = cbuffer.data(pg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_y_xxzz = cbuffer.data(pg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_y_xyyy = cbuffer.data(pg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_y_xyyz = cbuffer.data(pg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_y_xyzz = cbuffer.data(pg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_y_xzzz = cbuffer.data(pg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_y_yyyy = cbuffer.data(pg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_y_yyyz = cbuffer.data(pg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_y_yyzz = cbuffer.data(pg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_y_yzzz = cbuffer.data(pg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_y_zzzz = cbuffer.data(pg_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_z_xxxx = cbuffer.data(pg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_z_xxxy = cbuffer.data(pg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_z_xxxz = cbuffer.data(pg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_z_xxyy = cbuffer.data(pg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_z_xxyz = cbuffer.data(pg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_z_xxzz = cbuffer.data(pg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_z_xyyy = cbuffer.data(pg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_z_xyyz = cbuffer.data(pg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_z_xyzz = cbuffer.data(pg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_z_xzzz = cbuffer.data(pg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_z_yyyy = cbuffer.data(pg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_z_yyyz = cbuffer.data(pg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_z_yyzz = cbuffer.data(pg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_z_yzzz = cbuffer.data(pg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_z_zzzz = cbuffer.data(pg_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_x_xxxx = cbuffer.data(pg_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_x_xxxy = cbuffer.data(pg_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_x_xxxz = cbuffer.data(pg_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_x_xxyy = cbuffer.data(pg_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_x_xxyz = cbuffer.data(pg_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_x_xxzz = cbuffer.data(pg_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_x_xyyy = cbuffer.data(pg_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_x_xyyz = cbuffer.data(pg_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_x_xyzz = cbuffer.data(pg_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_x_xzzz = cbuffer.data(pg_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_x_yyyy = cbuffer.data(pg_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_x_yyyz = cbuffer.data(pg_geom_20_off + 56 * ccomps * dcomps);

            auto g_xy_0_x_yyzz = cbuffer.data(pg_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_x_yzzz = cbuffer.data(pg_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_x_zzzz = cbuffer.data(pg_geom_20_off + 59 * ccomps * dcomps);

            auto g_xy_0_y_xxxx = cbuffer.data(pg_geom_20_off + 60 * ccomps * dcomps);

            auto g_xy_0_y_xxxy = cbuffer.data(pg_geom_20_off + 61 * ccomps * dcomps);

            auto g_xy_0_y_xxxz = cbuffer.data(pg_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_y_xxyy = cbuffer.data(pg_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_y_xxyz = cbuffer.data(pg_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_y_xxzz = cbuffer.data(pg_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_y_xyyy = cbuffer.data(pg_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_y_xyyz = cbuffer.data(pg_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_y_xyzz = cbuffer.data(pg_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_y_xzzz = cbuffer.data(pg_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_y_yyyy = cbuffer.data(pg_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_y_yyyz = cbuffer.data(pg_geom_20_off + 71 * ccomps * dcomps);

            auto g_xy_0_y_yyzz = cbuffer.data(pg_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_y_yzzz = cbuffer.data(pg_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_y_zzzz = cbuffer.data(pg_geom_20_off + 74 * ccomps * dcomps);

            auto g_xy_0_z_xxxx = cbuffer.data(pg_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_z_xxxy = cbuffer.data(pg_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_z_xxxz = cbuffer.data(pg_geom_20_off + 77 * ccomps * dcomps);

            auto g_xy_0_z_xxyy = cbuffer.data(pg_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_z_xxyz = cbuffer.data(pg_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_z_xxzz = cbuffer.data(pg_geom_20_off + 80 * ccomps * dcomps);

            auto g_xy_0_z_xyyy = cbuffer.data(pg_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_z_xyyz = cbuffer.data(pg_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_z_xyzz = cbuffer.data(pg_geom_20_off + 83 * ccomps * dcomps);

            auto g_xy_0_z_xzzz = cbuffer.data(pg_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_z_yyyy = cbuffer.data(pg_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_z_yyyz = cbuffer.data(pg_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_z_yyzz = cbuffer.data(pg_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_z_yzzz = cbuffer.data(pg_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_z_zzzz = cbuffer.data(pg_geom_20_off + 89 * ccomps * dcomps);

            auto g_xz_0_x_xxxx = cbuffer.data(pg_geom_20_off + 90 * ccomps * dcomps);

            auto g_xz_0_x_xxxy = cbuffer.data(pg_geom_20_off + 91 * ccomps * dcomps);

            auto g_xz_0_x_xxxz = cbuffer.data(pg_geom_20_off + 92 * ccomps * dcomps);

            auto g_xz_0_x_xxyy = cbuffer.data(pg_geom_20_off + 93 * ccomps * dcomps);

            auto g_xz_0_x_xxyz = cbuffer.data(pg_geom_20_off + 94 * ccomps * dcomps);

            auto g_xz_0_x_xxzz = cbuffer.data(pg_geom_20_off + 95 * ccomps * dcomps);

            auto g_xz_0_x_xyyy = cbuffer.data(pg_geom_20_off + 96 * ccomps * dcomps);

            auto g_xz_0_x_xyyz = cbuffer.data(pg_geom_20_off + 97 * ccomps * dcomps);

            auto g_xz_0_x_xyzz = cbuffer.data(pg_geom_20_off + 98 * ccomps * dcomps);

            auto g_xz_0_x_xzzz = cbuffer.data(pg_geom_20_off + 99 * ccomps * dcomps);

            auto g_xz_0_x_yyyy = cbuffer.data(pg_geom_20_off + 100 * ccomps * dcomps);

            auto g_xz_0_x_yyyz = cbuffer.data(pg_geom_20_off + 101 * ccomps * dcomps);

            auto g_xz_0_x_yyzz = cbuffer.data(pg_geom_20_off + 102 * ccomps * dcomps);

            auto g_xz_0_x_yzzz = cbuffer.data(pg_geom_20_off + 103 * ccomps * dcomps);

            auto g_xz_0_x_zzzz = cbuffer.data(pg_geom_20_off + 104 * ccomps * dcomps);

            auto g_xz_0_y_xxxx = cbuffer.data(pg_geom_20_off + 105 * ccomps * dcomps);

            auto g_xz_0_y_xxxy = cbuffer.data(pg_geom_20_off + 106 * ccomps * dcomps);

            auto g_xz_0_y_xxxz = cbuffer.data(pg_geom_20_off + 107 * ccomps * dcomps);

            auto g_xz_0_y_xxyy = cbuffer.data(pg_geom_20_off + 108 * ccomps * dcomps);

            auto g_xz_0_y_xxyz = cbuffer.data(pg_geom_20_off + 109 * ccomps * dcomps);

            auto g_xz_0_y_xxzz = cbuffer.data(pg_geom_20_off + 110 * ccomps * dcomps);

            auto g_xz_0_y_xyyy = cbuffer.data(pg_geom_20_off + 111 * ccomps * dcomps);

            auto g_xz_0_y_xyyz = cbuffer.data(pg_geom_20_off + 112 * ccomps * dcomps);

            auto g_xz_0_y_xyzz = cbuffer.data(pg_geom_20_off + 113 * ccomps * dcomps);

            auto g_xz_0_y_xzzz = cbuffer.data(pg_geom_20_off + 114 * ccomps * dcomps);

            auto g_xz_0_y_yyyy = cbuffer.data(pg_geom_20_off + 115 * ccomps * dcomps);

            auto g_xz_0_y_yyyz = cbuffer.data(pg_geom_20_off + 116 * ccomps * dcomps);

            auto g_xz_0_y_yyzz = cbuffer.data(pg_geom_20_off + 117 * ccomps * dcomps);

            auto g_xz_0_y_yzzz = cbuffer.data(pg_geom_20_off + 118 * ccomps * dcomps);

            auto g_xz_0_y_zzzz = cbuffer.data(pg_geom_20_off + 119 * ccomps * dcomps);

            auto g_xz_0_z_xxxx = cbuffer.data(pg_geom_20_off + 120 * ccomps * dcomps);

            auto g_xz_0_z_xxxy = cbuffer.data(pg_geom_20_off + 121 * ccomps * dcomps);

            auto g_xz_0_z_xxxz = cbuffer.data(pg_geom_20_off + 122 * ccomps * dcomps);

            auto g_xz_0_z_xxyy = cbuffer.data(pg_geom_20_off + 123 * ccomps * dcomps);

            auto g_xz_0_z_xxyz = cbuffer.data(pg_geom_20_off + 124 * ccomps * dcomps);

            auto g_xz_0_z_xxzz = cbuffer.data(pg_geom_20_off + 125 * ccomps * dcomps);

            auto g_xz_0_z_xyyy = cbuffer.data(pg_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_z_xyyz = cbuffer.data(pg_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_z_xyzz = cbuffer.data(pg_geom_20_off + 128 * ccomps * dcomps);

            auto g_xz_0_z_xzzz = cbuffer.data(pg_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_z_yyyy = cbuffer.data(pg_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_z_yyyz = cbuffer.data(pg_geom_20_off + 131 * ccomps * dcomps);

            auto g_xz_0_z_yyzz = cbuffer.data(pg_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_z_yzzz = cbuffer.data(pg_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_z_zzzz = cbuffer.data(pg_geom_20_off + 134 * ccomps * dcomps);

            auto g_yy_0_x_xxxx = cbuffer.data(pg_geom_20_off + 135 * ccomps * dcomps);

            auto g_yy_0_x_xxxy = cbuffer.data(pg_geom_20_off + 136 * ccomps * dcomps);

            auto g_yy_0_x_xxxz = cbuffer.data(pg_geom_20_off + 137 * ccomps * dcomps);

            auto g_yy_0_x_xxyy = cbuffer.data(pg_geom_20_off + 138 * ccomps * dcomps);

            auto g_yy_0_x_xxyz = cbuffer.data(pg_geom_20_off + 139 * ccomps * dcomps);

            auto g_yy_0_x_xxzz = cbuffer.data(pg_geom_20_off + 140 * ccomps * dcomps);

            auto g_yy_0_x_xyyy = cbuffer.data(pg_geom_20_off + 141 * ccomps * dcomps);

            auto g_yy_0_x_xyyz = cbuffer.data(pg_geom_20_off + 142 * ccomps * dcomps);

            auto g_yy_0_x_xyzz = cbuffer.data(pg_geom_20_off + 143 * ccomps * dcomps);

            auto g_yy_0_x_xzzz = cbuffer.data(pg_geom_20_off + 144 * ccomps * dcomps);

            auto g_yy_0_x_yyyy = cbuffer.data(pg_geom_20_off + 145 * ccomps * dcomps);

            auto g_yy_0_x_yyyz = cbuffer.data(pg_geom_20_off + 146 * ccomps * dcomps);

            auto g_yy_0_x_yyzz = cbuffer.data(pg_geom_20_off + 147 * ccomps * dcomps);

            auto g_yy_0_x_yzzz = cbuffer.data(pg_geom_20_off + 148 * ccomps * dcomps);

            auto g_yy_0_x_zzzz = cbuffer.data(pg_geom_20_off + 149 * ccomps * dcomps);

            auto g_yy_0_y_xxxx = cbuffer.data(pg_geom_20_off + 150 * ccomps * dcomps);

            auto g_yy_0_y_xxxy = cbuffer.data(pg_geom_20_off + 151 * ccomps * dcomps);

            auto g_yy_0_y_xxxz = cbuffer.data(pg_geom_20_off + 152 * ccomps * dcomps);

            auto g_yy_0_y_xxyy = cbuffer.data(pg_geom_20_off + 153 * ccomps * dcomps);

            auto g_yy_0_y_xxyz = cbuffer.data(pg_geom_20_off + 154 * ccomps * dcomps);

            auto g_yy_0_y_xxzz = cbuffer.data(pg_geom_20_off + 155 * ccomps * dcomps);

            auto g_yy_0_y_xyyy = cbuffer.data(pg_geom_20_off + 156 * ccomps * dcomps);

            auto g_yy_0_y_xyyz = cbuffer.data(pg_geom_20_off + 157 * ccomps * dcomps);

            auto g_yy_0_y_xyzz = cbuffer.data(pg_geom_20_off + 158 * ccomps * dcomps);

            auto g_yy_0_y_xzzz = cbuffer.data(pg_geom_20_off + 159 * ccomps * dcomps);

            auto g_yy_0_y_yyyy = cbuffer.data(pg_geom_20_off + 160 * ccomps * dcomps);

            auto g_yy_0_y_yyyz = cbuffer.data(pg_geom_20_off + 161 * ccomps * dcomps);

            auto g_yy_0_y_yyzz = cbuffer.data(pg_geom_20_off + 162 * ccomps * dcomps);

            auto g_yy_0_y_yzzz = cbuffer.data(pg_geom_20_off + 163 * ccomps * dcomps);

            auto g_yy_0_y_zzzz = cbuffer.data(pg_geom_20_off + 164 * ccomps * dcomps);

            auto g_yy_0_z_xxxx = cbuffer.data(pg_geom_20_off + 165 * ccomps * dcomps);

            auto g_yy_0_z_xxxy = cbuffer.data(pg_geom_20_off + 166 * ccomps * dcomps);

            auto g_yy_0_z_xxxz = cbuffer.data(pg_geom_20_off + 167 * ccomps * dcomps);

            auto g_yy_0_z_xxyy = cbuffer.data(pg_geom_20_off + 168 * ccomps * dcomps);

            auto g_yy_0_z_xxyz = cbuffer.data(pg_geom_20_off + 169 * ccomps * dcomps);

            auto g_yy_0_z_xxzz = cbuffer.data(pg_geom_20_off + 170 * ccomps * dcomps);

            auto g_yy_0_z_xyyy = cbuffer.data(pg_geom_20_off + 171 * ccomps * dcomps);

            auto g_yy_0_z_xyyz = cbuffer.data(pg_geom_20_off + 172 * ccomps * dcomps);

            auto g_yy_0_z_xyzz = cbuffer.data(pg_geom_20_off + 173 * ccomps * dcomps);

            auto g_yy_0_z_xzzz = cbuffer.data(pg_geom_20_off + 174 * ccomps * dcomps);

            auto g_yy_0_z_yyyy = cbuffer.data(pg_geom_20_off + 175 * ccomps * dcomps);

            auto g_yy_0_z_yyyz = cbuffer.data(pg_geom_20_off + 176 * ccomps * dcomps);

            auto g_yy_0_z_yyzz = cbuffer.data(pg_geom_20_off + 177 * ccomps * dcomps);

            auto g_yy_0_z_yzzz = cbuffer.data(pg_geom_20_off + 178 * ccomps * dcomps);

            auto g_yy_0_z_zzzz = cbuffer.data(pg_geom_20_off + 179 * ccomps * dcomps);

            auto g_yz_0_x_xxxx = cbuffer.data(pg_geom_20_off + 180 * ccomps * dcomps);

            auto g_yz_0_x_xxxy = cbuffer.data(pg_geom_20_off + 181 * ccomps * dcomps);

            auto g_yz_0_x_xxxz = cbuffer.data(pg_geom_20_off + 182 * ccomps * dcomps);

            auto g_yz_0_x_xxyy = cbuffer.data(pg_geom_20_off + 183 * ccomps * dcomps);

            auto g_yz_0_x_xxyz = cbuffer.data(pg_geom_20_off + 184 * ccomps * dcomps);

            auto g_yz_0_x_xxzz = cbuffer.data(pg_geom_20_off + 185 * ccomps * dcomps);

            auto g_yz_0_x_xyyy = cbuffer.data(pg_geom_20_off + 186 * ccomps * dcomps);

            auto g_yz_0_x_xyyz = cbuffer.data(pg_geom_20_off + 187 * ccomps * dcomps);

            auto g_yz_0_x_xyzz = cbuffer.data(pg_geom_20_off + 188 * ccomps * dcomps);

            auto g_yz_0_x_xzzz = cbuffer.data(pg_geom_20_off + 189 * ccomps * dcomps);

            auto g_yz_0_x_yyyy = cbuffer.data(pg_geom_20_off + 190 * ccomps * dcomps);

            auto g_yz_0_x_yyyz = cbuffer.data(pg_geom_20_off + 191 * ccomps * dcomps);

            auto g_yz_0_x_yyzz = cbuffer.data(pg_geom_20_off + 192 * ccomps * dcomps);

            auto g_yz_0_x_yzzz = cbuffer.data(pg_geom_20_off + 193 * ccomps * dcomps);

            auto g_yz_0_x_zzzz = cbuffer.data(pg_geom_20_off + 194 * ccomps * dcomps);

            auto g_yz_0_y_xxxx = cbuffer.data(pg_geom_20_off + 195 * ccomps * dcomps);

            auto g_yz_0_y_xxxy = cbuffer.data(pg_geom_20_off + 196 * ccomps * dcomps);

            auto g_yz_0_y_xxxz = cbuffer.data(pg_geom_20_off + 197 * ccomps * dcomps);

            auto g_yz_0_y_xxyy = cbuffer.data(pg_geom_20_off + 198 * ccomps * dcomps);

            auto g_yz_0_y_xxyz = cbuffer.data(pg_geom_20_off + 199 * ccomps * dcomps);

            auto g_yz_0_y_xxzz = cbuffer.data(pg_geom_20_off + 200 * ccomps * dcomps);

            auto g_yz_0_y_xyyy = cbuffer.data(pg_geom_20_off + 201 * ccomps * dcomps);

            auto g_yz_0_y_xyyz = cbuffer.data(pg_geom_20_off + 202 * ccomps * dcomps);

            auto g_yz_0_y_xyzz = cbuffer.data(pg_geom_20_off + 203 * ccomps * dcomps);

            auto g_yz_0_y_xzzz = cbuffer.data(pg_geom_20_off + 204 * ccomps * dcomps);

            auto g_yz_0_y_yyyy = cbuffer.data(pg_geom_20_off + 205 * ccomps * dcomps);

            auto g_yz_0_y_yyyz = cbuffer.data(pg_geom_20_off + 206 * ccomps * dcomps);

            auto g_yz_0_y_yyzz = cbuffer.data(pg_geom_20_off + 207 * ccomps * dcomps);

            auto g_yz_0_y_yzzz = cbuffer.data(pg_geom_20_off + 208 * ccomps * dcomps);

            auto g_yz_0_y_zzzz = cbuffer.data(pg_geom_20_off + 209 * ccomps * dcomps);

            auto g_yz_0_z_xxxx = cbuffer.data(pg_geom_20_off + 210 * ccomps * dcomps);

            auto g_yz_0_z_xxxy = cbuffer.data(pg_geom_20_off + 211 * ccomps * dcomps);

            auto g_yz_0_z_xxxz = cbuffer.data(pg_geom_20_off + 212 * ccomps * dcomps);

            auto g_yz_0_z_xxyy = cbuffer.data(pg_geom_20_off + 213 * ccomps * dcomps);

            auto g_yz_0_z_xxyz = cbuffer.data(pg_geom_20_off + 214 * ccomps * dcomps);

            auto g_yz_0_z_xxzz = cbuffer.data(pg_geom_20_off + 215 * ccomps * dcomps);

            auto g_yz_0_z_xyyy = cbuffer.data(pg_geom_20_off + 216 * ccomps * dcomps);

            auto g_yz_0_z_xyyz = cbuffer.data(pg_geom_20_off + 217 * ccomps * dcomps);

            auto g_yz_0_z_xyzz = cbuffer.data(pg_geom_20_off + 218 * ccomps * dcomps);

            auto g_yz_0_z_xzzz = cbuffer.data(pg_geom_20_off + 219 * ccomps * dcomps);

            auto g_yz_0_z_yyyy = cbuffer.data(pg_geom_20_off + 220 * ccomps * dcomps);

            auto g_yz_0_z_yyyz = cbuffer.data(pg_geom_20_off + 221 * ccomps * dcomps);

            auto g_yz_0_z_yyzz = cbuffer.data(pg_geom_20_off + 222 * ccomps * dcomps);

            auto g_yz_0_z_yzzz = cbuffer.data(pg_geom_20_off + 223 * ccomps * dcomps);

            auto g_yz_0_z_zzzz = cbuffer.data(pg_geom_20_off + 224 * ccomps * dcomps);

            auto g_zz_0_x_xxxx = cbuffer.data(pg_geom_20_off + 225 * ccomps * dcomps);

            auto g_zz_0_x_xxxy = cbuffer.data(pg_geom_20_off + 226 * ccomps * dcomps);

            auto g_zz_0_x_xxxz = cbuffer.data(pg_geom_20_off + 227 * ccomps * dcomps);

            auto g_zz_0_x_xxyy = cbuffer.data(pg_geom_20_off + 228 * ccomps * dcomps);

            auto g_zz_0_x_xxyz = cbuffer.data(pg_geom_20_off + 229 * ccomps * dcomps);

            auto g_zz_0_x_xxzz = cbuffer.data(pg_geom_20_off + 230 * ccomps * dcomps);

            auto g_zz_0_x_xyyy = cbuffer.data(pg_geom_20_off + 231 * ccomps * dcomps);

            auto g_zz_0_x_xyyz = cbuffer.data(pg_geom_20_off + 232 * ccomps * dcomps);

            auto g_zz_0_x_xyzz = cbuffer.data(pg_geom_20_off + 233 * ccomps * dcomps);

            auto g_zz_0_x_xzzz = cbuffer.data(pg_geom_20_off + 234 * ccomps * dcomps);

            auto g_zz_0_x_yyyy = cbuffer.data(pg_geom_20_off + 235 * ccomps * dcomps);

            auto g_zz_0_x_yyyz = cbuffer.data(pg_geom_20_off + 236 * ccomps * dcomps);

            auto g_zz_0_x_yyzz = cbuffer.data(pg_geom_20_off + 237 * ccomps * dcomps);

            auto g_zz_0_x_yzzz = cbuffer.data(pg_geom_20_off + 238 * ccomps * dcomps);

            auto g_zz_0_x_zzzz = cbuffer.data(pg_geom_20_off + 239 * ccomps * dcomps);

            auto g_zz_0_y_xxxx = cbuffer.data(pg_geom_20_off + 240 * ccomps * dcomps);

            auto g_zz_0_y_xxxy = cbuffer.data(pg_geom_20_off + 241 * ccomps * dcomps);

            auto g_zz_0_y_xxxz = cbuffer.data(pg_geom_20_off + 242 * ccomps * dcomps);

            auto g_zz_0_y_xxyy = cbuffer.data(pg_geom_20_off + 243 * ccomps * dcomps);

            auto g_zz_0_y_xxyz = cbuffer.data(pg_geom_20_off + 244 * ccomps * dcomps);

            auto g_zz_0_y_xxzz = cbuffer.data(pg_geom_20_off + 245 * ccomps * dcomps);

            auto g_zz_0_y_xyyy = cbuffer.data(pg_geom_20_off + 246 * ccomps * dcomps);

            auto g_zz_0_y_xyyz = cbuffer.data(pg_geom_20_off + 247 * ccomps * dcomps);

            auto g_zz_0_y_xyzz = cbuffer.data(pg_geom_20_off + 248 * ccomps * dcomps);

            auto g_zz_0_y_xzzz = cbuffer.data(pg_geom_20_off + 249 * ccomps * dcomps);

            auto g_zz_0_y_yyyy = cbuffer.data(pg_geom_20_off + 250 * ccomps * dcomps);

            auto g_zz_0_y_yyyz = cbuffer.data(pg_geom_20_off + 251 * ccomps * dcomps);

            auto g_zz_0_y_yyzz = cbuffer.data(pg_geom_20_off + 252 * ccomps * dcomps);

            auto g_zz_0_y_yzzz = cbuffer.data(pg_geom_20_off + 253 * ccomps * dcomps);

            auto g_zz_0_y_zzzz = cbuffer.data(pg_geom_20_off + 254 * ccomps * dcomps);

            auto g_zz_0_z_xxxx = cbuffer.data(pg_geom_20_off + 255 * ccomps * dcomps);

            auto g_zz_0_z_xxxy = cbuffer.data(pg_geom_20_off + 256 * ccomps * dcomps);

            auto g_zz_0_z_xxxz = cbuffer.data(pg_geom_20_off + 257 * ccomps * dcomps);

            auto g_zz_0_z_xxyy = cbuffer.data(pg_geom_20_off + 258 * ccomps * dcomps);

            auto g_zz_0_z_xxyz = cbuffer.data(pg_geom_20_off + 259 * ccomps * dcomps);

            auto g_zz_0_z_xxzz = cbuffer.data(pg_geom_20_off + 260 * ccomps * dcomps);

            auto g_zz_0_z_xyyy = cbuffer.data(pg_geom_20_off + 261 * ccomps * dcomps);

            auto g_zz_0_z_xyyz = cbuffer.data(pg_geom_20_off + 262 * ccomps * dcomps);

            auto g_zz_0_z_xyzz = cbuffer.data(pg_geom_20_off + 263 * ccomps * dcomps);

            auto g_zz_0_z_xzzz = cbuffer.data(pg_geom_20_off + 264 * ccomps * dcomps);

            auto g_zz_0_z_yyyy = cbuffer.data(pg_geom_20_off + 265 * ccomps * dcomps);

            auto g_zz_0_z_yyyz = cbuffer.data(pg_geom_20_off + 266 * ccomps * dcomps);

            auto g_zz_0_z_yyzz = cbuffer.data(pg_geom_20_off + 267 * ccomps * dcomps);

            auto g_zz_0_z_yzzz = cbuffer.data(pg_geom_20_off + 268 * ccomps * dcomps);

            auto g_zz_0_z_zzzz = cbuffer.data(pg_geom_20_off + 269 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dfxx

            const auto df_geom_20_off = idx_geom_20_dfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xx_xxx = cbuffer.data(df_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xx_xxy = cbuffer.data(df_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xx_xxz = cbuffer.data(df_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xx_xyy = cbuffer.data(df_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xx_xyz = cbuffer.data(df_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xx_xzz = cbuffer.data(df_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xx_yyy = cbuffer.data(df_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xx_yyz = cbuffer.data(df_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xx_yzz = cbuffer.data(df_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xx_zzz = cbuffer.data(df_geom_20_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxx, g_x_0_x_xxy, g_x_0_x_xxz, g_x_0_x_xyy, g_x_0_x_xyz, g_x_0_x_xzz, g_x_0_x_yyy, g_x_0_x_yyz, g_x_0_x_yzz, g_x_0_x_zzz, g_xx_0_x_xxx, g_xx_0_x_xxxx, g_xx_0_x_xxxy, g_xx_0_x_xxxz, g_xx_0_x_xxy, g_xx_0_x_xxyy, g_xx_0_x_xxyz, g_xx_0_x_xxz, g_xx_0_x_xxzz, g_xx_0_x_xyy, g_xx_0_x_xyyy, g_xx_0_x_xyyz, g_xx_0_x_xyz, g_xx_0_x_xyzz, g_xx_0_x_xzz, g_xx_0_x_xzzz, g_xx_0_x_yyy, g_xx_0_x_yyz, g_xx_0_x_yzz, g_xx_0_x_zzz, g_xx_0_xx_xxx, g_xx_0_xx_xxy, g_xx_0_xx_xxz, g_xx_0_xx_xyy, g_xx_0_xx_xyz, g_xx_0_xx_xzz, g_xx_0_xx_yyy, g_xx_0_xx_yyz, g_xx_0_xx_yzz, g_xx_0_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xx_xxx[k] = -2.0 * g_x_0_x_xxx[k] - g_xx_0_x_xxx[k] * ab_x + g_xx_0_x_xxxx[k];

                g_xx_0_xx_xxy[k] = -2.0 * g_x_0_x_xxy[k] - g_xx_0_x_xxy[k] * ab_x + g_xx_0_x_xxxy[k];

                g_xx_0_xx_xxz[k] = -2.0 * g_x_0_x_xxz[k] - g_xx_0_x_xxz[k] * ab_x + g_xx_0_x_xxxz[k];

                g_xx_0_xx_xyy[k] = -2.0 * g_x_0_x_xyy[k] - g_xx_0_x_xyy[k] * ab_x + g_xx_0_x_xxyy[k];

                g_xx_0_xx_xyz[k] = -2.0 * g_x_0_x_xyz[k] - g_xx_0_x_xyz[k] * ab_x + g_xx_0_x_xxyz[k];

                g_xx_0_xx_xzz[k] = -2.0 * g_x_0_x_xzz[k] - g_xx_0_x_xzz[k] * ab_x + g_xx_0_x_xxzz[k];

                g_xx_0_xx_yyy[k] = -2.0 * g_x_0_x_yyy[k] - g_xx_0_x_yyy[k] * ab_x + g_xx_0_x_xyyy[k];

                g_xx_0_xx_yyz[k] = -2.0 * g_x_0_x_yyz[k] - g_xx_0_x_yyz[k] * ab_x + g_xx_0_x_xyyz[k];

                g_xx_0_xx_yzz[k] = -2.0 * g_x_0_x_yzz[k] - g_xx_0_x_yzz[k] * ab_x + g_xx_0_x_xyzz[k];

                g_xx_0_xx_zzz[k] = -2.0 * g_x_0_x_zzz[k] - g_xx_0_x_zzz[k] * ab_x + g_xx_0_x_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xy_xxx = cbuffer.data(df_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xy_xxy = cbuffer.data(df_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xy_xxz = cbuffer.data(df_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xy_xyy = cbuffer.data(df_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xy_xyz = cbuffer.data(df_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xy_xzz = cbuffer.data(df_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xy_yyy = cbuffer.data(df_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xy_yyz = cbuffer.data(df_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xy_yzz = cbuffer.data(df_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xy_zzz = cbuffer.data(df_geom_20_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_xxx, g_xx_0_x_xxxy, g_xx_0_x_xxy, g_xx_0_x_xxyy, g_xx_0_x_xxyz, g_xx_0_x_xxz, g_xx_0_x_xyy, g_xx_0_x_xyyy, g_xx_0_x_xyyz, g_xx_0_x_xyz, g_xx_0_x_xyzz, g_xx_0_x_xzz, g_xx_0_x_yyy, g_xx_0_x_yyyy, g_xx_0_x_yyyz, g_xx_0_x_yyz, g_xx_0_x_yyzz, g_xx_0_x_yzz, g_xx_0_x_yzzz, g_xx_0_x_zzz, g_xx_0_xy_xxx, g_xx_0_xy_xxy, g_xx_0_xy_xxz, g_xx_0_xy_xyy, g_xx_0_xy_xyz, g_xx_0_xy_xzz, g_xx_0_xy_yyy, g_xx_0_xy_yyz, g_xx_0_xy_yzz, g_xx_0_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xy_xxx[k] = -g_xx_0_x_xxx[k] * ab_y + g_xx_0_x_xxxy[k];

                g_xx_0_xy_xxy[k] = -g_xx_0_x_xxy[k] * ab_y + g_xx_0_x_xxyy[k];

                g_xx_0_xy_xxz[k] = -g_xx_0_x_xxz[k] * ab_y + g_xx_0_x_xxyz[k];

                g_xx_0_xy_xyy[k] = -g_xx_0_x_xyy[k] * ab_y + g_xx_0_x_xyyy[k];

                g_xx_0_xy_xyz[k] = -g_xx_0_x_xyz[k] * ab_y + g_xx_0_x_xyyz[k];

                g_xx_0_xy_xzz[k] = -g_xx_0_x_xzz[k] * ab_y + g_xx_0_x_xyzz[k];

                g_xx_0_xy_yyy[k] = -g_xx_0_x_yyy[k] * ab_y + g_xx_0_x_yyyy[k];

                g_xx_0_xy_yyz[k] = -g_xx_0_x_yyz[k] * ab_y + g_xx_0_x_yyyz[k];

                g_xx_0_xy_yzz[k] = -g_xx_0_x_yzz[k] * ab_y + g_xx_0_x_yyzz[k];

                g_xx_0_xy_zzz[k] = -g_xx_0_x_zzz[k] * ab_y + g_xx_0_x_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xz_xxx = cbuffer.data(df_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xz_xxy = cbuffer.data(df_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xz_xxz = cbuffer.data(df_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xz_xyy = cbuffer.data(df_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xz_xyz = cbuffer.data(df_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xz_xzz = cbuffer.data(df_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xz_yyy = cbuffer.data(df_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xz_yyz = cbuffer.data(df_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xz_yzz = cbuffer.data(df_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xz_zzz = cbuffer.data(df_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_xxx, g_xx_0_x_xxxz, g_xx_0_x_xxy, g_xx_0_x_xxyz, g_xx_0_x_xxz, g_xx_0_x_xxzz, g_xx_0_x_xyy, g_xx_0_x_xyyz, g_xx_0_x_xyz, g_xx_0_x_xyzz, g_xx_0_x_xzz, g_xx_0_x_xzzz, g_xx_0_x_yyy, g_xx_0_x_yyyz, g_xx_0_x_yyz, g_xx_0_x_yyzz, g_xx_0_x_yzz, g_xx_0_x_yzzz, g_xx_0_x_zzz, g_xx_0_x_zzzz, g_xx_0_xz_xxx, g_xx_0_xz_xxy, g_xx_0_xz_xxz, g_xx_0_xz_xyy, g_xx_0_xz_xyz, g_xx_0_xz_xzz, g_xx_0_xz_yyy, g_xx_0_xz_yyz, g_xx_0_xz_yzz, g_xx_0_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xz_xxx[k] = -g_xx_0_x_xxx[k] * ab_z + g_xx_0_x_xxxz[k];

                g_xx_0_xz_xxy[k] = -g_xx_0_x_xxy[k] * ab_z + g_xx_0_x_xxyz[k];

                g_xx_0_xz_xxz[k] = -g_xx_0_x_xxz[k] * ab_z + g_xx_0_x_xxzz[k];

                g_xx_0_xz_xyy[k] = -g_xx_0_x_xyy[k] * ab_z + g_xx_0_x_xyyz[k];

                g_xx_0_xz_xyz[k] = -g_xx_0_x_xyz[k] * ab_z + g_xx_0_x_xyzz[k];

                g_xx_0_xz_xzz[k] = -g_xx_0_x_xzz[k] * ab_z + g_xx_0_x_xzzz[k];

                g_xx_0_xz_yyy[k] = -g_xx_0_x_yyy[k] * ab_z + g_xx_0_x_yyyz[k];

                g_xx_0_xz_yyz[k] = -g_xx_0_x_yyz[k] * ab_z + g_xx_0_x_yyzz[k];

                g_xx_0_xz_yzz[k] = -g_xx_0_x_yzz[k] * ab_z + g_xx_0_x_yzzz[k];

                g_xx_0_xz_zzz[k] = -g_xx_0_x_zzz[k] * ab_z + g_xx_0_x_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yy_xxx = cbuffer.data(df_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_yy_xxy = cbuffer.data(df_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_yy_xxz = cbuffer.data(df_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_yy_xyy = cbuffer.data(df_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_yy_xyz = cbuffer.data(df_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_yy_xzz = cbuffer.data(df_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_yy_yyy = cbuffer.data(df_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_yy_yyz = cbuffer.data(df_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_yy_yzz = cbuffer.data(df_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_yy_zzz = cbuffer.data(df_geom_20_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_y_xxx, g_xx_0_y_xxxy, g_xx_0_y_xxy, g_xx_0_y_xxyy, g_xx_0_y_xxyz, g_xx_0_y_xxz, g_xx_0_y_xyy, g_xx_0_y_xyyy, g_xx_0_y_xyyz, g_xx_0_y_xyz, g_xx_0_y_xyzz, g_xx_0_y_xzz, g_xx_0_y_yyy, g_xx_0_y_yyyy, g_xx_0_y_yyyz, g_xx_0_y_yyz, g_xx_0_y_yyzz, g_xx_0_y_yzz, g_xx_0_y_yzzz, g_xx_0_y_zzz, g_xx_0_yy_xxx, g_xx_0_yy_xxy, g_xx_0_yy_xxz, g_xx_0_yy_xyy, g_xx_0_yy_xyz, g_xx_0_yy_xzz, g_xx_0_yy_yyy, g_xx_0_yy_yyz, g_xx_0_yy_yzz, g_xx_0_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yy_xxx[k] = -g_xx_0_y_xxx[k] * ab_y + g_xx_0_y_xxxy[k];

                g_xx_0_yy_xxy[k] = -g_xx_0_y_xxy[k] * ab_y + g_xx_0_y_xxyy[k];

                g_xx_0_yy_xxz[k] = -g_xx_0_y_xxz[k] * ab_y + g_xx_0_y_xxyz[k];

                g_xx_0_yy_xyy[k] = -g_xx_0_y_xyy[k] * ab_y + g_xx_0_y_xyyy[k];

                g_xx_0_yy_xyz[k] = -g_xx_0_y_xyz[k] * ab_y + g_xx_0_y_xyyz[k];

                g_xx_0_yy_xzz[k] = -g_xx_0_y_xzz[k] * ab_y + g_xx_0_y_xyzz[k];

                g_xx_0_yy_yyy[k] = -g_xx_0_y_yyy[k] * ab_y + g_xx_0_y_yyyy[k];

                g_xx_0_yy_yyz[k] = -g_xx_0_y_yyz[k] * ab_y + g_xx_0_y_yyyz[k];

                g_xx_0_yy_yzz[k] = -g_xx_0_y_yzz[k] * ab_y + g_xx_0_y_yyzz[k];

                g_xx_0_yy_zzz[k] = -g_xx_0_y_zzz[k] * ab_y + g_xx_0_y_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yz_xxx = cbuffer.data(df_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_yz_xxy = cbuffer.data(df_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_yz_xxz = cbuffer.data(df_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_yz_xyy = cbuffer.data(df_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_yz_xyz = cbuffer.data(df_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_yz_xzz = cbuffer.data(df_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_yz_yyy = cbuffer.data(df_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_yz_yyz = cbuffer.data(df_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_yz_yzz = cbuffer.data(df_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_yz_zzz = cbuffer.data(df_geom_20_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yz_xxx, g_xx_0_yz_xxy, g_xx_0_yz_xxz, g_xx_0_yz_xyy, g_xx_0_yz_xyz, g_xx_0_yz_xzz, g_xx_0_yz_yyy, g_xx_0_yz_yyz, g_xx_0_yz_yzz, g_xx_0_yz_zzz, g_xx_0_z_xxx, g_xx_0_z_xxxy, g_xx_0_z_xxy, g_xx_0_z_xxyy, g_xx_0_z_xxyz, g_xx_0_z_xxz, g_xx_0_z_xyy, g_xx_0_z_xyyy, g_xx_0_z_xyyz, g_xx_0_z_xyz, g_xx_0_z_xyzz, g_xx_0_z_xzz, g_xx_0_z_yyy, g_xx_0_z_yyyy, g_xx_0_z_yyyz, g_xx_0_z_yyz, g_xx_0_z_yyzz, g_xx_0_z_yzz, g_xx_0_z_yzzz, g_xx_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yz_xxx[k] = -g_xx_0_z_xxx[k] * ab_y + g_xx_0_z_xxxy[k];

                g_xx_0_yz_xxy[k] = -g_xx_0_z_xxy[k] * ab_y + g_xx_0_z_xxyy[k];

                g_xx_0_yz_xxz[k] = -g_xx_0_z_xxz[k] * ab_y + g_xx_0_z_xxyz[k];

                g_xx_0_yz_xyy[k] = -g_xx_0_z_xyy[k] * ab_y + g_xx_0_z_xyyy[k];

                g_xx_0_yz_xyz[k] = -g_xx_0_z_xyz[k] * ab_y + g_xx_0_z_xyyz[k];

                g_xx_0_yz_xzz[k] = -g_xx_0_z_xzz[k] * ab_y + g_xx_0_z_xyzz[k];

                g_xx_0_yz_yyy[k] = -g_xx_0_z_yyy[k] * ab_y + g_xx_0_z_yyyy[k];

                g_xx_0_yz_yyz[k] = -g_xx_0_z_yyz[k] * ab_y + g_xx_0_z_yyyz[k];

                g_xx_0_yz_yzz[k] = -g_xx_0_z_yzz[k] * ab_y + g_xx_0_z_yyzz[k];

                g_xx_0_yz_zzz[k] = -g_xx_0_z_zzz[k] * ab_y + g_xx_0_z_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zz_xxx = cbuffer.data(df_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_zz_xxy = cbuffer.data(df_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_zz_xxz = cbuffer.data(df_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_zz_xyy = cbuffer.data(df_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_zz_xyz = cbuffer.data(df_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_zz_xzz = cbuffer.data(df_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_zz_yyy = cbuffer.data(df_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_zz_yyz = cbuffer.data(df_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_zz_yzz = cbuffer.data(df_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_zz_zzz = cbuffer.data(df_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_z_xxx, g_xx_0_z_xxxz, g_xx_0_z_xxy, g_xx_0_z_xxyz, g_xx_0_z_xxz, g_xx_0_z_xxzz, g_xx_0_z_xyy, g_xx_0_z_xyyz, g_xx_0_z_xyz, g_xx_0_z_xyzz, g_xx_0_z_xzz, g_xx_0_z_xzzz, g_xx_0_z_yyy, g_xx_0_z_yyyz, g_xx_0_z_yyz, g_xx_0_z_yyzz, g_xx_0_z_yzz, g_xx_0_z_yzzz, g_xx_0_z_zzz, g_xx_0_z_zzzz, g_xx_0_zz_xxx, g_xx_0_zz_xxy, g_xx_0_zz_xxz, g_xx_0_zz_xyy, g_xx_0_zz_xyz, g_xx_0_zz_xzz, g_xx_0_zz_yyy, g_xx_0_zz_yyz, g_xx_0_zz_yzz, g_xx_0_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zz_xxx[k] = -g_xx_0_z_xxx[k] * ab_z + g_xx_0_z_xxxz[k];

                g_xx_0_zz_xxy[k] = -g_xx_0_z_xxy[k] * ab_z + g_xx_0_z_xxyz[k];

                g_xx_0_zz_xxz[k] = -g_xx_0_z_xxz[k] * ab_z + g_xx_0_z_xxzz[k];

                g_xx_0_zz_xyy[k] = -g_xx_0_z_xyy[k] * ab_z + g_xx_0_z_xyyz[k];

                g_xx_0_zz_xyz[k] = -g_xx_0_z_xyz[k] * ab_z + g_xx_0_z_xyzz[k];

                g_xx_0_zz_xzz[k] = -g_xx_0_z_xzz[k] * ab_z + g_xx_0_z_xzzz[k];

                g_xx_0_zz_yyy[k] = -g_xx_0_z_yyy[k] * ab_z + g_xx_0_z_yyyz[k];

                g_xx_0_zz_yyz[k] = -g_xx_0_z_yyz[k] * ab_z + g_xx_0_z_yyzz[k];

                g_xx_0_zz_yzz[k] = -g_xx_0_z_yzz[k] * ab_z + g_xx_0_z_yzzz[k];

                g_xx_0_zz_zzz[k] = -g_xx_0_z_zzz[k] * ab_z + g_xx_0_z_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xx_xxx = cbuffer.data(df_geom_20_off + 60 * ccomps * dcomps);

            auto g_xy_0_xx_xxy = cbuffer.data(df_geom_20_off + 61 * ccomps * dcomps);

            auto g_xy_0_xx_xxz = cbuffer.data(df_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_xx_xyy = cbuffer.data(df_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_xx_xyz = cbuffer.data(df_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_xx_xzz = cbuffer.data(df_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_xx_yyy = cbuffer.data(df_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_xx_yyz = cbuffer.data(df_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_xx_yzz = cbuffer.data(df_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_xx_zzz = cbuffer.data(df_geom_20_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_xxx, g_xy_0_x_xxxx, g_xy_0_x_xxxy, g_xy_0_x_xxxz, g_xy_0_x_xxy, g_xy_0_x_xxyy, g_xy_0_x_xxyz, g_xy_0_x_xxz, g_xy_0_x_xxzz, g_xy_0_x_xyy, g_xy_0_x_xyyy, g_xy_0_x_xyyz, g_xy_0_x_xyz, g_xy_0_x_xyzz, g_xy_0_x_xzz, g_xy_0_x_xzzz, g_xy_0_x_yyy, g_xy_0_x_yyz, g_xy_0_x_yzz, g_xy_0_x_zzz, g_xy_0_xx_xxx, g_xy_0_xx_xxy, g_xy_0_xx_xxz, g_xy_0_xx_xyy, g_xy_0_xx_xyz, g_xy_0_xx_xzz, g_xy_0_xx_yyy, g_xy_0_xx_yyz, g_xy_0_xx_yzz, g_xy_0_xx_zzz, g_y_0_x_xxx, g_y_0_x_xxy, g_y_0_x_xxz, g_y_0_x_xyy, g_y_0_x_xyz, g_y_0_x_xzz, g_y_0_x_yyy, g_y_0_x_yyz, g_y_0_x_yzz, g_y_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xx_xxx[k] = -g_y_0_x_xxx[k] - g_xy_0_x_xxx[k] * ab_x + g_xy_0_x_xxxx[k];

                g_xy_0_xx_xxy[k] = -g_y_0_x_xxy[k] - g_xy_0_x_xxy[k] * ab_x + g_xy_0_x_xxxy[k];

                g_xy_0_xx_xxz[k] = -g_y_0_x_xxz[k] - g_xy_0_x_xxz[k] * ab_x + g_xy_0_x_xxxz[k];

                g_xy_0_xx_xyy[k] = -g_y_0_x_xyy[k] - g_xy_0_x_xyy[k] * ab_x + g_xy_0_x_xxyy[k];

                g_xy_0_xx_xyz[k] = -g_y_0_x_xyz[k] - g_xy_0_x_xyz[k] * ab_x + g_xy_0_x_xxyz[k];

                g_xy_0_xx_xzz[k] = -g_y_0_x_xzz[k] - g_xy_0_x_xzz[k] * ab_x + g_xy_0_x_xxzz[k];

                g_xy_0_xx_yyy[k] = -g_y_0_x_yyy[k] - g_xy_0_x_yyy[k] * ab_x + g_xy_0_x_xyyy[k];

                g_xy_0_xx_yyz[k] = -g_y_0_x_yyz[k] - g_xy_0_x_yyz[k] * ab_x + g_xy_0_x_xyyz[k];

                g_xy_0_xx_yzz[k] = -g_y_0_x_yzz[k] - g_xy_0_x_yzz[k] * ab_x + g_xy_0_x_xyzz[k];

                g_xy_0_xx_zzz[k] = -g_y_0_x_zzz[k] - g_xy_0_x_zzz[k] * ab_x + g_xy_0_x_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xy_xxx = cbuffer.data(df_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_xy_xxy = cbuffer.data(df_geom_20_off + 71 * ccomps * dcomps);

            auto g_xy_0_xy_xxz = cbuffer.data(df_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_xy_xyy = cbuffer.data(df_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_xy_xyz = cbuffer.data(df_geom_20_off + 74 * ccomps * dcomps);

            auto g_xy_0_xy_xzz = cbuffer.data(df_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_xy_yyy = cbuffer.data(df_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_xy_yyz = cbuffer.data(df_geom_20_off + 77 * ccomps * dcomps);

            auto g_xy_0_xy_yzz = cbuffer.data(df_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_xy_zzz = cbuffer.data(df_geom_20_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xy_xxx, g_xy_0_xy_xxy, g_xy_0_xy_xxz, g_xy_0_xy_xyy, g_xy_0_xy_xyz, g_xy_0_xy_xzz, g_xy_0_xy_yyy, g_xy_0_xy_yyz, g_xy_0_xy_yzz, g_xy_0_xy_zzz, g_xy_0_y_xxx, g_xy_0_y_xxxx, g_xy_0_y_xxxy, g_xy_0_y_xxxz, g_xy_0_y_xxy, g_xy_0_y_xxyy, g_xy_0_y_xxyz, g_xy_0_y_xxz, g_xy_0_y_xxzz, g_xy_0_y_xyy, g_xy_0_y_xyyy, g_xy_0_y_xyyz, g_xy_0_y_xyz, g_xy_0_y_xyzz, g_xy_0_y_xzz, g_xy_0_y_xzzz, g_xy_0_y_yyy, g_xy_0_y_yyz, g_xy_0_y_yzz, g_xy_0_y_zzz, g_y_0_y_xxx, g_y_0_y_xxy, g_y_0_y_xxz, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xzz, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yzz, g_y_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xy_xxx[k] = -g_y_0_y_xxx[k] - g_xy_0_y_xxx[k] * ab_x + g_xy_0_y_xxxx[k];

                g_xy_0_xy_xxy[k] = -g_y_0_y_xxy[k] - g_xy_0_y_xxy[k] * ab_x + g_xy_0_y_xxxy[k];

                g_xy_0_xy_xxz[k] = -g_y_0_y_xxz[k] - g_xy_0_y_xxz[k] * ab_x + g_xy_0_y_xxxz[k];

                g_xy_0_xy_xyy[k] = -g_y_0_y_xyy[k] - g_xy_0_y_xyy[k] * ab_x + g_xy_0_y_xxyy[k];

                g_xy_0_xy_xyz[k] = -g_y_0_y_xyz[k] - g_xy_0_y_xyz[k] * ab_x + g_xy_0_y_xxyz[k];

                g_xy_0_xy_xzz[k] = -g_y_0_y_xzz[k] - g_xy_0_y_xzz[k] * ab_x + g_xy_0_y_xxzz[k];

                g_xy_0_xy_yyy[k] = -g_y_0_y_yyy[k] - g_xy_0_y_yyy[k] * ab_x + g_xy_0_y_xyyy[k];

                g_xy_0_xy_yyz[k] = -g_y_0_y_yyz[k] - g_xy_0_y_yyz[k] * ab_x + g_xy_0_y_xyyz[k];

                g_xy_0_xy_yzz[k] = -g_y_0_y_yzz[k] - g_xy_0_y_yzz[k] * ab_x + g_xy_0_y_xyzz[k];

                g_xy_0_xy_zzz[k] = -g_y_0_y_zzz[k] - g_xy_0_y_zzz[k] * ab_x + g_xy_0_y_xzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xz_xxx = cbuffer.data(df_geom_20_off + 80 * ccomps * dcomps);

            auto g_xy_0_xz_xxy = cbuffer.data(df_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_xz_xxz = cbuffer.data(df_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_xz_xyy = cbuffer.data(df_geom_20_off + 83 * ccomps * dcomps);

            auto g_xy_0_xz_xyz = cbuffer.data(df_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_xz_xzz = cbuffer.data(df_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_xz_yyy = cbuffer.data(df_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_xz_yyz = cbuffer.data(df_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_xz_yzz = cbuffer.data(df_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_xz_zzz = cbuffer.data(df_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_xxx, g_xy_0_x_xxxz, g_xy_0_x_xxy, g_xy_0_x_xxyz, g_xy_0_x_xxz, g_xy_0_x_xxzz, g_xy_0_x_xyy, g_xy_0_x_xyyz, g_xy_0_x_xyz, g_xy_0_x_xyzz, g_xy_0_x_xzz, g_xy_0_x_xzzz, g_xy_0_x_yyy, g_xy_0_x_yyyz, g_xy_0_x_yyz, g_xy_0_x_yyzz, g_xy_0_x_yzz, g_xy_0_x_yzzz, g_xy_0_x_zzz, g_xy_0_x_zzzz, g_xy_0_xz_xxx, g_xy_0_xz_xxy, g_xy_0_xz_xxz, g_xy_0_xz_xyy, g_xy_0_xz_xyz, g_xy_0_xz_xzz, g_xy_0_xz_yyy, g_xy_0_xz_yyz, g_xy_0_xz_yzz, g_xy_0_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xz_xxx[k] = -g_xy_0_x_xxx[k] * ab_z + g_xy_0_x_xxxz[k];

                g_xy_0_xz_xxy[k] = -g_xy_0_x_xxy[k] * ab_z + g_xy_0_x_xxyz[k];

                g_xy_0_xz_xxz[k] = -g_xy_0_x_xxz[k] * ab_z + g_xy_0_x_xxzz[k];

                g_xy_0_xz_xyy[k] = -g_xy_0_x_xyy[k] * ab_z + g_xy_0_x_xyyz[k];

                g_xy_0_xz_xyz[k] = -g_xy_0_x_xyz[k] * ab_z + g_xy_0_x_xyzz[k];

                g_xy_0_xz_xzz[k] = -g_xy_0_x_xzz[k] * ab_z + g_xy_0_x_xzzz[k];

                g_xy_0_xz_yyy[k] = -g_xy_0_x_yyy[k] * ab_z + g_xy_0_x_yyyz[k];

                g_xy_0_xz_yyz[k] = -g_xy_0_x_yyz[k] * ab_z + g_xy_0_x_yyzz[k];

                g_xy_0_xz_yzz[k] = -g_xy_0_x_yzz[k] * ab_z + g_xy_0_x_yzzz[k];

                g_xy_0_xz_zzz[k] = -g_xy_0_x_zzz[k] * ab_z + g_xy_0_x_zzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yy_xxx = cbuffer.data(df_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_yy_xxy = cbuffer.data(df_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_yy_xxz = cbuffer.data(df_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_yy_xyy = cbuffer.data(df_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_yy_xyz = cbuffer.data(df_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_yy_xzz = cbuffer.data(df_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_yy_yyy = cbuffer.data(df_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_yy_yyz = cbuffer.data(df_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_yy_yzz = cbuffer.data(df_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_yy_zzz = cbuffer.data(df_geom_20_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxx, g_x_0_y_xxy, g_x_0_y_xxz, g_x_0_y_xyy, g_x_0_y_xyz, g_x_0_y_xzz, g_x_0_y_yyy, g_x_0_y_yyz, g_x_0_y_yzz, g_x_0_y_zzz, g_xy_0_y_xxx, g_xy_0_y_xxxy, g_xy_0_y_xxy, g_xy_0_y_xxyy, g_xy_0_y_xxyz, g_xy_0_y_xxz, g_xy_0_y_xyy, g_xy_0_y_xyyy, g_xy_0_y_xyyz, g_xy_0_y_xyz, g_xy_0_y_xyzz, g_xy_0_y_xzz, g_xy_0_y_yyy, g_xy_0_y_yyyy, g_xy_0_y_yyyz, g_xy_0_y_yyz, g_xy_0_y_yyzz, g_xy_0_y_yzz, g_xy_0_y_yzzz, g_xy_0_y_zzz, g_xy_0_yy_xxx, g_xy_0_yy_xxy, g_xy_0_yy_xxz, g_xy_0_yy_xyy, g_xy_0_yy_xyz, g_xy_0_yy_xzz, g_xy_0_yy_yyy, g_xy_0_yy_yyz, g_xy_0_yy_yzz, g_xy_0_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yy_xxx[k] = -g_x_0_y_xxx[k] - g_xy_0_y_xxx[k] * ab_y + g_xy_0_y_xxxy[k];

                g_xy_0_yy_xxy[k] = -g_x_0_y_xxy[k] - g_xy_0_y_xxy[k] * ab_y + g_xy_0_y_xxyy[k];

                g_xy_0_yy_xxz[k] = -g_x_0_y_xxz[k] - g_xy_0_y_xxz[k] * ab_y + g_xy_0_y_xxyz[k];

                g_xy_0_yy_xyy[k] = -g_x_0_y_xyy[k] - g_xy_0_y_xyy[k] * ab_y + g_xy_0_y_xyyy[k];

                g_xy_0_yy_xyz[k] = -g_x_0_y_xyz[k] - g_xy_0_y_xyz[k] * ab_y + g_xy_0_y_xyyz[k];

                g_xy_0_yy_xzz[k] = -g_x_0_y_xzz[k] - g_xy_0_y_xzz[k] * ab_y + g_xy_0_y_xyzz[k];

                g_xy_0_yy_yyy[k] = -g_x_0_y_yyy[k] - g_xy_0_y_yyy[k] * ab_y + g_xy_0_y_yyyy[k];

                g_xy_0_yy_yyz[k] = -g_x_0_y_yyz[k] - g_xy_0_y_yyz[k] * ab_y + g_xy_0_y_yyyz[k];

                g_xy_0_yy_yzz[k] = -g_x_0_y_yzz[k] - g_xy_0_y_yzz[k] * ab_y + g_xy_0_y_yyzz[k];

                g_xy_0_yy_zzz[k] = -g_x_0_y_zzz[k] - g_xy_0_y_zzz[k] * ab_y + g_xy_0_y_yzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yz_xxx = cbuffer.data(df_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_yz_xxy = cbuffer.data(df_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_yz_xxz = cbuffer.data(df_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_yz_xyy = cbuffer.data(df_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_yz_xyz = cbuffer.data(df_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_yz_xzz = cbuffer.data(df_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_yz_yyy = cbuffer.data(df_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_yz_yyz = cbuffer.data(df_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_yz_yzz = cbuffer.data(df_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_yz_zzz = cbuffer.data(df_geom_20_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_y_xxx, g_xy_0_y_xxxz, g_xy_0_y_xxy, g_xy_0_y_xxyz, g_xy_0_y_xxz, g_xy_0_y_xxzz, g_xy_0_y_xyy, g_xy_0_y_xyyz, g_xy_0_y_xyz, g_xy_0_y_xyzz, g_xy_0_y_xzz, g_xy_0_y_xzzz, g_xy_0_y_yyy, g_xy_0_y_yyyz, g_xy_0_y_yyz, g_xy_0_y_yyzz, g_xy_0_y_yzz, g_xy_0_y_yzzz, g_xy_0_y_zzz, g_xy_0_y_zzzz, g_xy_0_yz_xxx, g_xy_0_yz_xxy, g_xy_0_yz_xxz, g_xy_0_yz_xyy, g_xy_0_yz_xyz, g_xy_0_yz_xzz, g_xy_0_yz_yyy, g_xy_0_yz_yyz, g_xy_0_yz_yzz, g_xy_0_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yz_xxx[k] = -g_xy_0_y_xxx[k] * ab_z + g_xy_0_y_xxxz[k];

                g_xy_0_yz_xxy[k] = -g_xy_0_y_xxy[k] * ab_z + g_xy_0_y_xxyz[k];

                g_xy_0_yz_xxz[k] = -g_xy_0_y_xxz[k] * ab_z + g_xy_0_y_xxzz[k];

                g_xy_0_yz_xyy[k] = -g_xy_0_y_xyy[k] * ab_z + g_xy_0_y_xyyz[k];

                g_xy_0_yz_xyz[k] = -g_xy_0_y_xyz[k] * ab_z + g_xy_0_y_xyzz[k];

                g_xy_0_yz_xzz[k] = -g_xy_0_y_xzz[k] * ab_z + g_xy_0_y_xzzz[k];

                g_xy_0_yz_yyy[k] = -g_xy_0_y_yyy[k] * ab_z + g_xy_0_y_yyyz[k];

                g_xy_0_yz_yyz[k] = -g_xy_0_y_yyz[k] * ab_z + g_xy_0_y_yyzz[k];

                g_xy_0_yz_yzz[k] = -g_xy_0_y_yzz[k] * ab_z + g_xy_0_y_yzzz[k];

                g_xy_0_yz_zzz[k] = -g_xy_0_y_zzz[k] * ab_z + g_xy_0_y_zzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zz_xxx = cbuffer.data(df_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_zz_xxy = cbuffer.data(df_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_zz_xxz = cbuffer.data(df_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_zz_xyy = cbuffer.data(df_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_zz_xyz = cbuffer.data(df_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_zz_xzz = cbuffer.data(df_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_zz_yyy = cbuffer.data(df_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_zz_yyz = cbuffer.data(df_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_zz_yzz = cbuffer.data(df_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_zz_zzz = cbuffer.data(df_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_z_xxx, g_xy_0_z_xxxz, g_xy_0_z_xxy, g_xy_0_z_xxyz, g_xy_0_z_xxz, g_xy_0_z_xxzz, g_xy_0_z_xyy, g_xy_0_z_xyyz, g_xy_0_z_xyz, g_xy_0_z_xyzz, g_xy_0_z_xzz, g_xy_0_z_xzzz, g_xy_0_z_yyy, g_xy_0_z_yyyz, g_xy_0_z_yyz, g_xy_0_z_yyzz, g_xy_0_z_yzz, g_xy_0_z_yzzz, g_xy_0_z_zzz, g_xy_0_z_zzzz, g_xy_0_zz_xxx, g_xy_0_zz_xxy, g_xy_0_zz_xxz, g_xy_0_zz_xyy, g_xy_0_zz_xyz, g_xy_0_zz_xzz, g_xy_0_zz_yyy, g_xy_0_zz_yyz, g_xy_0_zz_yzz, g_xy_0_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zz_xxx[k] = -g_xy_0_z_xxx[k] * ab_z + g_xy_0_z_xxxz[k];

                g_xy_0_zz_xxy[k] = -g_xy_0_z_xxy[k] * ab_z + g_xy_0_z_xxyz[k];

                g_xy_0_zz_xxz[k] = -g_xy_0_z_xxz[k] * ab_z + g_xy_0_z_xxzz[k];

                g_xy_0_zz_xyy[k] = -g_xy_0_z_xyy[k] * ab_z + g_xy_0_z_xyyz[k];

                g_xy_0_zz_xyz[k] = -g_xy_0_z_xyz[k] * ab_z + g_xy_0_z_xyzz[k];

                g_xy_0_zz_xzz[k] = -g_xy_0_z_xzz[k] * ab_z + g_xy_0_z_xzzz[k];

                g_xy_0_zz_yyy[k] = -g_xy_0_z_yyy[k] * ab_z + g_xy_0_z_yyyz[k];

                g_xy_0_zz_yyz[k] = -g_xy_0_z_yyz[k] * ab_z + g_xy_0_z_yyzz[k];

                g_xy_0_zz_yzz[k] = -g_xy_0_z_yzz[k] * ab_z + g_xy_0_z_yzzz[k];

                g_xy_0_zz_zzz[k] = -g_xy_0_z_zzz[k] * ab_z + g_xy_0_z_zzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xx_xxx = cbuffer.data(df_geom_20_off + 120 * ccomps * dcomps);

            auto g_xz_0_xx_xxy = cbuffer.data(df_geom_20_off + 121 * ccomps * dcomps);

            auto g_xz_0_xx_xxz = cbuffer.data(df_geom_20_off + 122 * ccomps * dcomps);

            auto g_xz_0_xx_xyy = cbuffer.data(df_geom_20_off + 123 * ccomps * dcomps);

            auto g_xz_0_xx_xyz = cbuffer.data(df_geom_20_off + 124 * ccomps * dcomps);

            auto g_xz_0_xx_xzz = cbuffer.data(df_geom_20_off + 125 * ccomps * dcomps);

            auto g_xz_0_xx_yyy = cbuffer.data(df_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_xx_yyz = cbuffer.data(df_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_xx_yzz = cbuffer.data(df_geom_20_off + 128 * ccomps * dcomps);

            auto g_xz_0_xx_zzz = cbuffer.data(df_geom_20_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_xxx, g_xz_0_x_xxxx, g_xz_0_x_xxxy, g_xz_0_x_xxxz, g_xz_0_x_xxy, g_xz_0_x_xxyy, g_xz_0_x_xxyz, g_xz_0_x_xxz, g_xz_0_x_xxzz, g_xz_0_x_xyy, g_xz_0_x_xyyy, g_xz_0_x_xyyz, g_xz_0_x_xyz, g_xz_0_x_xyzz, g_xz_0_x_xzz, g_xz_0_x_xzzz, g_xz_0_x_yyy, g_xz_0_x_yyz, g_xz_0_x_yzz, g_xz_0_x_zzz, g_xz_0_xx_xxx, g_xz_0_xx_xxy, g_xz_0_xx_xxz, g_xz_0_xx_xyy, g_xz_0_xx_xyz, g_xz_0_xx_xzz, g_xz_0_xx_yyy, g_xz_0_xx_yyz, g_xz_0_xx_yzz, g_xz_0_xx_zzz, g_z_0_x_xxx, g_z_0_x_xxy, g_z_0_x_xxz, g_z_0_x_xyy, g_z_0_x_xyz, g_z_0_x_xzz, g_z_0_x_yyy, g_z_0_x_yyz, g_z_0_x_yzz, g_z_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xx_xxx[k] = -g_z_0_x_xxx[k] - g_xz_0_x_xxx[k] * ab_x + g_xz_0_x_xxxx[k];

                g_xz_0_xx_xxy[k] = -g_z_0_x_xxy[k] - g_xz_0_x_xxy[k] * ab_x + g_xz_0_x_xxxy[k];

                g_xz_0_xx_xxz[k] = -g_z_0_x_xxz[k] - g_xz_0_x_xxz[k] * ab_x + g_xz_0_x_xxxz[k];

                g_xz_0_xx_xyy[k] = -g_z_0_x_xyy[k] - g_xz_0_x_xyy[k] * ab_x + g_xz_0_x_xxyy[k];

                g_xz_0_xx_xyz[k] = -g_z_0_x_xyz[k] - g_xz_0_x_xyz[k] * ab_x + g_xz_0_x_xxyz[k];

                g_xz_0_xx_xzz[k] = -g_z_0_x_xzz[k] - g_xz_0_x_xzz[k] * ab_x + g_xz_0_x_xxzz[k];

                g_xz_0_xx_yyy[k] = -g_z_0_x_yyy[k] - g_xz_0_x_yyy[k] * ab_x + g_xz_0_x_xyyy[k];

                g_xz_0_xx_yyz[k] = -g_z_0_x_yyz[k] - g_xz_0_x_yyz[k] * ab_x + g_xz_0_x_xyyz[k];

                g_xz_0_xx_yzz[k] = -g_z_0_x_yzz[k] - g_xz_0_x_yzz[k] * ab_x + g_xz_0_x_xyzz[k];

                g_xz_0_xx_zzz[k] = -g_z_0_x_zzz[k] - g_xz_0_x_zzz[k] * ab_x + g_xz_0_x_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xy_xxx = cbuffer.data(df_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_xy_xxy = cbuffer.data(df_geom_20_off + 131 * ccomps * dcomps);

            auto g_xz_0_xy_xxz = cbuffer.data(df_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_xy_xyy = cbuffer.data(df_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_xy_xyz = cbuffer.data(df_geom_20_off + 134 * ccomps * dcomps);

            auto g_xz_0_xy_xzz = cbuffer.data(df_geom_20_off + 135 * ccomps * dcomps);

            auto g_xz_0_xy_yyy = cbuffer.data(df_geom_20_off + 136 * ccomps * dcomps);

            auto g_xz_0_xy_yyz = cbuffer.data(df_geom_20_off + 137 * ccomps * dcomps);

            auto g_xz_0_xy_yzz = cbuffer.data(df_geom_20_off + 138 * ccomps * dcomps);

            auto g_xz_0_xy_zzz = cbuffer.data(df_geom_20_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_xxx, g_xz_0_x_xxxy, g_xz_0_x_xxy, g_xz_0_x_xxyy, g_xz_0_x_xxyz, g_xz_0_x_xxz, g_xz_0_x_xyy, g_xz_0_x_xyyy, g_xz_0_x_xyyz, g_xz_0_x_xyz, g_xz_0_x_xyzz, g_xz_0_x_xzz, g_xz_0_x_yyy, g_xz_0_x_yyyy, g_xz_0_x_yyyz, g_xz_0_x_yyz, g_xz_0_x_yyzz, g_xz_0_x_yzz, g_xz_0_x_yzzz, g_xz_0_x_zzz, g_xz_0_xy_xxx, g_xz_0_xy_xxy, g_xz_0_xy_xxz, g_xz_0_xy_xyy, g_xz_0_xy_xyz, g_xz_0_xy_xzz, g_xz_0_xy_yyy, g_xz_0_xy_yyz, g_xz_0_xy_yzz, g_xz_0_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xy_xxx[k] = -g_xz_0_x_xxx[k] * ab_y + g_xz_0_x_xxxy[k];

                g_xz_0_xy_xxy[k] = -g_xz_0_x_xxy[k] * ab_y + g_xz_0_x_xxyy[k];

                g_xz_0_xy_xxz[k] = -g_xz_0_x_xxz[k] * ab_y + g_xz_0_x_xxyz[k];

                g_xz_0_xy_xyy[k] = -g_xz_0_x_xyy[k] * ab_y + g_xz_0_x_xyyy[k];

                g_xz_0_xy_xyz[k] = -g_xz_0_x_xyz[k] * ab_y + g_xz_0_x_xyyz[k];

                g_xz_0_xy_xzz[k] = -g_xz_0_x_xzz[k] * ab_y + g_xz_0_x_xyzz[k];

                g_xz_0_xy_yyy[k] = -g_xz_0_x_yyy[k] * ab_y + g_xz_0_x_yyyy[k];

                g_xz_0_xy_yyz[k] = -g_xz_0_x_yyz[k] * ab_y + g_xz_0_x_yyyz[k];

                g_xz_0_xy_yzz[k] = -g_xz_0_x_yzz[k] * ab_y + g_xz_0_x_yyzz[k];

                g_xz_0_xy_zzz[k] = -g_xz_0_x_zzz[k] * ab_y + g_xz_0_x_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xz_xxx = cbuffer.data(df_geom_20_off + 140 * ccomps * dcomps);

            auto g_xz_0_xz_xxy = cbuffer.data(df_geom_20_off + 141 * ccomps * dcomps);

            auto g_xz_0_xz_xxz = cbuffer.data(df_geom_20_off + 142 * ccomps * dcomps);

            auto g_xz_0_xz_xyy = cbuffer.data(df_geom_20_off + 143 * ccomps * dcomps);

            auto g_xz_0_xz_xyz = cbuffer.data(df_geom_20_off + 144 * ccomps * dcomps);

            auto g_xz_0_xz_xzz = cbuffer.data(df_geom_20_off + 145 * ccomps * dcomps);

            auto g_xz_0_xz_yyy = cbuffer.data(df_geom_20_off + 146 * ccomps * dcomps);

            auto g_xz_0_xz_yyz = cbuffer.data(df_geom_20_off + 147 * ccomps * dcomps);

            auto g_xz_0_xz_yzz = cbuffer.data(df_geom_20_off + 148 * ccomps * dcomps);

            auto g_xz_0_xz_zzz = cbuffer.data(df_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xz_xxx, g_xz_0_xz_xxy, g_xz_0_xz_xxz, g_xz_0_xz_xyy, g_xz_0_xz_xyz, g_xz_0_xz_xzz, g_xz_0_xz_yyy, g_xz_0_xz_yyz, g_xz_0_xz_yzz, g_xz_0_xz_zzz, g_xz_0_z_xxx, g_xz_0_z_xxxx, g_xz_0_z_xxxy, g_xz_0_z_xxxz, g_xz_0_z_xxy, g_xz_0_z_xxyy, g_xz_0_z_xxyz, g_xz_0_z_xxz, g_xz_0_z_xxzz, g_xz_0_z_xyy, g_xz_0_z_xyyy, g_xz_0_z_xyyz, g_xz_0_z_xyz, g_xz_0_z_xyzz, g_xz_0_z_xzz, g_xz_0_z_xzzz, g_xz_0_z_yyy, g_xz_0_z_yyz, g_xz_0_z_yzz, g_xz_0_z_zzz, g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xz_xxx[k] = -g_z_0_z_xxx[k] - g_xz_0_z_xxx[k] * ab_x + g_xz_0_z_xxxx[k];

                g_xz_0_xz_xxy[k] = -g_z_0_z_xxy[k] - g_xz_0_z_xxy[k] * ab_x + g_xz_0_z_xxxy[k];

                g_xz_0_xz_xxz[k] = -g_z_0_z_xxz[k] - g_xz_0_z_xxz[k] * ab_x + g_xz_0_z_xxxz[k];

                g_xz_0_xz_xyy[k] = -g_z_0_z_xyy[k] - g_xz_0_z_xyy[k] * ab_x + g_xz_0_z_xxyy[k];

                g_xz_0_xz_xyz[k] = -g_z_0_z_xyz[k] - g_xz_0_z_xyz[k] * ab_x + g_xz_0_z_xxyz[k];

                g_xz_0_xz_xzz[k] = -g_z_0_z_xzz[k] - g_xz_0_z_xzz[k] * ab_x + g_xz_0_z_xxzz[k];

                g_xz_0_xz_yyy[k] = -g_z_0_z_yyy[k] - g_xz_0_z_yyy[k] * ab_x + g_xz_0_z_xyyy[k];

                g_xz_0_xz_yyz[k] = -g_z_0_z_yyz[k] - g_xz_0_z_yyz[k] * ab_x + g_xz_0_z_xyyz[k];

                g_xz_0_xz_yzz[k] = -g_z_0_z_yzz[k] - g_xz_0_z_yzz[k] * ab_x + g_xz_0_z_xyzz[k];

                g_xz_0_xz_zzz[k] = -g_z_0_z_zzz[k] - g_xz_0_z_zzz[k] * ab_x + g_xz_0_z_xzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yy_xxx = cbuffer.data(df_geom_20_off + 150 * ccomps * dcomps);

            auto g_xz_0_yy_xxy = cbuffer.data(df_geom_20_off + 151 * ccomps * dcomps);

            auto g_xz_0_yy_xxz = cbuffer.data(df_geom_20_off + 152 * ccomps * dcomps);

            auto g_xz_0_yy_xyy = cbuffer.data(df_geom_20_off + 153 * ccomps * dcomps);

            auto g_xz_0_yy_xyz = cbuffer.data(df_geom_20_off + 154 * ccomps * dcomps);

            auto g_xz_0_yy_xzz = cbuffer.data(df_geom_20_off + 155 * ccomps * dcomps);

            auto g_xz_0_yy_yyy = cbuffer.data(df_geom_20_off + 156 * ccomps * dcomps);

            auto g_xz_0_yy_yyz = cbuffer.data(df_geom_20_off + 157 * ccomps * dcomps);

            auto g_xz_0_yy_yzz = cbuffer.data(df_geom_20_off + 158 * ccomps * dcomps);

            auto g_xz_0_yy_zzz = cbuffer.data(df_geom_20_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_y_xxx, g_xz_0_y_xxxy, g_xz_0_y_xxy, g_xz_0_y_xxyy, g_xz_0_y_xxyz, g_xz_0_y_xxz, g_xz_0_y_xyy, g_xz_0_y_xyyy, g_xz_0_y_xyyz, g_xz_0_y_xyz, g_xz_0_y_xyzz, g_xz_0_y_xzz, g_xz_0_y_yyy, g_xz_0_y_yyyy, g_xz_0_y_yyyz, g_xz_0_y_yyz, g_xz_0_y_yyzz, g_xz_0_y_yzz, g_xz_0_y_yzzz, g_xz_0_y_zzz, g_xz_0_yy_xxx, g_xz_0_yy_xxy, g_xz_0_yy_xxz, g_xz_0_yy_xyy, g_xz_0_yy_xyz, g_xz_0_yy_xzz, g_xz_0_yy_yyy, g_xz_0_yy_yyz, g_xz_0_yy_yzz, g_xz_0_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yy_xxx[k] = -g_xz_0_y_xxx[k] * ab_y + g_xz_0_y_xxxy[k];

                g_xz_0_yy_xxy[k] = -g_xz_0_y_xxy[k] * ab_y + g_xz_0_y_xxyy[k];

                g_xz_0_yy_xxz[k] = -g_xz_0_y_xxz[k] * ab_y + g_xz_0_y_xxyz[k];

                g_xz_0_yy_xyy[k] = -g_xz_0_y_xyy[k] * ab_y + g_xz_0_y_xyyy[k];

                g_xz_0_yy_xyz[k] = -g_xz_0_y_xyz[k] * ab_y + g_xz_0_y_xyyz[k];

                g_xz_0_yy_xzz[k] = -g_xz_0_y_xzz[k] * ab_y + g_xz_0_y_xyzz[k];

                g_xz_0_yy_yyy[k] = -g_xz_0_y_yyy[k] * ab_y + g_xz_0_y_yyyy[k];

                g_xz_0_yy_yyz[k] = -g_xz_0_y_yyz[k] * ab_y + g_xz_0_y_yyyz[k];

                g_xz_0_yy_yzz[k] = -g_xz_0_y_yzz[k] * ab_y + g_xz_0_y_yyzz[k];

                g_xz_0_yy_zzz[k] = -g_xz_0_y_zzz[k] * ab_y + g_xz_0_y_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yz_xxx = cbuffer.data(df_geom_20_off + 160 * ccomps * dcomps);

            auto g_xz_0_yz_xxy = cbuffer.data(df_geom_20_off + 161 * ccomps * dcomps);

            auto g_xz_0_yz_xxz = cbuffer.data(df_geom_20_off + 162 * ccomps * dcomps);

            auto g_xz_0_yz_xyy = cbuffer.data(df_geom_20_off + 163 * ccomps * dcomps);

            auto g_xz_0_yz_xyz = cbuffer.data(df_geom_20_off + 164 * ccomps * dcomps);

            auto g_xz_0_yz_xzz = cbuffer.data(df_geom_20_off + 165 * ccomps * dcomps);

            auto g_xz_0_yz_yyy = cbuffer.data(df_geom_20_off + 166 * ccomps * dcomps);

            auto g_xz_0_yz_yyz = cbuffer.data(df_geom_20_off + 167 * ccomps * dcomps);

            auto g_xz_0_yz_yzz = cbuffer.data(df_geom_20_off + 168 * ccomps * dcomps);

            auto g_xz_0_yz_zzz = cbuffer.data(df_geom_20_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yz_xxx, g_xz_0_yz_xxy, g_xz_0_yz_xxz, g_xz_0_yz_xyy, g_xz_0_yz_xyz, g_xz_0_yz_xzz, g_xz_0_yz_yyy, g_xz_0_yz_yyz, g_xz_0_yz_yzz, g_xz_0_yz_zzz, g_xz_0_z_xxx, g_xz_0_z_xxxy, g_xz_0_z_xxy, g_xz_0_z_xxyy, g_xz_0_z_xxyz, g_xz_0_z_xxz, g_xz_0_z_xyy, g_xz_0_z_xyyy, g_xz_0_z_xyyz, g_xz_0_z_xyz, g_xz_0_z_xyzz, g_xz_0_z_xzz, g_xz_0_z_yyy, g_xz_0_z_yyyy, g_xz_0_z_yyyz, g_xz_0_z_yyz, g_xz_0_z_yyzz, g_xz_0_z_yzz, g_xz_0_z_yzzz, g_xz_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yz_xxx[k] = -g_xz_0_z_xxx[k] * ab_y + g_xz_0_z_xxxy[k];

                g_xz_0_yz_xxy[k] = -g_xz_0_z_xxy[k] * ab_y + g_xz_0_z_xxyy[k];

                g_xz_0_yz_xxz[k] = -g_xz_0_z_xxz[k] * ab_y + g_xz_0_z_xxyz[k];

                g_xz_0_yz_xyy[k] = -g_xz_0_z_xyy[k] * ab_y + g_xz_0_z_xyyy[k];

                g_xz_0_yz_xyz[k] = -g_xz_0_z_xyz[k] * ab_y + g_xz_0_z_xyyz[k];

                g_xz_0_yz_xzz[k] = -g_xz_0_z_xzz[k] * ab_y + g_xz_0_z_xyzz[k];

                g_xz_0_yz_yyy[k] = -g_xz_0_z_yyy[k] * ab_y + g_xz_0_z_yyyy[k];

                g_xz_0_yz_yyz[k] = -g_xz_0_z_yyz[k] * ab_y + g_xz_0_z_yyyz[k];

                g_xz_0_yz_yzz[k] = -g_xz_0_z_yzz[k] * ab_y + g_xz_0_z_yyzz[k];

                g_xz_0_yz_zzz[k] = -g_xz_0_z_zzz[k] * ab_y + g_xz_0_z_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zz_xxx = cbuffer.data(df_geom_20_off + 170 * ccomps * dcomps);

            auto g_xz_0_zz_xxy = cbuffer.data(df_geom_20_off + 171 * ccomps * dcomps);

            auto g_xz_0_zz_xxz = cbuffer.data(df_geom_20_off + 172 * ccomps * dcomps);

            auto g_xz_0_zz_xyy = cbuffer.data(df_geom_20_off + 173 * ccomps * dcomps);

            auto g_xz_0_zz_xyz = cbuffer.data(df_geom_20_off + 174 * ccomps * dcomps);

            auto g_xz_0_zz_xzz = cbuffer.data(df_geom_20_off + 175 * ccomps * dcomps);

            auto g_xz_0_zz_yyy = cbuffer.data(df_geom_20_off + 176 * ccomps * dcomps);

            auto g_xz_0_zz_yyz = cbuffer.data(df_geom_20_off + 177 * ccomps * dcomps);

            auto g_xz_0_zz_yzz = cbuffer.data(df_geom_20_off + 178 * ccomps * dcomps);

            auto g_xz_0_zz_zzz = cbuffer.data(df_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxx, g_x_0_z_xxy, g_x_0_z_xxz, g_x_0_z_xyy, g_x_0_z_xyz, g_x_0_z_xzz, g_x_0_z_yyy, g_x_0_z_yyz, g_x_0_z_yzz, g_x_0_z_zzz, g_xz_0_z_xxx, g_xz_0_z_xxxz, g_xz_0_z_xxy, g_xz_0_z_xxyz, g_xz_0_z_xxz, g_xz_0_z_xxzz, g_xz_0_z_xyy, g_xz_0_z_xyyz, g_xz_0_z_xyz, g_xz_0_z_xyzz, g_xz_0_z_xzz, g_xz_0_z_xzzz, g_xz_0_z_yyy, g_xz_0_z_yyyz, g_xz_0_z_yyz, g_xz_0_z_yyzz, g_xz_0_z_yzz, g_xz_0_z_yzzz, g_xz_0_z_zzz, g_xz_0_z_zzzz, g_xz_0_zz_xxx, g_xz_0_zz_xxy, g_xz_0_zz_xxz, g_xz_0_zz_xyy, g_xz_0_zz_xyz, g_xz_0_zz_xzz, g_xz_0_zz_yyy, g_xz_0_zz_yyz, g_xz_0_zz_yzz, g_xz_0_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zz_xxx[k] = -g_x_0_z_xxx[k] - g_xz_0_z_xxx[k] * ab_z + g_xz_0_z_xxxz[k];

                g_xz_0_zz_xxy[k] = -g_x_0_z_xxy[k] - g_xz_0_z_xxy[k] * ab_z + g_xz_0_z_xxyz[k];

                g_xz_0_zz_xxz[k] = -g_x_0_z_xxz[k] - g_xz_0_z_xxz[k] * ab_z + g_xz_0_z_xxzz[k];

                g_xz_0_zz_xyy[k] = -g_x_0_z_xyy[k] - g_xz_0_z_xyy[k] * ab_z + g_xz_0_z_xyyz[k];

                g_xz_0_zz_xyz[k] = -g_x_0_z_xyz[k] - g_xz_0_z_xyz[k] * ab_z + g_xz_0_z_xyzz[k];

                g_xz_0_zz_xzz[k] = -g_x_0_z_xzz[k] - g_xz_0_z_xzz[k] * ab_z + g_xz_0_z_xzzz[k];

                g_xz_0_zz_yyy[k] = -g_x_0_z_yyy[k] - g_xz_0_z_yyy[k] * ab_z + g_xz_0_z_yyyz[k];

                g_xz_0_zz_yyz[k] = -g_x_0_z_yyz[k] - g_xz_0_z_yyz[k] * ab_z + g_xz_0_z_yyzz[k];

                g_xz_0_zz_yzz[k] = -g_x_0_z_yzz[k] - g_xz_0_z_yzz[k] * ab_z + g_xz_0_z_yzzz[k];

                g_xz_0_zz_zzz[k] = -g_x_0_z_zzz[k] - g_xz_0_z_zzz[k] * ab_z + g_xz_0_z_zzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xx_xxx = cbuffer.data(df_geom_20_off + 180 * ccomps * dcomps);

            auto g_yy_0_xx_xxy = cbuffer.data(df_geom_20_off + 181 * ccomps * dcomps);

            auto g_yy_0_xx_xxz = cbuffer.data(df_geom_20_off + 182 * ccomps * dcomps);

            auto g_yy_0_xx_xyy = cbuffer.data(df_geom_20_off + 183 * ccomps * dcomps);

            auto g_yy_0_xx_xyz = cbuffer.data(df_geom_20_off + 184 * ccomps * dcomps);

            auto g_yy_0_xx_xzz = cbuffer.data(df_geom_20_off + 185 * ccomps * dcomps);

            auto g_yy_0_xx_yyy = cbuffer.data(df_geom_20_off + 186 * ccomps * dcomps);

            auto g_yy_0_xx_yyz = cbuffer.data(df_geom_20_off + 187 * ccomps * dcomps);

            auto g_yy_0_xx_yzz = cbuffer.data(df_geom_20_off + 188 * ccomps * dcomps);

            auto g_yy_0_xx_zzz = cbuffer.data(df_geom_20_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_x_xxx, g_yy_0_x_xxxx, g_yy_0_x_xxxy, g_yy_0_x_xxxz, g_yy_0_x_xxy, g_yy_0_x_xxyy, g_yy_0_x_xxyz, g_yy_0_x_xxz, g_yy_0_x_xxzz, g_yy_0_x_xyy, g_yy_0_x_xyyy, g_yy_0_x_xyyz, g_yy_0_x_xyz, g_yy_0_x_xyzz, g_yy_0_x_xzz, g_yy_0_x_xzzz, g_yy_0_x_yyy, g_yy_0_x_yyz, g_yy_0_x_yzz, g_yy_0_x_zzz, g_yy_0_xx_xxx, g_yy_0_xx_xxy, g_yy_0_xx_xxz, g_yy_0_xx_xyy, g_yy_0_xx_xyz, g_yy_0_xx_xzz, g_yy_0_xx_yyy, g_yy_0_xx_yyz, g_yy_0_xx_yzz, g_yy_0_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xx_xxx[k] = -g_yy_0_x_xxx[k] * ab_x + g_yy_0_x_xxxx[k];

                g_yy_0_xx_xxy[k] = -g_yy_0_x_xxy[k] * ab_x + g_yy_0_x_xxxy[k];

                g_yy_0_xx_xxz[k] = -g_yy_0_x_xxz[k] * ab_x + g_yy_0_x_xxxz[k];

                g_yy_0_xx_xyy[k] = -g_yy_0_x_xyy[k] * ab_x + g_yy_0_x_xxyy[k];

                g_yy_0_xx_xyz[k] = -g_yy_0_x_xyz[k] * ab_x + g_yy_0_x_xxyz[k];

                g_yy_0_xx_xzz[k] = -g_yy_0_x_xzz[k] * ab_x + g_yy_0_x_xxzz[k];

                g_yy_0_xx_yyy[k] = -g_yy_0_x_yyy[k] * ab_x + g_yy_0_x_xyyy[k];

                g_yy_0_xx_yyz[k] = -g_yy_0_x_yyz[k] * ab_x + g_yy_0_x_xyyz[k];

                g_yy_0_xx_yzz[k] = -g_yy_0_x_yzz[k] * ab_x + g_yy_0_x_xyzz[k];

                g_yy_0_xx_zzz[k] = -g_yy_0_x_zzz[k] * ab_x + g_yy_0_x_xzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xy_xxx = cbuffer.data(df_geom_20_off + 190 * ccomps * dcomps);

            auto g_yy_0_xy_xxy = cbuffer.data(df_geom_20_off + 191 * ccomps * dcomps);

            auto g_yy_0_xy_xxz = cbuffer.data(df_geom_20_off + 192 * ccomps * dcomps);

            auto g_yy_0_xy_xyy = cbuffer.data(df_geom_20_off + 193 * ccomps * dcomps);

            auto g_yy_0_xy_xyz = cbuffer.data(df_geom_20_off + 194 * ccomps * dcomps);

            auto g_yy_0_xy_xzz = cbuffer.data(df_geom_20_off + 195 * ccomps * dcomps);

            auto g_yy_0_xy_yyy = cbuffer.data(df_geom_20_off + 196 * ccomps * dcomps);

            auto g_yy_0_xy_yyz = cbuffer.data(df_geom_20_off + 197 * ccomps * dcomps);

            auto g_yy_0_xy_yzz = cbuffer.data(df_geom_20_off + 198 * ccomps * dcomps);

            auto g_yy_0_xy_zzz = cbuffer.data(df_geom_20_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xy_xxx, g_yy_0_xy_xxy, g_yy_0_xy_xxz, g_yy_0_xy_xyy, g_yy_0_xy_xyz, g_yy_0_xy_xzz, g_yy_0_xy_yyy, g_yy_0_xy_yyz, g_yy_0_xy_yzz, g_yy_0_xy_zzz, g_yy_0_y_xxx, g_yy_0_y_xxxx, g_yy_0_y_xxxy, g_yy_0_y_xxxz, g_yy_0_y_xxy, g_yy_0_y_xxyy, g_yy_0_y_xxyz, g_yy_0_y_xxz, g_yy_0_y_xxzz, g_yy_0_y_xyy, g_yy_0_y_xyyy, g_yy_0_y_xyyz, g_yy_0_y_xyz, g_yy_0_y_xyzz, g_yy_0_y_xzz, g_yy_0_y_xzzz, g_yy_0_y_yyy, g_yy_0_y_yyz, g_yy_0_y_yzz, g_yy_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xy_xxx[k] = -g_yy_0_y_xxx[k] * ab_x + g_yy_0_y_xxxx[k];

                g_yy_0_xy_xxy[k] = -g_yy_0_y_xxy[k] * ab_x + g_yy_0_y_xxxy[k];

                g_yy_0_xy_xxz[k] = -g_yy_0_y_xxz[k] * ab_x + g_yy_0_y_xxxz[k];

                g_yy_0_xy_xyy[k] = -g_yy_0_y_xyy[k] * ab_x + g_yy_0_y_xxyy[k];

                g_yy_0_xy_xyz[k] = -g_yy_0_y_xyz[k] * ab_x + g_yy_0_y_xxyz[k];

                g_yy_0_xy_xzz[k] = -g_yy_0_y_xzz[k] * ab_x + g_yy_0_y_xxzz[k];

                g_yy_0_xy_yyy[k] = -g_yy_0_y_yyy[k] * ab_x + g_yy_0_y_xyyy[k];

                g_yy_0_xy_yyz[k] = -g_yy_0_y_yyz[k] * ab_x + g_yy_0_y_xyyz[k];

                g_yy_0_xy_yzz[k] = -g_yy_0_y_yzz[k] * ab_x + g_yy_0_y_xyzz[k];

                g_yy_0_xy_zzz[k] = -g_yy_0_y_zzz[k] * ab_x + g_yy_0_y_xzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xz_xxx = cbuffer.data(df_geom_20_off + 200 * ccomps * dcomps);

            auto g_yy_0_xz_xxy = cbuffer.data(df_geom_20_off + 201 * ccomps * dcomps);

            auto g_yy_0_xz_xxz = cbuffer.data(df_geom_20_off + 202 * ccomps * dcomps);

            auto g_yy_0_xz_xyy = cbuffer.data(df_geom_20_off + 203 * ccomps * dcomps);

            auto g_yy_0_xz_xyz = cbuffer.data(df_geom_20_off + 204 * ccomps * dcomps);

            auto g_yy_0_xz_xzz = cbuffer.data(df_geom_20_off + 205 * ccomps * dcomps);

            auto g_yy_0_xz_yyy = cbuffer.data(df_geom_20_off + 206 * ccomps * dcomps);

            auto g_yy_0_xz_yyz = cbuffer.data(df_geom_20_off + 207 * ccomps * dcomps);

            auto g_yy_0_xz_yzz = cbuffer.data(df_geom_20_off + 208 * ccomps * dcomps);

            auto g_yy_0_xz_zzz = cbuffer.data(df_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xz_xxx, g_yy_0_xz_xxy, g_yy_0_xz_xxz, g_yy_0_xz_xyy, g_yy_0_xz_xyz, g_yy_0_xz_xzz, g_yy_0_xz_yyy, g_yy_0_xz_yyz, g_yy_0_xz_yzz, g_yy_0_xz_zzz, g_yy_0_z_xxx, g_yy_0_z_xxxx, g_yy_0_z_xxxy, g_yy_0_z_xxxz, g_yy_0_z_xxy, g_yy_0_z_xxyy, g_yy_0_z_xxyz, g_yy_0_z_xxz, g_yy_0_z_xxzz, g_yy_0_z_xyy, g_yy_0_z_xyyy, g_yy_0_z_xyyz, g_yy_0_z_xyz, g_yy_0_z_xyzz, g_yy_0_z_xzz, g_yy_0_z_xzzz, g_yy_0_z_yyy, g_yy_0_z_yyz, g_yy_0_z_yzz, g_yy_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xz_xxx[k] = -g_yy_0_z_xxx[k] * ab_x + g_yy_0_z_xxxx[k];

                g_yy_0_xz_xxy[k] = -g_yy_0_z_xxy[k] * ab_x + g_yy_0_z_xxxy[k];

                g_yy_0_xz_xxz[k] = -g_yy_0_z_xxz[k] * ab_x + g_yy_0_z_xxxz[k];

                g_yy_0_xz_xyy[k] = -g_yy_0_z_xyy[k] * ab_x + g_yy_0_z_xxyy[k];

                g_yy_0_xz_xyz[k] = -g_yy_0_z_xyz[k] * ab_x + g_yy_0_z_xxyz[k];

                g_yy_0_xz_xzz[k] = -g_yy_0_z_xzz[k] * ab_x + g_yy_0_z_xxzz[k];

                g_yy_0_xz_yyy[k] = -g_yy_0_z_yyy[k] * ab_x + g_yy_0_z_xyyy[k];

                g_yy_0_xz_yyz[k] = -g_yy_0_z_yyz[k] * ab_x + g_yy_0_z_xyyz[k];

                g_yy_0_xz_yzz[k] = -g_yy_0_z_yzz[k] * ab_x + g_yy_0_z_xyzz[k];

                g_yy_0_xz_zzz[k] = -g_yy_0_z_zzz[k] * ab_x + g_yy_0_z_xzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yy_xxx = cbuffer.data(df_geom_20_off + 210 * ccomps * dcomps);

            auto g_yy_0_yy_xxy = cbuffer.data(df_geom_20_off + 211 * ccomps * dcomps);

            auto g_yy_0_yy_xxz = cbuffer.data(df_geom_20_off + 212 * ccomps * dcomps);

            auto g_yy_0_yy_xyy = cbuffer.data(df_geom_20_off + 213 * ccomps * dcomps);

            auto g_yy_0_yy_xyz = cbuffer.data(df_geom_20_off + 214 * ccomps * dcomps);

            auto g_yy_0_yy_xzz = cbuffer.data(df_geom_20_off + 215 * ccomps * dcomps);

            auto g_yy_0_yy_yyy = cbuffer.data(df_geom_20_off + 216 * ccomps * dcomps);

            auto g_yy_0_yy_yyz = cbuffer.data(df_geom_20_off + 217 * ccomps * dcomps);

            auto g_yy_0_yy_yzz = cbuffer.data(df_geom_20_off + 218 * ccomps * dcomps);

            auto g_yy_0_yy_zzz = cbuffer.data(df_geom_20_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxx, g_y_0_y_xxy, g_y_0_y_xxz, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xzz, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yzz, g_y_0_y_zzz, g_yy_0_y_xxx, g_yy_0_y_xxxy, g_yy_0_y_xxy, g_yy_0_y_xxyy, g_yy_0_y_xxyz, g_yy_0_y_xxz, g_yy_0_y_xyy, g_yy_0_y_xyyy, g_yy_0_y_xyyz, g_yy_0_y_xyz, g_yy_0_y_xyzz, g_yy_0_y_xzz, g_yy_0_y_yyy, g_yy_0_y_yyyy, g_yy_0_y_yyyz, g_yy_0_y_yyz, g_yy_0_y_yyzz, g_yy_0_y_yzz, g_yy_0_y_yzzz, g_yy_0_y_zzz, g_yy_0_yy_xxx, g_yy_0_yy_xxy, g_yy_0_yy_xxz, g_yy_0_yy_xyy, g_yy_0_yy_xyz, g_yy_0_yy_xzz, g_yy_0_yy_yyy, g_yy_0_yy_yyz, g_yy_0_yy_yzz, g_yy_0_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yy_xxx[k] = -2.0 * g_y_0_y_xxx[k] - g_yy_0_y_xxx[k] * ab_y + g_yy_0_y_xxxy[k];

                g_yy_0_yy_xxy[k] = -2.0 * g_y_0_y_xxy[k] - g_yy_0_y_xxy[k] * ab_y + g_yy_0_y_xxyy[k];

                g_yy_0_yy_xxz[k] = -2.0 * g_y_0_y_xxz[k] - g_yy_0_y_xxz[k] * ab_y + g_yy_0_y_xxyz[k];

                g_yy_0_yy_xyy[k] = -2.0 * g_y_0_y_xyy[k] - g_yy_0_y_xyy[k] * ab_y + g_yy_0_y_xyyy[k];

                g_yy_0_yy_xyz[k] = -2.0 * g_y_0_y_xyz[k] - g_yy_0_y_xyz[k] * ab_y + g_yy_0_y_xyyz[k];

                g_yy_0_yy_xzz[k] = -2.0 * g_y_0_y_xzz[k] - g_yy_0_y_xzz[k] * ab_y + g_yy_0_y_xyzz[k];

                g_yy_0_yy_yyy[k] = -2.0 * g_y_0_y_yyy[k] - g_yy_0_y_yyy[k] * ab_y + g_yy_0_y_yyyy[k];

                g_yy_0_yy_yyz[k] = -2.0 * g_y_0_y_yyz[k] - g_yy_0_y_yyz[k] * ab_y + g_yy_0_y_yyyz[k];

                g_yy_0_yy_yzz[k] = -2.0 * g_y_0_y_yzz[k] - g_yy_0_y_yzz[k] * ab_y + g_yy_0_y_yyzz[k];

                g_yy_0_yy_zzz[k] = -2.0 * g_y_0_y_zzz[k] - g_yy_0_y_zzz[k] * ab_y + g_yy_0_y_yzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yz_xxx = cbuffer.data(df_geom_20_off + 220 * ccomps * dcomps);

            auto g_yy_0_yz_xxy = cbuffer.data(df_geom_20_off + 221 * ccomps * dcomps);

            auto g_yy_0_yz_xxz = cbuffer.data(df_geom_20_off + 222 * ccomps * dcomps);

            auto g_yy_0_yz_xyy = cbuffer.data(df_geom_20_off + 223 * ccomps * dcomps);

            auto g_yy_0_yz_xyz = cbuffer.data(df_geom_20_off + 224 * ccomps * dcomps);

            auto g_yy_0_yz_xzz = cbuffer.data(df_geom_20_off + 225 * ccomps * dcomps);

            auto g_yy_0_yz_yyy = cbuffer.data(df_geom_20_off + 226 * ccomps * dcomps);

            auto g_yy_0_yz_yyz = cbuffer.data(df_geom_20_off + 227 * ccomps * dcomps);

            auto g_yy_0_yz_yzz = cbuffer.data(df_geom_20_off + 228 * ccomps * dcomps);

            auto g_yy_0_yz_zzz = cbuffer.data(df_geom_20_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_y_xxx, g_yy_0_y_xxxz, g_yy_0_y_xxy, g_yy_0_y_xxyz, g_yy_0_y_xxz, g_yy_0_y_xxzz, g_yy_0_y_xyy, g_yy_0_y_xyyz, g_yy_0_y_xyz, g_yy_0_y_xyzz, g_yy_0_y_xzz, g_yy_0_y_xzzz, g_yy_0_y_yyy, g_yy_0_y_yyyz, g_yy_0_y_yyz, g_yy_0_y_yyzz, g_yy_0_y_yzz, g_yy_0_y_yzzz, g_yy_0_y_zzz, g_yy_0_y_zzzz, g_yy_0_yz_xxx, g_yy_0_yz_xxy, g_yy_0_yz_xxz, g_yy_0_yz_xyy, g_yy_0_yz_xyz, g_yy_0_yz_xzz, g_yy_0_yz_yyy, g_yy_0_yz_yyz, g_yy_0_yz_yzz, g_yy_0_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yz_xxx[k] = -g_yy_0_y_xxx[k] * ab_z + g_yy_0_y_xxxz[k];

                g_yy_0_yz_xxy[k] = -g_yy_0_y_xxy[k] * ab_z + g_yy_0_y_xxyz[k];

                g_yy_0_yz_xxz[k] = -g_yy_0_y_xxz[k] * ab_z + g_yy_0_y_xxzz[k];

                g_yy_0_yz_xyy[k] = -g_yy_0_y_xyy[k] * ab_z + g_yy_0_y_xyyz[k];

                g_yy_0_yz_xyz[k] = -g_yy_0_y_xyz[k] * ab_z + g_yy_0_y_xyzz[k];

                g_yy_0_yz_xzz[k] = -g_yy_0_y_xzz[k] * ab_z + g_yy_0_y_xzzz[k];

                g_yy_0_yz_yyy[k] = -g_yy_0_y_yyy[k] * ab_z + g_yy_0_y_yyyz[k];

                g_yy_0_yz_yyz[k] = -g_yy_0_y_yyz[k] * ab_z + g_yy_0_y_yyzz[k];

                g_yy_0_yz_yzz[k] = -g_yy_0_y_yzz[k] * ab_z + g_yy_0_y_yzzz[k];

                g_yy_0_yz_zzz[k] = -g_yy_0_y_zzz[k] * ab_z + g_yy_0_y_zzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zz_xxx = cbuffer.data(df_geom_20_off + 230 * ccomps * dcomps);

            auto g_yy_0_zz_xxy = cbuffer.data(df_geom_20_off + 231 * ccomps * dcomps);

            auto g_yy_0_zz_xxz = cbuffer.data(df_geom_20_off + 232 * ccomps * dcomps);

            auto g_yy_0_zz_xyy = cbuffer.data(df_geom_20_off + 233 * ccomps * dcomps);

            auto g_yy_0_zz_xyz = cbuffer.data(df_geom_20_off + 234 * ccomps * dcomps);

            auto g_yy_0_zz_xzz = cbuffer.data(df_geom_20_off + 235 * ccomps * dcomps);

            auto g_yy_0_zz_yyy = cbuffer.data(df_geom_20_off + 236 * ccomps * dcomps);

            auto g_yy_0_zz_yyz = cbuffer.data(df_geom_20_off + 237 * ccomps * dcomps);

            auto g_yy_0_zz_yzz = cbuffer.data(df_geom_20_off + 238 * ccomps * dcomps);

            auto g_yy_0_zz_zzz = cbuffer.data(df_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_z_xxx, g_yy_0_z_xxxz, g_yy_0_z_xxy, g_yy_0_z_xxyz, g_yy_0_z_xxz, g_yy_0_z_xxzz, g_yy_0_z_xyy, g_yy_0_z_xyyz, g_yy_0_z_xyz, g_yy_0_z_xyzz, g_yy_0_z_xzz, g_yy_0_z_xzzz, g_yy_0_z_yyy, g_yy_0_z_yyyz, g_yy_0_z_yyz, g_yy_0_z_yyzz, g_yy_0_z_yzz, g_yy_0_z_yzzz, g_yy_0_z_zzz, g_yy_0_z_zzzz, g_yy_0_zz_xxx, g_yy_0_zz_xxy, g_yy_0_zz_xxz, g_yy_0_zz_xyy, g_yy_0_zz_xyz, g_yy_0_zz_xzz, g_yy_0_zz_yyy, g_yy_0_zz_yyz, g_yy_0_zz_yzz, g_yy_0_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zz_xxx[k] = -g_yy_0_z_xxx[k] * ab_z + g_yy_0_z_xxxz[k];

                g_yy_0_zz_xxy[k] = -g_yy_0_z_xxy[k] * ab_z + g_yy_0_z_xxyz[k];

                g_yy_0_zz_xxz[k] = -g_yy_0_z_xxz[k] * ab_z + g_yy_0_z_xxzz[k];

                g_yy_0_zz_xyy[k] = -g_yy_0_z_xyy[k] * ab_z + g_yy_0_z_xyyz[k];

                g_yy_0_zz_xyz[k] = -g_yy_0_z_xyz[k] * ab_z + g_yy_0_z_xyzz[k];

                g_yy_0_zz_xzz[k] = -g_yy_0_z_xzz[k] * ab_z + g_yy_0_z_xzzz[k];

                g_yy_0_zz_yyy[k] = -g_yy_0_z_yyy[k] * ab_z + g_yy_0_z_yyyz[k];

                g_yy_0_zz_yyz[k] = -g_yy_0_z_yyz[k] * ab_z + g_yy_0_z_yyzz[k];

                g_yy_0_zz_yzz[k] = -g_yy_0_z_yzz[k] * ab_z + g_yy_0_z_yzzz[k];

                g_yy_0_zz_zzz[k] = -g_yy_0_z_zzz[k] * ab_z + g_yy_0_z_zzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xx_xxx = cbuffer.data(df_geom_20_off + 240 * ccomps * dcomps);

            auto g_yz_0_xx_xxy = cbuffer.data(df_geom_20_off + 241 * ccomps * dcomps);

            auto g_yz_0_xx_xxz = cbuffer.data(df_geom_20_off + 242 * ccomps * dcomps);

            auto g_yz_0_xx_xyy = cbuffer.data(df_geom_20_off + 243 * ccomps * dcomps);

            auto g_yz_0_xx_xyz = cbuffer.data(df_geom_20_off + 244 * ccomps * dcomps);

            auto g_yz_0_xx_xzz = cbuffer.data(df_geom_20_off + 245 * ccomps * dcomps);

            auto g_yz_0_xx_yyy = cbuffer.data(df_geom_20_off + 246 * ccomps * dcomps);

            auto g_yz_0_xx_yyz = cbuffer.data(df_geom_20_off + 247 * ccomps * dcomps);

            auto g_yz_0_xx_yzz = cbuffer.data(df_geom_20_off + 248 * ccomps * dcomps);

            auto g_yz_0_xx_zzz = cbuffer.data(df_geom_20_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_x_xxx, g_yz_0_x_xxxx, g_yz_0_x_xxxy, g_yz_0_x_xxxz, g_yz_0_x_xxy, g_yz_0_x_xxyy, g_yz_0_x_xxyz, g_yz_0_x_xxz, g_yz_0_x_xxzz, g_yz_0_x_xyy, g_yz_0_x_xyyy, g_yz_0_x_xyyz, g_yz_0_x_xyz, g_yz_0_x_xyzz, g_yz_0_x_xzz, g_yz_0_x_xzzz, g_yz_0_x_yyy, g_yz_0_x_yyz, g_yz_0_x_yzz, g_yz_0_x_zzz, g_yz_0_xx_xxx, g_yz_0_xx_xxy, g_yz_0_xx_xxz, g_yz_0_xx_xyy, g_yz_0_xx_xyz, g_yz_0_xx_xzz, g_yz_0_xx_yyy, g_yz_0_xx_yyz, g_yz_0_xx_yzz, g_yz_0_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xx_xxx[k] = -g_yz_0_x_xxx[k] * ab_x + g_yz_0_x_xxxx[k];

                g_yz_0_xx_xxy[k] = -g_yz_0_x_xxy[k] * ab_x + g_yz_0_x_xxxy[k];

                g_yz_0_xx_xxz[k] = -g_yz_0_x_xxz[k] * ab_x + g_yz_0_x_xxxz[k];

                g_yz_0_xx_xyy[k] = -g_yz_0_x_xyy[k] * ab_x + g_yz_0_x_xxyy[k];

                g_yz_0_xx_xyz[k] = -g_yz_0_x_xyz[k] * ab_x + g_yz_0_x_xxyz[k];

                g_yz_0_xx_xzz[k] = -g_yz_0_x_xzz[k] * ab_x + g_yz_0_x_xxzz[k];

                g_yz_0_xx_yyy[k] = -g_yz_0_x_yyy[k] * ab_x + g_yz_0_x_xyyy[k];

                g_yz_0_xx_yyz[k] = -g_yz_0_x_yyz[k] * ab_x + g_yz_0_x_xyyz[k];

                g_yz_0_xx_yzz[k] = -g_yz_0_x_yzz[k] * ab_x + g_yz_0_x_xyzz[k];

                g_yz_0_xx_zzz[k] = -g_yz_0_x_zzz[k] * ab_x + g_yz_0_x_xzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xy_xxx = cbuffer.data(df_geom_20_off + 250 * ccomps * dcomps);

            auto g_yz_0_xy_xxy = cbuffer.data(df_geom_20_off + 251 * ccomps * dcomps);

            auto g_yz_0_xy_xxz = cbuffer.data(df_geom_20_off + 252 * ccomps * dcomps);

            auto g_yz_0_xy_xyy = cbuffer.data(df_geom_20_off + 253 * ccomps * dcomps);

            auto g_yz_0_xy_xyz = cbuffer.data(df_geom_20_off + 254 * ccomps * dcomps);

            auto g_yz_0_xy_xzz = cbuffer.data(df_geom_20_off + 255 * ccomps * dcomps);

            auto g_yz_0_xy_yyy = cbuffer.data(df_geom_20_off + 256 * ccomps * dcomps);

            auto g_yz_0_xy_yyz = cbuffer.data(df_geom_20_off + 257 * ccomps * dcomps);

            auto g_yz_0_xy_yzz = cbuffer.data(df_geom_20_off + 258 * ccomps * dcomps);

            auto g_yz_0_xy_zzz = cbuffer.data(df_geom_20_off + 259 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xy_xxx, g_yz_0_xy_xxy, g_yz_0_xy_xxz, g_yz_0_xy_xyy, g_yz_0_xy_xyz, g_yz_0_xy_xzz, g_yz_0_xy_yyy, g_yz_0_xy_yyz, g_yz_0_xy_yzz, g_yz_0_xy_zzz, g_yz_0_y_xxx, g_yz_0_y_xxxx, g_yz_0_y_xxxy, g_yz_0_y_xxxz, g_yz_0_y_xxy, g_yz_0_y_xxyy, g_yz_0_y_xxyz, g_yz_0_y_xxz, g_yz_0_y_xxzz, g_yz_0_y_xyy, g_yz_0_y_xyyy, g_yz_0_y_xyyz, g_yz_0_y_xyz, g_yz_0_y_xyzz, g_yz_0_y_xzz, g_yz_0_y_xzzz, g_yz_0_y_yyy, g_yz_0_y_yyz, g_yz_0_y_yzz, g_yz_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xy_xxx[k] = -g_yz_0_y_xxx[k] * ab_x + g_yz_0_y_xxxx[k];

                g_yz_0_xy_xxy[k] = -g_yz_0_y_xxy[k] * ab_x + g_yz_0_y_xxxy[k];

                g_yz_0_xy_xxz[k] = -g_yz_0_y_xxz[k] * ab_x + g_yz_0_y_xxxz[k];

                g_yz_0_xy_xyy[k] = -g_yz_0_y_xyy[k] * ab_x + g_yz_0_y_xxyy[k];

                g_yz_0_xy_xyz[k] = -g_yz_0_y_xyz[k] * ab_x + g_yz_0_y_xxyz[k];

                g_yz_0_xy_xzz[k] = -g_yz_0_y_xzz[k] * ab_x + g_yz_0_y_xxzz[k];

                g_yz_0_xy_yyy[k] = -g_yz_0_y_yyy[k] * ab_x + g_yz_0_y_xyyy[k];

                g_yz_0_xy_yyz[k] = -g_yz_0_y_yyz[k] * ab_x + g_yz_0_y_xyyz[k];

                g_yz_0_xy_yzz[k] = -g_yz_0_y_yzz[k] * ab_x + g_yz_0_y_xyzz[k];

                g_yz_0_xy_zzz[k] = -g_yz_0_y_zzz[k] * ab_x + g_yz_0_y_xzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xz_xxx = cbuffer.data(df_geom_20_off + 260 * ccomps * dcomps);

            auto g_yz_0_xz_xxy = cbuffer.data(df_geom_20_off + 261 * ccomps * dcomps);

            auto g_yz_0_xz_xxz = cbuffer.data(df_geom_20_off + 262 * ccomps * dcomps);

            auto g_yz_0_xz_xyy = cbuffer.data(df_geom_20_off + 263 * ccomps * dcomps);

            auto g_yz_0_xz_xyz = cbuffer.data(df_geom_20_off + 264 * ccomps * dcomps);

            auto g_yz_0_xz_xzz = cbuffer.data(df_geom_20_off + 265 * ccomps * dcomps);

            auto g_yz_0_xz_yyy = cbuffer.data(df_geom_20_off + 266 * ccomps * dcomps);

            auto g_yz_0_xz_yyz = cbuffer.data(df_geom_20_off + 267 * ccomps * dcomps);

            auto g_yz_0_xz_yzz = cbuffer.data(df_geom_20_off + 268 * ccomps * dcomps);

            auto g_yz_0_xz_zzz = cbuffer.data(df_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xz_xxx, g_yz_0_xz_xxy, g_yz_0_xz_xxz, g_yz_0_xz_xyy, g_yz_0_xz_xyz, g_yz_0_xz_xzz, g_yz_0_xz_yyy, g_yz_0_xz_yyz, g_yz_0_xz_yzz, g_yz_0_xz_zzz, g_yz_0_z_xxx, g_yz_0_z_xxxx, g_yz_0_z_xxxy, g_yz_0_z_xxxz, g_yz_0_z_xxy, g_yz_0_z_xxyy, g_yz_0_z_xxyz, g_yz_0_z_xxz, g_yz_0_z_xxzz, g_yz_0_z_xyy, g_yz_0_z_xyyy, g_yz_0_z_xyyz, g_yz_0_z_xyz, g_yz_0_z_xyzz, g_yz_0_z_xzz, g_yz_0_z_xzzz, g_yz_0_z_yyy, g_yz_0_z_yyz, g_yz_0_z_yzz, g_yz_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xz_xxx[k] = -g_yz_0_z_xxx[k] * ab_x + g_yz_0_z_xxxx[k];

                g_yz_0_xz_xxy[k] = -g_yz_0_z_xxy[k] * ab_x + g_yz_0_z_xxxy[k];

                g_yz_0_xz_xxz[k] = -g_yz_0_z_xxz[k] * ab_x + g_yz_0_z_xxxz[k];

                g_yz_0_xz_xyy[k] = -g_yz_0_z_xyy[k] * ab_x + g_yz_0_z_xxyy[k];

                g_yz_0_xz_xyz[k] = -g_yz_0_z_xyz[k] * ab_x + g_yz_0_z_xxyz[k];

                g_yz_0_xz_xzz[k] = -g_yz_0_z_xzz[k] * ab_x + g_yz_0_z_xxzz[k];

                g_yz_0_xz_yyy[k] = -g_yz_0_z_yyy[k] * ab_x + g_yz_0_z_xyyy[k];

                g_yz_0_xz_yyz[k] = -g_yz_0_z_yyz[k] * ab_x + g_yz_0_z_xyyz[k];

                g_yz_0_xz_yzz[k] = -g_yz_0_z_yzz[k] * ab_x + g_yz_0_z_xyzz[k];

                g_yz_0_xz_zzz[k] = -g_yz_0_z_zzz[k] * ab_x + g_yz_0_z_xzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yy_xxx = cbuffer.data(df_geom_20_off + 270 * ccomps * dcomps);

            auto g_yz_0_yy_xxy = cbuffer.data(df_geom_20_off + 271 * ccomps * dcomps);

            auto g_yz_0_yy_xxz = cbuffer.data(df_geom_20_off + 272 * ccomps * dcomps);

            auto g_yz_0_yy_xyy = cbuffer.data(df_geom_20_off + 273 * ccomps * dcomps);

            auto g_yz_0_yy_xyz = cbuffer.data(df_geom_20_off + 274 * ccomps * dcomps);

            auto g_yz_0_yy_xzz = cbuffer.data(df_geom_20_off + 275 * ccomps * dcomps);

            auto g_yz_0_yy_yyy = cbuffer.data(df_geom_20_off + 276 * ccomps * dcomps);

            auto g_yz_0_yy_yyz = cbuffer.data(df_geom_20_off + 277 * ccomps * dcomps);

            auto g_yz_0_yy_yzz = cbuffer.data(df_geom_20_off + 278 * ccomps * dcomps);

            auto g_yz_0_yy_zzz = cbuffer.data(df_geom_20_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_y_xxx, g_yz_0_y_xxxy, g_yz_0_y_xxy, g_yz_0_y_xxyy, g_yz_0_y_xxyz, g_yz_0_y_xxz, g_yz_0_y_xyy, g_yz_0_y_xyyy, g_yz_0_y_xyyz, g_yz_0_y_xyz, g_yz_0_y_xyzz, g_yz_0_y_xzz, g_yz_0_y_yyy, g_yz_0_y_yyyy, g_yz_0_y_yyyz, g_yz_0_y_yyz, g_yz_0_y_yyzz, g_yz_0_y_yzz, g_yz_0_y_yzzz, g_yz_0_y_zzz, g_yz_0_yy_xxx, g_yz_0_yy_xxy, g_yz_0_yy_xxz, g_yz_0_yy_xyy, g_yz_0_yy_xyz, g_yz_0_yy_xzz, g_yz_0_yy_yyy, g_yz_0_yy_yyz, g_yz_0_yy_yzz, g_yz_0_yy_zzz, g_z_0_y_xxx, g_z_0_y_xxy, g_z_0_y_xxz, g_z_0_y_xyy, g_z_0_y_xyz, g_z_0_y_xzz, g_z_0_y_yyy, g_z_0_y_yyz, g_z_0_y_yzz, g_z_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yy_xxx[k] = -g_z_0_y_xxx[k] - g_yz_0_y_xxx[k] * ab_y + g_yz_0_y_xxxy[k];

                g_yz_0_yy_xxy[k] = -g_z_0_y_xxy[k] - g_yz_0_y_xxy[k] * ab_y + g_yz_0_y_xxyy[k];

                g_yz_0_yy_xxz[k] = -g_z_0_y_xxz[k] - g_yz_0_y_xxz[k] * ab_y + g_yz_0_y_xxyz[k];

                g_yz_0_yy_xyy[k] = -g_z_0_y_xyy[k] - g_yz_0_y_xyy[k] * ab_y + g_yz_0_y_xyyy[k];

                g_yz_0_yy_xyz[k] = -g_z_0_y_xyz[k] - g_yz_0_y_xyz[k] * ab_y + g_yz_0_y_xyyz[k];

                g_yz_0_yy_xzz[k] = -g_z_0_y_xzz[k] - g_yz_0_y_xzz[k] * ab_y + g_yz_0_y_xyzz[k];

                g_yz_0_yy_yyy[k] = -g_z_0_y_yyy[k] - g_yz_0_y_yyy[k] * ab_y + g_yz_0_y_yyyy[k];

                g_yz_0_yy_yyz[k] = -g_z_0_y_yyz[k] - g_yz_0_y_yyz[k] * ab_y + g_yz_0_y_yyyz[k];

                g_yz_0_yy_yzz[k] = -g_z_0_y_yzz[k] - g_yz_0_y_yzz[k] * ab_y + g_yz_0_y_yyzz[k];

                g_yz_0_yy_zzz[k] = -g_z_0_y_zzz[k] - g_yz_0_y_zzz[k] * ab_y + g_yz_0_y_yzzz[k];
            }

            /// Set up 280-290 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yz_xxx = cbuffer.data(df_geom_20_off + 280 * ccomps * dcomps);

            auto g_yz_0_yz_xxy = cbuffer.data(df_geom_20_off + 281 * ccomps * dcomps);

            auto g_yz_0_yz_xxz = cbuffer.data(df_geom_20_off + 282 * ccomps * dcomps);

            auto g_yz_0_yz_xyy = cbuffer.data(df_geom_20_off + 283 * ccomps * dcomps);

            auto g_yz_0_yz_xyz = cbuffer.data(df_geom_20_off + 284 * ccomps * dcomps);

            auto g_yz_0_yz_xzz = cbuffer.data(df_geom_20_off + 285 * ccomps * dcomps);

            auto g_yz_0_yz_yyy = cbuffer.data(df_geom_20_off + 286 * ccomps * dcomps);

            auto g_yz_0_yz_yyz = cbuffer.data(df_geom_20_off + 287 * ccomps * dcomps);

            auto g_yz_0_yz_yzz = cbuffer.data(df_geom_20_off + 288 * ccomps * dcomps);

            auto g_yz_0_yz_zzz = cbuffer.data(df_geom_20_off + 289 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yz_xxx, g_yz_0_yz_xxy, g_yz_0_yz_xxz, g_yz_0_yz_xyy, g_yz_0_yz_xyz, g_yz_0_yz_xzz, g_yz_0_yz_yyy, g_yz_0_yz_yyz, g_yz_0_yz_yzz, g_yz_0_yz_zzz, g_yz_0_z_xxx, g_yz_0_z_xxxy, g_yz_0_z_xxy, g_yz_0_z_xxyy, g_yz_0_z_xxyz, g_yz_0_z_xxz, g_yz_0_z_xyy, g_yz_0_z_xyyy, g_yz_0_z_xyyz, g_yz_0_z_xyz, g_yz_0_z_xyzz, g_yz_0_z_xzz, g_yz_0_z_yyy, g_yz_0_z_yyyy, g_yz_0_z_yyyz, g_yz_0_z_yyz, g_yz_0_z_yyzz, g_yz_0_z_yzz, g_yz_0_z_yzzz, g_yz_0_z_zzz, g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yz_xxx[k] = -g_z_0_z_xxx[k] - g_yz_0_z_xxx[k] * ab_y + g_yz_0_z_xxxy[k];

                g_yz_0_yz_xxy[k] = -g_z_0_z_xxy[k] - g_yz_0_z_xxy[k] * ab_y + g_yz_0_z_xxyy[k];

                g_yz_0_yz_xxz[k] = -g_z_0_z_xxz[k] - g_yz_0_z_xxz[k] * ab_y + g_yz_0_z_xxyz[k];

                g_yz_0_yz_xyy[k] = -g_z_0_z_xyy[k] - g_yz_0_z_xyy[k] * ab_y + g_yz_0_z_xyyy[k];

                g_yz_0_yz_xyz[k] = -g_z_0_z_xyz[k] - g_yz_0_z_xyz[k] * ab_y + g_yz_0_z_xyyz[k];

                g_yz_0_yz_xzz[k] = -g_z_0_z_xzz[k] - g_yz_0_z_xzz[k] * ab_y + g_yz_0_z_xyzz[k];

                g_yz_0_yz_yyy[k] = -g_z_0_z_yyy[k] - g_yz_0_z_yyy[k] * ab_y + g_yz_0_z_yyyy[k];

                g_yz_0_yz_yyz[k] = -g_z_0_z_yyz[k] - g_yz_0_z_yyz[k] * ab_y + g_yz_0_z_yyyz[k];

                g_yz_0_yz_yzz[k] = -g_z_0_z_yzz[k] - g_yz_0_z_yzz[k] * ab_y + g_yz_0_z_yyzz[k];

                g_yz_0_yz_zzz[k] = -g_z_0_z_zzz[k] - g_yz_0_z_zzz[k] * ab_y + g_yz_0_z_yzzz[k];
            }

            /// Set up 290-300 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zz_xxx = cbuffer.data(df_geom_20_off + 290 * ccomps * dcomps);

            auto g_yz_0_zz_xxy = cbuffer.data(df_geom_20_off + 291 * ccomps * dcomps);

            auto g_yz_0_zz_xxz = cbuffer.data(df_geom_20_off + 292 * ccomps * dcomps);

            auto g_yz_0_zz_xyy = cbuffer.data(df_geom_20_off + 293 * ccomps * dcomps);

            auto g_yz_0_zz_xyz = cbuffer.data(df_geom_20_off + 294 * ccomps * dcomps);

            auto g_yz_0_zz_xzz = cbuffer.data(df_geom_20_off + 295 * ccomps * dcomps);

            auto g_yz_0_zz_yyy = cbuffer.data(df_geom_20_off + 296 * ccomps * dcomps);

            auto g_yz_0_zz_yyz = cbuffer.data(df_geom_20_off + 297 * ccomps * dcomps);

            auto g_yz_0_zz_yzz = cbuffer.data(df_geom_20_off + 298 * ccomps * dcomps);

            auto g_yz_0_zz_zzz = cbuffer.data(df_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxx, g_y_0_z_xxy, g_y_0_z_xxz, g_y_0_z_xyy, g_y_0_z_xyz, g_y_0_z_xzz, g_y_0_z_yyy, g_y_0_z_yyz, g_y_0_z_yzz, g_y_0_z_zzz, g_yz_0_z_xxx, g_yz_0_z_xxxz, g_yz_0_z_xxy, g_yz_0_z_xxyz, g_yz_0_z_xxz, g_yz_0_z_xxzz, g_yz_0_z_xyy, g_yz_0_z_xyyz, g_yz_0_z_xyz, g_yz_0_z_xyzz, g_yz_0_z_xzz, g_yz_0_z_xzzz, g_yz_0_z_yyy, g_yz_0_z_yyyz, g_yz_0_z_yyz, g_yz_0_z_yyzz, g_yz_0_z_yzz, g_yz_0_z_yzzz, g_yz_0_z_zzz, g_yz_0_z_zzzz, g_yz_0_zz_xxx, g_yz_0_zz_xxy, g_yz_0_zz_xxz, g_yz_0_zz_xyy, g_yz_0_zz_xyz, g_yz_0_zz_xzz, g_yz_0_zz_yyy, g_yz_0_zz_yyz, g_yz_0_zz_yzz, g_yz_0_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zz_xxx[k] = -g_y_0_z_xxx[k] - g_yz_0_z_xxx[k] * ab_z + g_yz_0_z_xxxz[k];

                g_yz_0_zz_xxy[k] = -g_y_0_z_xxy[k] - g_yz_0_z_xxy[k] * ab_z + g_yz_0_z_xxyz[k];

                g_yz_0_zz_xxz[k] = -g_y_0_z_xxz[k] - g_yz_0_z_xxz[k] * ab_z + g_yz_0_z_xxzz[k];

                g_yz_0_zz_xyy[k] = -g_y_0_z_xyy[k] - g_yz_0_z_xyy[k] * ab_z + g_yz_0_z_xyyz[k];

                g_yz_0_zz_xyz[k] = -g_y_0_z_xyz[k] - g_yz_0_z_xyz[k] * ab_z + g_yz_0_z_xyzz[k];

                g_yz_0_zz_xzz[k] = -g_y_0_z_xzz[k] - g_yz_0_z_xzz[k] * ab_z + g_yz_0_z_xzzz[k];

                g_yz_0_zz_yyy[k] = -g_y_0_z_yyy[k] - g_yz_0_z_yyy[k] * ab_z + g_yz_0_z_yyyz[k];

                g_yz_0_zz_yyz[k] = -g_y_0_z_yyz[k] - g_yz_0_z_yyz[k] * ab_z + g_yz_0_z_yyzz[k];

                g_yz_0_zz_yzz[k] = -g_y_0_z_yzz[k] - g_yz_0_z_yzz[k] * ab_z + g_yz_0_z_yzzz[k];

                g_yz_0_zz_zzz[k] = -g_y_0_z_zzz[k] - g_yz_0_z_zzz[k] * ab_z + g_yz_0_z_zzzz[k];
            }

            /// Set up 300-310 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xx_xxx = cbuffer.data(df_geom_20_off + 300 * ccomps * dcomps);

            auto g_zz_0_xx_xxy = cbuffer.data(df_geom_20_off + 301 * ccomps * dcomps);

            auto g_zz_0_xx_xxz = cbuffer.data(df_geom_20_off + 302 * ccomps * dcomps);

            auto g_zz_0_xx_xyy = cbuffer.data(df_geom_20_off + 303 * ccomps * dcomps);

            auto g_zz_0_xx_xyz = cbuffer.data(df_geom_20_off + 304 * ccomps * dcomps);

            auto g_zz_0_xx_xzz = cbuffer.data(df_geom_20_off + 305 * ccomps * dcomps);

            auto g_zz_0_xx_yyy = cbuffer.data(df_geom_20_off + 306 * ccomps * dcomps);

            auto g_zz_0_xx_yyz = cbuffer.data(df_geom_20_off + 307 * ccomps * dcomps);

            auto g_zz_0_xx_yzz = cbuffer.data(df_geom_20_off + 308 * ccomps * dcomps);

            auto g_zz_0_xx_zzz = cbuffer.data(df_geom_20_off + 309 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_x_xxx, g_zz_0_x_xxxx, g_zz_0_x_xxxy, g_zz_0_x_xxxz, g_zz_0_x_xxy, g_zz_0_x_xxyy, g_zz_0_x_xxyz, g_zz_0_x_xxz, g_zz_0_x_xxzz, g_zz_0_x_xyy, g_zz_0_x_xyyy, g_zz_0_x_xyyz, g_zz_0_x_xyz, g_zz_0_x_xyzz, g_zz_0_x_xzz, g_zz_0_x_xzzz, g_zz_0_x_yyy, g_zz_0_x_yyz, g_zz_0_x_yzz, g_zz_0_x_zzz, g_zz_0_xx_xxx, g_zz_0_xx_xxy, g_zz_0_xx_xxz, g_zz_0_xx_xyy, g_zz_0_xx_xyz, g_zz_0_xx_xzz, g_zz_0_xx_yyy, g_zz_0_xx_yyz, g_zz_0_xx_yzz, g_zz_0_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xx_xxx[k] = -g_zz_0_x_xxx[k] * ab_x + g_zz_0_x_xxxx[k];

                g_zz_0_xx_xxy[k] = -g_zz_0_x_xxy[k] * ab_x + g_zz_0_x_xxxy[k];

                g_zz_0_xx_xxz[k] = -g_zz_0_x_xxz[k] * ab_x + g_zz_0_x_xxxz[k];

                g_zz_0_xx_xyy[k] = -g_zz_0_x_xyy[k] * ab_x + g_zz_0_x_xxyy[k];

                g_zz_0_xx_xyz[k] = -g_zz_0_x_xyz[k] * ab_x + g_zz_0_x_xxyz[k];

                g_zz_0_xx_xzz[k] = -g_zz_0_x_xzz[k] * ab_x + g_zz_0_x_xxzz[k];

                g_zz_0_xx_yyy[k] = -g_zz_0_x_yyy[k] * ab_x + g_zz_0_x_xyyy[k];

                g_zz_0_xx_yyz[k] = -g_zz_0_x_yyz[k] * ab_x + g_zz_0_x_xyyz[k];

                g_zz_0_xx_yzz[k] = -g_zz_0_x_yzz[k] * ab_x + g_zz_0_x_xyzz[k];

                g_zz_0_xx_zzz[k] = -g_zz_0_x_zzz[k] * ab_x + g_zz_0_x_xzzz[k];
            }

            /// Set up 310-320 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xy_xxx = cbuffer.data(df_geom_20_off + 310 * ccomps * dcomps);

            auto g_zz_0_xy_xxy = cbuffer.data(df_geom_20_off + 311 * ccomps * dcomps);

            auto g_zz_0_xy_xxz = cbuffer.data(df_geom_20_off + 312 * ccomps * dcomps);

            auto g_zz_0_xy_xyy = cbuffer.data(df_geom_20_off + 313 * ccomps * dcomps);

            auto g_zz_0_xy_xyz = cbuffer.data(df_geom_20_off + 314 * ccomps * dcomps);

            auto g_zz_0_xy_xzz = cbuffer.data(df_geom_20_off + 315 * ccomps * dcomps);

            auto g_zz_0_xy_yyy = cbuffer.data(df_geom_20_off + 316 * ccomps * dcomps);

            auto g_zz_0_xy_yyz = cbuffer.data(df_geom_20_off + 317 * ccomps * dcomps);

            auto g_zz_0_xy_yzz = cbuffer.data(df_geom_20_off + 318 * ccomps * dcomps);

            auto g_zz_0_xy_zzz = cbuffer.data(df_geom_20_off + 319 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xy_xxx, g_zz_0_xy_xxy, g_zz_0_xy_xxz, g_zz_0_xy_xyy, g_zz_0_xy_xyz, g_zz_0_xy_xzz, g_zz_0_xy_yyy, g_zz_0_xy_yyz, g_zz_0_xy_yzz, g_zz_0_xy_zzz, g_zz_0_y_xxx, g_zz_0_y_xxxx, g_zz_0_y_xxxy, g_zz_0_y_xxxz, g_zz_0_y_xxy, g_zz_0_y_xxyy, g_zz_0_y_xxyz, g_zz_0_y_xxz, g_zz_0_y_xxzz, g_zz_0_y_xyy, g_zz_0_y_xyyy, g_zz_0_y_xyyz, g_zz_0_y_xyz, g_zz_0_y_xyzz, g_zz_0_y_xzz, g_zz_0_y_xzzz, g_zz_0_y_yyy, g_zz_0_y_yyz, g_zz_0_y_yzz, g_zz_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xy_xxx[k] = -g_zz_0_y_xxx[k] * ab_x + g_zz_0_y_xxxx[k];

                g_zz_0_xy_xxy[k] = -g_zz_0_y_xxy[k] * ab_x + g_zz_0_y_xxxy[k];

                g_zz_0_xy_xxz[k] = -g_zz_0_y_xxz[k] * ab_x + g_zz_0_y_xxxz[k];

                g_zz_0_xy_xyy[k] = -g_zz_0_y_xyy[k] * ab_x + g_zz_0_y_xxyy[k];

                g_zz_0_xy_xyz[k] = -g_zz_0_y_xyz[k] * ab_x + g_zz_0_y_xxyz[k];

                g_zz_0_xy_xzz[k] = -g_zz_0_y_xzz[k] * ab_x + g_zz_0_y_xxzz[k];

                g_zz_0_xy_yyy[k] = -g_zz_0_y_yyy[k] * ab_x + g_zz_0_y_xyyy[k];

                g_zz_0_xy_yyz[k] = -g_zz_0_y_yyz[k] * ab_x + g_zz_0_y_xyyz[k];

                g_zz_0_xy_yzz[k] = -g_zz_0_y_yzz[k] * ab_x + g_zz_0_y_xyzz[k];

                g_zz_0_xy_zzz[k] = -g_zz_0_y_zzz[k] * ab_x + g_zz_0_y_xzzz[k];
            }

            /// Set up 320-330 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xz_xxx = cbuffer.data(df_geom_20_off + 320 * ccomps * dcomps);

            auto g_zz_0_xz_xxy = cbuffer.data(df_geom_20_off + 321 * ccomps * dcomps);

            auto g_zz_0_xz_xxz = cbuffer.data(df_geom_20_off + 322 * ccomps * dcomps);

            auto g_zz_0_xz_xyy = cbuffer.data(df_geom_20_off + 323 * ccomps * dcomps);

            auto g_zz_0_xz_xyz = cbuffer.data(df_geom_20_off + 324 * ccomps * dcomps);

            auto g_zz_0_xz_xzz = cbuffer.data(df_geom_20_off + 325 * ccomps * dcomps);

            auto g_zz_0_xz_yyy = cbuffer.data(df_geom_20_off + 326 * ccomps * dcomps);

            auto g_zz_0_xz_yyz = cbuffer.data(df_geom_20_off + 327 * ccomps * dcomps);

            auto g_zz_0_xz_yzz = cbuffer.data(df_geom_20_off + 328 * ccomps * dcomps);

            auto g_zz_0_xz_zzz = cbuffer.data(df_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xz_xxx, g_zz_0_xz_xxy, g_zz_0_xz_xxz, g_zz_0_xz_xyy, g_zz_0_xz_xyz, g_zz_0_xz_xzz, g_zz_0_xz_yyy, g_zz_0_xz_yyz, g_zz_0_xz_yzz, g_zz_0_xz_zzz, g_zz_0_z_xxx, g_zz_0_z_xxxx, g_zz_0_z_xxxy, g_zz_0_z_xxxz, g_zz_0_z_xxy, g_zz_0_z_xxyy, g_zz_0_z_xxyz, g_zz_0_z_xxz, g_zz_0_z_xxzz, g_zz_0_z_xyy, g_zz_0_z_xyyy, g_zz_0_z_xyyz, g_zz_0_z_xyz, g_zz_0_z_xyzz, g_zz_0_z_xzz, g_zz_0_z_xzzz, g_zz_0_z_yyy, g_zz_0_z_yyz, g_zz_0_z_yzz, g_zz_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xz_xxx[k] = -g_zz_0_z_xxx[k] * ab_x + g_zz_0_z_xxxx[k];

                g_zz_0_xz_xxy[k] = -g_zz_0_z_xxy[k] * ab_x + g_zz_0_z_xxxy[k];

                g_zz_0_xz_xxz[k] = -g_zz_0_z_xxz[k] * ab_x + g_zz_0_z_xxxz[k];

                g_zz_0_xz_xyy[k] = -g_zz_0_z_xyy[k] * ab_x + g_zz_0_z_xxyy[k];

                g_zz_0_xz_xyz[k] = -g_zz_0_z_xyz[k] * ab_x + g_zz_0_z_xxyz[k];

                g_zz_0_xz_xzz[k] = -g_zz_0_z_xzz[k] * ab_x + g_zz_0_z_xxzz[k];

                g_zz_0_xz_yyy[k] = -g_zz_0_z_yyy[k] * ab_x + g_zz_0_z_xyyy[k];

                g_zz_0_xz_yyz[k] = -g_zz_0_z_yyz[k] * ab_x + g_zz_0_z_xyyz[k];

                g_zz_0_xz_yzz[k] = -g_zz_0_z_yzz[k] * ab_x + g_zz_0_z_xyzz[k];

                g_zz_0_xz_zzz[k] = -g_zz_0_z_zzz[k] * ab_x + g_zz_0_z_xzzz[k];
            }

            /// Set up 330-340 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yy_xxx = cbuffer.data(df_geom_20_off + 330 * ccomps * dcomps);

            auto g_zz_0_yy_xxy = cbuffer.data(df_geom_20_off + 331 * ccomps * dcomps);

            auto g_zz_0_yy_xxz = cbuffer.data(df_geom_20_off + 332 * ccomps * dcomps);

            auto g_zz_0_yy_xyy = cbuffer.data(df_geom_20_off + 333 * ccomps * dcomps);

            auto g_zz_0_yy_xyz = cbuffer.data(df_geom_20_off + 334 * ccomps * dcomps);

            auto g_zz_0_yy_xzz = cbuffer.data(df_geom_20_off + 335 * ccomps * dcomps);

            auto g_zz_0_yy_yyy = cbuffer.data(df_geom_20_off + 336 * ccomps * dcomps);

            auto g_zz_0_yy_yyz = cbuffer.data(df_geom_20_off + 337 * ccomps * dcomps);

            auto g_zz_0_yy_yzz = cbuffer.data(df_geom_20_off + 338 * ccomps * dcomps);

            auto g_zz_0_yy_zzz = cbuffer.data(df_geom_20_off + 339 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_y_xxx, g_zz_0_y_xxxy, g_zz_0_y_xxy, g_zz_0_y_xxyy, g_zz_0_y_xxyz, g_zz_0_y_xxz, g_zz_0_y_xyy, g_zz_0_y_xyyy, g_zz_0_y_xyyz, g_zz_0_y_xyz, g_zz_0_y_xyzz, g_zz_0_y_xzz, g_zz_0_y_yyy, g_zz_0_y_yyyy, g_zz_0_y_yyyz, g_zz_0_y_yyz, g_zz_0_y_yyzz, g_zz_0_y_yzz, g_zz_0_y_yzzz, g_zz_0_y_zzz, g_zz_0_yy_xxx, g_zz_0_yy_xxy, g_zz_0_yy_xxz, g_zz_0_yy_xyy, g_zz_0_yy_xyz, g_zz_0_yy_xzz, g_zz_0_yy_yyy, g_zz_0_yy_yyz, g_zz_0_yy_yzz, g_zz_0_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yy_xxx[k] = -g_zz_0_y_xxx[k] * ab_y + g_zz_0_y_xxxy[k];

                g_zz_0_yy_xxy[k] = -g_zz_0_y_xxy[k] * ab_y + g_zz_0_y_xxyy[k];

                g_zz_0_yy_xxz[k] = -g_zz_0_y_xxz[k] * ab_y + g_zz_0_y_xxyz[k];

                g_zz_0_yy_xyy[k] = -g_zz_0_y_xyy[k] * ab_y + g_zz_0_y_xyyy[k];

                g_zz_0_yy_xyz[k] = -g_zz_0_y_xyz[k] * ab_y + g_zz_0_y_xyyz[k];

                g_zz_0_yy_xzz[k] = -g_zz_0_y_xzz[k] * ab_y + g_zz_0_y_xyzz[k];

                g_zz_0_yy_yyy[k] = -g_zz_0_y_yyy[k] * ab_y + g_zz_0_y_yyyy[k];

                g_zz_0_yy_yyz[k] = -g_zz_0_y_yyz[k] * ab_y + g_zz_0_y_yyyz[k];

                g_zz_0_yy_yzz[k] = -g_zz_0_y_yzz[k] * ab_y + g_zz_0_y_yyzz[k];

                g_zz_0_yy_zzz[k] = -g_zz_0_y_zzz[k] * ab_y + g_zz_0_y_yzzz[k];
            }

            /// Set up 340-350 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yz_xxx = cbuffer.data(df_geom_20_off + 340 * ccomps * dcomps);

            auto g_zz_0_yz_xxy = cbuffer.data(df_geom_20_off + 341 * ccomps * dcomps);

            auto g_zz_0_yz_xxz = cbuffer.data(df_geom_20_off + 342 * ccomps * dcomps);

            auto g_zz_0_yz_xyy = cbuffer.data(df_geom_20_off + 343 * ccomps * dcomps);

            auto g_zz_0_yz_xyz = cbuffer.data(df_geom_20_off + 344 * ccomps * dcomps);

            auto g_zz_0_yz_xzz = cbuffer.data(df_geom_20_off + 345 * ccomps * dcomps);

            auto g_zz_0_yz_yyy = cbuffer.data(df_geom_20_off + 346 * ccomps * dcomps);

            auto g_zz_0_yz_yyz = cbuffer.data(df_geom_20_off + 347 * ccomps * dcomps);

            auto g_zz_0_yz_yzz = cbuffer.data(df_geom_20_off + 348 * ccomps * dcomps);

            auto g_zz_0_yz_zzz = cbuffer.data(df_geom_20_off + 349 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yz_xxx, g_zz_0_yz_xxy, g_zz_0_yz_xxz, g_zz_0_yz_xyy, g_zz_0_yz_xyz, g_zz_0_yz_xzz, g_zz_0_yz_yyy, g_zz_0_yz_yyz, g_zz_0_yz_yzz, g_zz_0_yz_zzz, g_zz_0_z_xxx, g_zz_0_z_xxxy, g_zz_0_z_xxy, g_zz_0_z_xxyy, g_zz_0_z_xxyz, g_zz_0_z_xxz, g_zz_0_z_xyy, g_zz_0_z_xyyy, g_zz_0_z_xyyz, g_zz_0_z_xyz, g_zz_0_z_xyzz, g_zz_0_z_xzz, g_zz_0_z_yyy, g_zz_0_z_yyyy, g_zz_0_z_yyyz, g_zz_0_z_yyz, g_zz_0_z_yyzz, g_zz_0_z_yzz, g_zz_0_z_yzzz, g_zz_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yz_xxx[k] = -g_zz_0_z_xxx[k] * ab_y + g_zz_0_z_xxxy[k];

                g_zz_0_yz_xxy[k] = -g_zz_0_z_xxy[k] * ab_y + g_zz_0_z_xxyy[k];

                g_zz_0_yz_xxz[k] = -g_zz_0_z_xxz[k] * ab_y + g_zz_0_z_xxyz[k];

                g_zz_0_yz_xyy[k] = -g_zz_0_z_xyy[k] * ab_y + g_zz_0_z_xyyy[k];

                g_zz_0_yz_xyz[k] = -g_zz_0_z_xyz[k] * ab_y + g_zz_0_z_xyyz[k];

                g_zz_0_yz_xzz[k] = -g_zz_0_z_xzz[k] * ab_y + g_zz_0_z_xyzz[k];

                g_zz_0_yz_yyy[k] = -g_zz_0_z_yyy[k] * ab_y + g_zz_0_z_yyyy[k];

                g_zz_0_yz_yyz[k] = -g_zz_0_z_yyz[k] * ab_y + g_zz_0_z_yyyz[k];

                g_zz_0_yz_yzz[k] = -g_zz_0_z_yzz[k] * ab_y + g_zz_0_z_yyzz[k];

                g_zz_0_yz_zzz[k] = -g_zz_0_z_zzz[k] * ab_y + g_zz_0_z_yzzz[k];
            }

            /// Set up 350-360 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zz_xxx = cbuffer.data(df_geom_20_off + 350 * ccomps * dcomps);

            auto g_zz_0_zz_xxy = cbuffer.data(df_geom_20_off + 351 * ccomps * dcomps);

            auto g_zz_0_zz_xxz = cbuffer.data(df_geom_20_off + 352 * ccomps * dcomps);

            auto g_zz_0_zz_xyy = cbuffer.data(df_geom_20_off + 353 * ccomps * dcomps);

            auto g_zz_0_zz_xyz = cbuffer.data(df_geom_20_off + 354 * ccomps * dcomps);

            auto g_zz_0_zz_xzz = cbuffer.data(df_geom_20_off + 355 * ccomps * dcomps);

            auto g_zz_0_zz_yyy = cbuffer.data(df_geom_20_off + 356 * ccomps * dcomps);

            auto g_zz_0_zz_yyz = cbuffer.data(df_geom_20_off + 357 * ccomps * dcomps);

            auto g_zz_0_zz_yzz = cbuffer.data(df_geom_20_off + 358 * ccomps * dcomps);

            auto g_zz_0_zz_zzz = cbuffer.data(df_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz, g_zz_0_z_xxx, g_zz_0_z_xxxz, g_zz_0_z_xxy, g_zz_0_z_xxyz, g_zz_0_z_xxz, g_zz_0_z_xxzz, g_zz_0_z_xyy, g_zz_0_z_xyyz, g_zz_0_z_xyz, g_zz_0_z_xyzz, g_zz_0_z_xzz, g_zz_0_z_xzzz, g_zz_0_z_yyy, g_zz_0_z_yyyz, g_zz_0_z_yyz, g_zz_0_z_yyzz, g_zz_0_z_yzz, g_zz_0_z_yzzz, g_zz_0_z_zzz, g_zz_0_z_zzzz, g_zz_0_zz_xxx, g_zz_0_zz_xxy, g_zz_0_zz_xxz, g_zz_0_zz_xyy, g_zz_0_zz_xyz, g_zz_0_zz_xzz, g_zz_0_zz_yyy, g_zz_0_zz_yyz, g_zz_0_zz_yzz, g_zz_0_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zz_xxx[k] = -2.0 * g_z_0_z_xxx[k] - g_zz_0_z_xxx[k] * ab_z + g_zz_0_z_xxxz[k];

                g_zz_0_zz_xxy[k] = -2.0 * g_z_0_z_xxy[k] - g_zz_0_z_xxy[k] * ab_z + g_zz_0_z_xxyz[k];

                g_zz_0_zz_xxz[k] = -2.0 * g_z_0_z_xxz[k] - g_zz_0_z_xxz[k] * ab_z + g_zz_0_z_xxzz[k];

                g_zz_0_zz_xyy[k] = -2.0 * g_z_0_z_xyy[k] - g_zz_0_z_xyy[k] * ab_z + g_zz_0_z_xyyz[k];

                g_zz_0_zz_xyz[k] = -2.0 * g_z_0_z_xyz[k] - g_zz_0_z_xyz[k] * ab_z + g_zz_0_z_xyzz[k];

                g_zz_0_zz_xzz[k] = -2.0 * g_z_0_z_xzz[k] - g_zz_0_z_xzz[k] * ab_z + g_zz_0_z_xzzz[k];

                g_zz_0_zz_yyy[k] = -2.0 * g_z_0_z_yyy[k] - g_zz_0_z_yyy[k] * ab_z + g_zz_0_z_yyyz[k];

                g_zz_0_zz_yyz[k] = -2.0 * g_z_0_z_yyz[k] - g_zz_0_z_yyz[k] * ab_z + g_zz_0_z_yyzz[k];

                g_zz_0_zz_yzz[k] = -2.0 * g_z_0_z_yzz[k] - g_zz_0_z_yzz[k] * ab_z + g_zz_0_z_yzzz[k];

                g_zz_0_zz_zzz[k] = -2.0 * g_z_0_z_zzz[k] - g_zz_0_z_zzz[k] * ab_z + g_zz_0_z_zzzz[k];
            }
        }
    }
}

} // erirec namespace

