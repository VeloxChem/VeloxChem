#include "ElectronRepulsionGeom1100ContrRecDFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_dfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_dfxx,
                                            const size_t idx_geom_01_pfxx,
                                            const size_t idx_geom_10_pfxx,
                                            const size_t idx_geom_11_pfxx,
                                            const size_t idx_geom_11_pgxx,
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

            const auto pf_geom_01_off = idx_geom_01_pfxx + i * dcomps + j;

            auto g_0_x_x_xxx = cbuffer.data(pf_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xxy = cbuffer.data(pf_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xxz = cbuffer.data(pf_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_xyy = cbuffer.data(pf_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_xyz = cbuffer.data(pf_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_xzz = cbuffer.data(pf_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_x_yyy = cbuffer.data(pf_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_x_yyz = cbuffer.data(pf_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_x_yzz = cbuffer.data(pf_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_x_zzz = cbuffer.data(pf_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_y_xxx = cbuffer.data(pf_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_y_xxy = cbuffer.data(pf_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_y_xxz = cbuffer.data(pf_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_y_xyy = cbuffer.data(pf_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_y_xyz = cbuffer.data(pf_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_y_xzz = cbuffer.data(pf_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_y_yyy = cbuffer.data(pf_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_y_yyz = cbuffer.data(pf_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_y_yzz = cbuffer.data(pf_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_y_zzz = cbuffer.data(pf_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_z_xxx = cbuffer.data(pf_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_z_xxy = cbuffer.data(pf_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_z_xxz = cbuffer.data(pf_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_z_xyy = cbuffer.data(pf_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_z_xyz = cbuffer.data(pf_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_z_xzz = cbuffer.data(pf_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_z_yyy = cbuffer.data(pf_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_z_yyz = cbuffer.data(pf_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_z_yzz = cbuffer.data(pf_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_z_zzz = cbuffer.data(pf_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_x_xxx = cbuffer.data(pf_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_x_xxy = cbuffer.data(pf_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_x_xxz = cbuffer.data(pf_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_x_xyy = cbuffer.data(pf_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_x_xyz = cbuffer.data(pf_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_x_xzz = cbuffer.data(pf_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_y_x_yyy = cbuffer.data(pf_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_y_x_yyz = cbuffer.data(pf_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_y_x_yzz = cbuffer.data(pf_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_y_x_zzz = cbuffer.data(pf_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_y_y_xxx = cbuffer.data(pf_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_y_y_xxy = cbuffer.data(pf_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_y_y_xxz = cbuffer.data(pf_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_y_y_xyy = cbuffer.data(pf_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_y_y_xyz = cbuffer.data(pf_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_y_xzz = cbuffer.data(pf_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_y_yyy = cbuffer.data(pf_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_y_yyz = cbuffer.data(pf_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_y_yzz = cbuffer.data(pf_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_y_zzz = cbuffer.data(pf_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_z_xxx = cbuffer.data(pf_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_z_xxy = cbuffer.data(pf_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_z_xxz = cbuffer.data(pf_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_z_xyy = cbuffer.data(pf_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_z_xyz = cbuffer.data(pf_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_z_xzz = cbuffer.data(pf_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_y_z_yyy = cbuffer.data(pf_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_y_z_yyz = cbuffer.data(pf_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_y_z_yzz = cbuffer.data(pf_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_y_z_zzz = cbuffer.data(pf_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_z_x_xxx = cbuffer.data(pf_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_z_x_xxy = cbuffer.data(pf_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_z_x_xxz = cbuffer.data(pf_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_z_x_xyy = cbuffer.data(pf_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_z_x_xyz = cbuffer.data(pf_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_z_x_xzz = cbuffer.data(pf_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_z_x_yyy = cbuffer.data(pf_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_z_x_yyz = cbuffer.data(pf_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_z_x_yzz = cbuffer.data(pf_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_z_x_zzz = cbuffer.data(pf_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_z_y_xxx = cbuffer.data(pf_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_z_y_xxy = cbuffer.data(pf_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_z_y_xxz = cbuffer.data(pf_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_z_y_xyy = cbuffer.data(pf_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_z_y_xyz = cbuffer.data(pf_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_z_y_xzz = cbuffer.data(pf_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_z_y_yyy = cbuffer.data(pf_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_z_y_yyz = cbuffer.data(pf_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_z_y_yzz = cbuffer.data(pf_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_z_y_zzz = cbuffer.data(pf_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_z_z_xxx = cbuffer.data(pf_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_z_z_xxy = cbuffer.data(pf_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_z_z_xxz = cbuffer.data(pf_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_z_z_xyy = cbuffer.data(pf_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_z_z_xyz = cbuffer.data(pf_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_z_z_xzz = cbuffer.data(pf_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_z_z_yyy = cbuffer.data(pf_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_z_z_yyz = cbuffer.data(pf_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_z_z_yzz = cbuffer.data(pf_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_z_z_zzz = cbuffer.data(pf_geom_01_off + 89 * ccomps * dcomps);

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

            const auto pf_geom_11_off = idx_geom_11_pfxx + i * dcomps + j;

            auto g_x_x_x_xxx = cbuffer.data(pf_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xxy = cbuffer.data(pf_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xxz = cbuffer.data(pf_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_xyy = cbuffer.data(pf_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_xyz = cbuffer.data(pf_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_xzz = cbuffer.data(pf_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_x_yyy = cbuffer.data(pf_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_x_yyz = cbuffer.data(pf_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_x_yzz = cbuffer.data(pf_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_x_zzz = cbuffer.data(pf_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_y_xxx = cbuffer.data(pf_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_y_xxy = cbuffer.data(pf_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_y_xxz = cbuffer.data(pf_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_y_xyy = cbuffer.data(pf_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_y_xyz = cbuffer.data(pf_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_y_xzz = cbuffer.data(pf_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_y_yyy = cbuffer.data(pf_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_y_yyz = cbuffer.data(pf_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_y_yzz = cbuffer.data(pf_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_y_zzz = cbuffer.data(pf_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_z_xxx = cbuffer.data(pf_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_z_xxy = cbuffer.data(pf_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_z_xxz = cbuffer.data(pf_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_z_xyy = cbuffer.data(pf_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_z_xyz = cbuffer.data(pf_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_z_xzz = cbuffer.data(pf_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_z_yyy = cbuffer.data(pf_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_z_yyz = cbuffer.data(pf_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_z_yzz = cbuffer.data(pf_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_z_zzz = cbuffer.data(pf_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_y_x_xxx = cbuffer.data(pf_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_x_xxy = cbuffer.data(pf_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_x_xxz = cbuffer.data(pf_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_y_x_xyy = cbuffer.data(pf_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_x_xyz = cbuffer.data(pf_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_x_xzz = cbuffer.data(pf_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_y_x_yyy = cbuffer.data(pf_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_x_yyz = cbuffer.data(pf_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_x_yzz = cbuffer.data(pf_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_y_x_zzz = cbuffer.data(pf_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_y_y_xxx = cbuffer.data(pf_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_y_xxy = cbuffer.data(pf_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_y_y_xxz = cbuffer.data(pf_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_y_y_xyy = cbuffer.data(pf_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_y_y_xyz = cbuffer.data(pf_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_y_xzz = cbuffer.data(pf_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_y_yyy = cbuffer.data(pf_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_y_yyz = cbuffer.data(pf_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_y_yzz = cbuffer.data(pf_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_y_zzz = cbuffer.data(pf_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_z_xxx = cbuffer.data(pf_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_z_xxy = cbuffer.data(pf_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_z_xxz = cbuffer.data(pf_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_z_xyy = cbuffer.data(pf_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_z_xyz = cbuffer.data(pf_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_z_xzz = cbuffer.data(pf_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_z_yyy = cbuffer.data(pf_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_y_z_yyz = cbuffer.data(pf_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_z_yzz = cbuffer.data(pf_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_z_zzz = cbuffer.data(pf_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_z_x_xxx = cbuffer.data(pf_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_z_x_xxy = cbuffer.data(pf_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_z_x_xxz = cbuffer.data(pf_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_z_x_xyy = cbuffer.data(pf_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_z_x_xyz = cbuffer.data(pf_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_z_x_xzz = cbuffer.data(pf_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_z_x_yyy = cbuffer.data(pf_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_z_x_yyz = cbuffer.data(pf_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_z_x_yzz = cbuffer.data(pf_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_z_x_zzz = cbuffer.data(pf_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_z_y_xxx = cbuffer.data(pf_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_z_y_xxy = cbuffer.data(pf_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_z_y_xxz = cbuffer.data(pf_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_z_y_xyy = cbuffer.data(pf_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_y_xyz = cbuffer.data(pf_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_z_y_xzz = cbuffer.data(pf_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_z_y_yyy = cbuffer.data(pf_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_y_yyz = cbuffer.data(pf_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_z_y_yzz = cbuffer.data(pf_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_z_y_zzz = cbuffer.data(pf_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_z_z_xxx = cbuffer.data(pf_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_z_z_xxy = cbuffer.data(pf_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_z_z_xxz = cbuffer.data(pf_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_z_z_xyy = cbuffer.data(pf_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_z_z_xyz = cbuffer.data(pf_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_z_z_xzz = cbuffer.data(pf_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_z_z_yyy = cbuffer.data(pf_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_z_z_yyz = cbuffer.data(pf_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_z_z_yzz = cbuffer.data(pf_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_z_z_zzz = cbuffer.data(pf_geom_11_off + 89 * ccomps * dcomps);

            auto g_y_x_x_xxx = cbuffer.data(pf_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_x_x_xxy = cbuffer.data(pf_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_x_x_xxz = cbuffer.data(pf_geom_11_off + 92 * ccomps * dcomps);

            auto g_y_x_x_xyy = cbuffer.data(pf_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_x_x_xyz = cbuffer.data(pf_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_x_x_xzz = cbuffer.data(pf_geom_11_off + 95 * ccomps * dcomps);

            auto g_y_x_x_yyy = cbuffer.data(pf_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_x_x_yyz = cbuffer.data(pf_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_x_x_yzz = cbuffer.data(pf_geom_11_off + 98 * ccomps * dcomps);

            auto g_y_x_x_zzz = cbuffer.data(pf_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_x_y_xxx = cbuffer.data(pf_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_x_y_xxy = cbuffer.data(pf_geom_11_off + 101 * ccomps * dcomps);

            auto g_y_x_y_xxz = cbuffer.data(pf_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_x_y_xyy = cbuffer.data(pf_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_x_y_xyz = cbuffer.data(pf_geom_11_off + 104 * ccomps * dcomps);

            auto g_y_x_y_xzz = cbuffer.data(pf_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_x_y_yyy = cbuffer.data(pf_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_x_y_yyz = cbuffer.data(pf_geom_11_off + 107 * ccomps * dcomps);

            auto g_y_x_y_yzz = cbuffer.data(pf_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_x_y_zzz = cbuffer.data(pf_geom_11_off + 109 * ccomps * dcomps);

            auto g_y_x_z_xxx = cbuffer.data(pf_geom_11_off + 110 * ccomps * dcomps);

            auto g_y_x_z_xxy = cbuffer.data(pf_geom_11_off + 111 * ccomps * dcomps);

            auto g_y_x_z_xxz = cbuffer.data(pf_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_x_z_xyy = cbuffer.data(pf_geom_11_off + 113 * ccomps * dcomps);

            auto g_y_x_z_xyz = cbuffer.data(pf_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_x_z_xzz = cbuffer.data(pf_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_x_z_yyy = cbuffer.data(pf_geom_11_off + 116 * ccomps * dcomps);

            auto g_y_x_z_yyz = cbuffer.data(pf_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_x_z_yzz = cbuffer.data(pf_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_x_z_zzz = cbuffer.data(pf_geom_11_off + 119 * ccomps * dcomps);

            auto g_y_y_x_xxx = cbuffer.data(pf_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_y_x_xxy = cbuffer.data(pf_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_y_x_xxz = cbuffer.data(pf_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_y_x_xyy = cbuffer.data(pf_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_y_x_xyz = cbuffer.data(pf_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_y_x_xzz = cbuffer.data(pf_geom_11_off + 125 * ccomps * dcomps);

            auto g_y_y_x_yyy = cbuffer.data(pf_geom_11_off + 126 * ccomps * dcomps);

            auto g_y_y_x_yyz = cbuffer.data(pf_geom_11_off + 127 * ccomps * dcomps);

            auto g_y_y_x_yzz = cbuffer.data(pf_geom_11_off + 128 * ccomps * dcomps);

            auto g_y_y_x_zzz = cbuffer.data(pf_geom_11_off + 129 * ccomps * dcomps);

            auto g_y_y_y_xxx = cbuffer.data(pf_geom_11_off + 130 * ccomps * dcomps);

            auto g_y_y_y_xxy = cbuffer.data(pf_geom_11_off + 131 * ccomps * dcomps);

            auto g_y_y_y_xxz = cbuffer.data(pf_geom_11_off + 132 * ccomps * dcomps);

            auto g_y_y_y_xyy = cbuffer.data(pf_geom_11_off + 133 * ccomps * dcomps);

            auto g_y_y_y_xyz = cbuffer.data(pf_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_y_y_xzz = cbuffer.data(pf_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_y_y_yyy = cbuffer.data(pf_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_y_y_yyz = cbuffer.data(pf_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_y_y_yzz = cbuffer.data(pf_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_y_y_zzz = cbuffer.data(pf_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_y_z_xxx = cbuffer.data(pf_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_y_z_xxy = cbuffer.data(pf_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_y_z_xxz = cbuffer.data(pf_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_y_z_xyy = cbuffer.data(pf_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_y_z_xyz = cbuffer.data(pf_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_y_z_xzz = cbuffer.data(pf_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_y_z_yyy = cbuffer.data(pf_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_y_z_yyz = cbuffer.data(pf_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_y_z_yzz = cbuffer.data(pf_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_y_z_zzz = cbuffer.data(pf_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_z_x_xxx = cbuffer.data(pf_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_z_x_xxy = cbuffer.data(pf_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_z_x_xxz = cbuffer.data(pf_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_z_x_xyy = cbuffer.data(pf_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_z_x_xyz = cbuffer.data(pf_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_z_x_xzz = cbuffer.data(pf_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_z_x_yyy = cbuffer.data(pf_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_z_x_yyz = cbuffer.data(pf_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_z_x_yzz = cbuffer.data(pf_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_z_x_zzz = cbuffer.data(pf_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_z_y_xxx = cbuffer.data(pf_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_z_y_xxy = cbuffer.data(pf_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_z_y_xxz = cbuffer.data(pf_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_z_y_xyy = cbuffer.data(pf_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_z_y_xyz = cbuffer.data(pf_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_z_y_xzz = cbuffer.data(pf_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_z_y_yyy = cbuffer.data(pf_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_z_y_yyz = cbuffer.data(pf_geom_11_off + 167 * ccomps * dcomps);

            auto g_y_z_y_yzz = cbuffer.data(pf_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_z_y_zzz = cbuffer.data(pf_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_z_z_xxx = cbuffer.data(pf_geom_11_off + 170 * ccomps * dcomps);

            auto g_y_z_z_xxy = cbuffer.data(pf_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_z_z_xxz = cbuffer.data(pf_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_z_z_xyy = cbuffer.data(pf_geom_11_off + 173 * ccomps * dcomps);

            auto g_y_z_z_xyz = cbuffer.data(pf_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_z_z_xzz = cbuffer.data(pf_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_z_z_yyy = cbuffer.data(pf_geom_11_off + 176 * ccomps * dcomps);

            auto g_y_z_z_yyz = cbuffer.data(pf_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_z_z_yzz = cbuffer.data(pf_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_z_z_zzz = cbuffer.data(pf_geom_11_off + 179 * ccomps * dcomps);

            auto g_z_x_x_xxx = cbuffer.data(pf_geom_11_off + 180 * ccomps * dcomps);

            auto g_z_x_x_xxy = cbuffer.data(pf_geom_11_off + 181 * ccomps * dcomps);

            auto g_z_x_x_xxz = cbuffer.data(pf_geom_11_off + 182 * ccomps * dcomps);

            auto g_z_x_x_xyy = cbuffer.data(pf_geom_11_off + 183 * ccomps * dcomps);

            auto g_z_x_x_xyz = cbuffer.data(pf_geom_11_off + 184 * ccomps * dcomps);

            auto g_z_x_x_xzz = cbuffer.data(pf_geom_11_off + 185 * ccomps * dcomps);

            auto g_z_x_x_yyy = cbuffer.data(pf_geom_11_off + 186 * ccomps * dcomps);

            auto g_z_x_x_yyz = cbuffer.data(pf_geom_11_off + 187 * ccomps * dcomps);

            auto g_z_x_x_yzz = cbuffer.data(pf_geom_11_off + 188 * ccomps * dcomps);

            auto g_z_x_x_zzz = cbuffer.data(pf_geom_11_off + 189 * ccomps * dcomps);

            auto g_z_x_y_xxx = cbuffer.data(pf_geom_11_off + 190 * ccomps * dcomps);

            auto g_z_x_y_xxy = cbuffer.data(pf_geom_11_off + 191 * ccomps * dcomps);

            auto g_z_x_y_xxz = cbuffer.data(pf_geom_11_off + 192 * ccomps * dcomps);

            auto g_z_x_y_xyy = cbuffer.data(pf_geom_11_off + 193 * ccomps * dcomps);

            auto g_z_x_y_xyz = cbuffer.data(pf_geom_11_off + 194 * ccomps * dcomps);

            auto g_z_x_y_xzz = cbuffer.data(pf_geom_11_off + 195 * ccomps * dcomps);

            auto g_z_x_y_yyy = cbuffer.data(pf_geom_11_off + 196 * ccomps * dcomps);

            auto g_z_x_y_yyz = cbuffer.data(pf_geom_11_off + 197 * ccomps * dcomps);

            auto g_z_x_y_yzz = cbuffer.data(pf_geom_11_off + 198 * ccomps * dcomps);

            auto g_z_x_y_zzz = cbuffer.data(pf_geom_11_off + 199 * ccomps * dcomps);

            auto g_z_x_z_xxx = cbuffer.data(pf_geom_11_off + 200 * ccomps * dcomps);

            auto g_z_x_z_xxy = cbuffer.data(pf_geom_11_off + 201 * ccomps * dcomps);

            auto g_z_x_z_xxz = cbuffer.data(pf_geom_11_off + 202 * ccomps * dcomps);

            auto g_z_x_z_xyy = cbuffer.data(pf_geom_11_off + 203 * ccomps * dcomps);

            auto g_z_x_z_xyz = cbuffer.data(pf_geom_11_off + 204 * ccomps * dcomps);

            auto g_z_x_z_xzz = cbuffer.data(pf_geom_11_off + 205 * ccomps * dcomps);

            auto g_z_x_z_yyy = cbuffer.data(pf_geom_11_off + 206 * ccomps * dcomps);

            auto g_z_x_z_yyz = cbuffer.data(pf_geom_11_off + 207 * ccomps * dcomps);

            auto g_z_x_z_yzz = cbuffer.data(pf_geom_11_off + 208 * ccomps * dcomps);

            auto g_z_x_z_zzz = cbuffer.data(pf_geom_11_off + 209 * ccomps * dcomps);

            auto g_z_y_x_xxx = cbuffer.data(pf_geom_11_off + 210 * ccomps * dcomps);

            auto g_z_y_x_xxy = cbuffer.data(pf_geom_11_off + 211 * ccomps * dcomps);

            auto g_z_y_x_xxz = cbuffer.data(pf_geom_11_off + 212 * ccomps * dcomps);

            auto g_z_y_x_xyy = cbuffer.data(pf_geom_11_off + 213 * ccomps * dcomps);

            auto g_z_y_x_xyz = cbuffer.data(pf_geom_11_off + 214 * ccomps * dcomps);

            auto g_z_y_x_xzz = cbuffer.data(pf_geom_11_off + 215 * ccomps * dcomps);

            auto g_z_y_x_yyy = cbuffer.data(pf_geom_11_off + 216 * ccomps * dcomps);

            auto g_z_y_x_yyz = cbuffer.data(pf_geom_11_off + 217 * ccomps * dcomps);

            auto g_z_y_x_yzz = cbuffer.data(pf_geom_11_off + 218 * ccomps * dcomps);

            auto g_z_y_x_zzz = cbuffer.data(pf_geom_11_off + 219 * ccomps * dcomps);

            auto g_z_y_y_xxx = cbuffer.data(pf_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_y_y_xxy = cbuffer.data(pf_geom_11_off + 221 * ccomps * dcomps);

            auto g_z_y_y_xxz = cbuffer.data(pf_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_y_y_xyy = cbuffer.data(pf_geom_11_off + 223 * ccomps * dcomps);

            auto g_z_y_y_xyz = cbuffer.data(pf_geom_11_off + 224 * ccomps * dcomps);

            auto g_z_y_y_xzz = cbuffer.data(pf_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_y_y_yyy = cbuffer.data(pf_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_y_y_yyz = cbuffer.data(pf_geom_11_off + 227 * ccomps * dcomps);

            auto g_z_y_y_yzz = cbuffer.data(pf_geom_11_off + 228 * ccomps * dcomps);

            auto g_z_y_y_zzz = cbuffer.data(pf_geom_11_off + 229 * ccomps * dcomps);

            auto g_z_y_z_xxx = cbuffer.data(pf_geom_11_off + 230 * ccomps * dcomps);

            auto g_z_y_z_xxy = cbuffer.data(pf_geom_11_off + 231 * ccomps * dcomps);

            auto g_z_y_z_xxz = cbuffer.data(pf_geom_11_off + 232 * ccomps * dcomps);

            auto g_z_y_z_xyy = cbuffer.data(pf_geom_11_off + 233 * ccomps * dcomps);

            auto g_z_y_z_xyz = cbuffer.data(pf_geom_11_off + 234 * ccomps * dcomps);

            auto g_z_y_z_xzz = cbuffer.data(pf_geom_11_off + 235 * ccomps * dcomps);

            auto g_z_y_z_yyy = cbuffer.data(pf_geom_11_off + 236 * ccomps * dcomps);

            auto g_z_y_z_yyz = cbuffer.data(pf_geom_11_off + 237 * ccomps * dcomps);

            auto g_z_y_z_yzz = cbuffer.data(pf_geom_11_off + 238 * ccomps * dcomps);

            auto g_z_y_z_zzz = cbuffer.data(pf_geom_11_off + 239 * ccomps * dcomps);

            auto g_z_z_x_xxx = cbuffer.data(pf_geom_11_off + 240 * ccomps * dcomps);

            auto g_z_z_x_xxy = cbuffer.data(pf_geom_11_off + 241 * ccomps * dcomps);

            auto g_z_z_x_xxz = cbuffer.data(pf_geom_11_off + 242 * ccomps * dcomps);

            auto g_z_z_x_xyy = cbuffer.data(pf_geom_11_off + 243 * ccomps * dcomps);

            auto g_z_z_x_xyz = cbuffer.data(pf_geom_11_off + 244 * ccomps * dcomps);

            auto g_z_z_x_xzz = cbuffer.data(pf_geom_11_off + 245 * ccomps * dcomps);

            auto g_z_z_x_yyy = cbuffer.data(pf_geom_11_off + 246 * ccomps * dcomps);

            auto g_z_z_x_yyz = cbuffer.data(pf_geom_11_off + 247 * ccomps * dcomps);

            auto g_z_z_x_yzz = cbuffer.data(pf_geom_11_off + 248 * ccomps * dcomps);

            auto g_z_z_x_zzz = cbuffer.data(pf_geom_11_off + 249 * ccomps * dcomps);

            auto g_z_z_y_xxx = cbuffer.data(pf_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_z_y_xxy = cbuffer.data(pf_geom_11_off + 251 * ccomps * dcomps);

            auto g_z_z_y_xxz = cbuffer.data(pf_geom_11_off + 252 * ccomps * dcomps);

            auto g_z_z_y_xyy = cbuffer.data(pf_geom_11_off + 253 * ccomps * dcomps);

            auto g_z_z_y_xyz = cbuffer.data(pf_geom_11_off + 254 * ccomps * dcomps);

            auto g_z_z_y_xzz = cbuffer.data(pf_geom_11_off + 255 * ccomps * dcomps);

            auto g_z_z_y_yyy = cbuffer.data(pf_geom_11_off + 256 * ccomps * dcomps);

            auto g_z_z_y_yyz = cbuffer.data(pf_geom_11_off + 257 * ccomps * dcomps);

            auto g_z_z_y_yzz = cbuffer.data(pf_geom_11_off + 258 * ccomps * dcomps);

            auto g_z_z_y_zzz = cbuffer.data(pf_geom_11_off + 259 * ccomps * dcomps);

            auto g_z_z_z_xxx = cbuffer.data(pf_geom_11_off + 260 * ccomps * dcomps);

            auto g_z_z_z_xxy = cbuffer.data(pf_geom_11_off + 261 * ccomps * dcomps);

            auto g_z_z_z_xxz = cbuffer.data(pf_geom_11_off + 262 * ccomps * dcomps);

            auto g_z_z_z_xyy = cbuffer.data(pf_geom_11_off + 263 * ccomps * dcomps);

            auto g_z_z_z_xyz = cbuffer.data(pf_geom_11_off + 264 * ccomps * dcomps);

            auto g_z_z_z_xzz = cbuffer.data(pf_geom_11_off + 265 * ccomps * dcomps);

            auto g_z_z_z_yyy = cbuffer.data(pf_geom_11_off + 266 * ccomps * dcomps);

            auto g_z_z_z_yyz = cbuffer.data(pf_geom_11_off + 267 * ccomps * dcomps);

            auto g_z_z_z_yzz = cbuffer.data(pf_geom_11_off + 268 * ccomps * dcomps);

            auto g_z_z_z_zzz = cbuffer.data(pf_geom_11_off + 269 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PGSS

            const auto pg_geom_11_off = idx_geom_11_pgxx + i * dcomps + j;

            auto g_x_x_x_xxxx = cbuffer.data(pg_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xxxy = cbuffer.data(pg_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xxxz = cbuffer.data(pg_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_xxyy = cbuffer.data(pg_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_xxyz = cbuffer.data(pg_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_xxzz = cbuffer.data(pg_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_x_xyyy = cbuffer.data(pg_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_x_xyyz = cbuffer.data(pg_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_x_xyzz = cbuffer.data(pg_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_x_xzzz = cbuffer.data(pg_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_x_yyyy = cbuffer.data(pg_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_x_yyyz = cbuffer.data(pg_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_x_yyzz = cbuffer.data(pg_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_x_yzzz = cbuffer.data(pg_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_x_zzzz = cbuffer.data(pg_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_y_xxxx = cbuffer.data(pg_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_y_xxxy = cbuffer.data(pg_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_y_xxxz = cbuffer.data(pg_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_y_xxyy = cbuffer.data(pg_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_y_xxyz = cbuffer.data(pg_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_y_xxzz = cbuffer.data(pg_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_y_xyyy = cbuffer.data(pg_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_y_xyyz = cbuffer.data(pg_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_y_xyzz = cbuffer.data(pg_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_y_xzzz = cbuffer.data(pg_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_y_yyyy = cbuffer.data(pg_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_y_yyyz = cbuffer.data(pg_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_y_yyzz = cbuffer.data(pg_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_y_yzzz = cbuffer.data(pg_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_y_zzzz = cbuffer.data(pg_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_z_xxxx = cbuffer.data(pg_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_z_xxxy = cbuffer.data(pg_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_z_xxxz = cbuffer.data(pg_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_z_xxyy = cbuffer.data(pg_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_z_xxyz = cbuffer.data(pg_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_z_xxzz = cbuffer.data(pg_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_z_xyyy = cbuffer.data(pg_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_z_xyyz = cbuffer.data(pg_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_z_xyzz = cbuffer.data(pg_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_z_xzzz = cbuffer.data(pg_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_z_yyyy = cbuffer.data(pg_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_z_yyyz = cbuffer.data(pg_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_z_yyzz = cbuffer.data(pg_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_z_yzzz = cbuffer.data(pg_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_z_zzzz = cbuffer.data(pg_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_x_xxxx = cbuffer.data(pg_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_x_xxxy = cbuffer.data(pg_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_x_xxxz = cbuffer.data(pg_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_x_xxyy = cbuffer.data(pg_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_x_xxyz = cbuffer.data(pg_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_x_xxzz = cbuffer.data(pg_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_x_xyyy = cbuffer.data(pg_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_x_xyyz = cbuffer.data(pg_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_x_xyzz = cbuffer.data(pg_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_x_xzzz = cbuffer.data(pg_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_x_yyyy = cbuffer.data(pg_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_x_yyyz = cbuffer.data(pg_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_y_x_yyzz = cbuffer.data(pg_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_x_yzzz = cbuffer.data(pg_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_x_zzzz = cbuffer.data(pg_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_y_y_xxxx = cbuffer.data(pg_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_y_xxxy = cbuffer.data(pg_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_y_xxxz = cbuffer.data(pg_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_y_xxyy = cbuffer.data(pg_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_y_xxyz = cbuffer.data(pg_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_y_xxzz = cbuffer.data(pg_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_y_xyyy = cbuffer.data(pg_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_y_xyyz = cbuffer.data(pg_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_y_xyzz = cbuffer.data(pg_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_y_xzzz = cbuffer.data(pg_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_y_yyyy = cbuffer.data(pg_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_y_yyyz = cbuffer.data(pg_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_y_y_yyzz = cbuffer.data(pg_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_y_yzzz = cbuffer.data(pg_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_y_zzzz = cbuffer.data(pg_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_y_z_xxxx = cbuffer.data(pg_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_z_xxxy = cbuffer.data(pg_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_z_xxxz = cbuffer.data(pg_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_y_z_xxyy = cbuffer.data(pg_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_z_xxyz = cbuffer.data(pg_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_y_z_xxzz = cbuffer.data(pg_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_y_z_xyyy = cbuffer.data(pg_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_z_xyyz = cbuffer.data(pg_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_z_xyzz = cbuffer.data(pg_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_y_z_xzzz = cbuffer.data(pg_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_z_yyyy = cbuffer.data(pg_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_z_yyyz = cbuffer.data(pg_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_z_yyzz = cbuffer.data(pg_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_z_yzzz = cbuffer.data(pg_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_z_zzzz = cbuffer.data(pg_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_z_x_xxxx = cbuffer.data(pg_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_z_x_xxxy = cbuffer.data(pg_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_z_x_xxxz = cbuffer.data(pg_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_z_x_xxyy = cbuffer.data(pg_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_z_x_xxyz = cbuffer.data(pg_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_z_x_xxzz = cbuffer.data(pg_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_z_x_xyyy = cbuffer.data(pg_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_z_x_xyyz = cbuffer.data(pg_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_z_x_xyzz = cbuffer.data(pg_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_z_x_xzzz = cbuffer.data(pg_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_z_x_yyyy = cbuffer.data(pg_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_z_x_yyyz = cbuffer.data(pg_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_z_x_yyzz = cbuffer.data(pg_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_z_x_yzzz = cbuffer.data(pg_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_z_x_zzzz = cbuffer.data(pg_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_z_y_xxxx = cbuffer.data(pg_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_z_y_xxxy = cbuffer.data(pg_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_z_y_xxxz = cbuffer.data(pg_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_z_y_xxyy = cbuffer.data(pg_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_z_y_xxyz = cbuffer.data(pg_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_z_y_xxzz = cbuffer.data(pg_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_z_y_xyyy = cbuffer.data(pg_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_z_y_xyyz = cbuffer.data(pg_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_z_y_xyzz = cbuffer.data(pg_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_z_y_xzzz = cbuffer.data(pg_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_z_y_yyyy = cbuffer.data(pg_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_z_y_yyyz = cbuffer.data(pg_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_z_y_yyzz = cbuffer.data(pg_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_z_y_yzzz = cbuffer.data(pg_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_z_y_zzzz = cbuffer.data(pg_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_z_z_xxxx = cbuffer.data(pg_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_z_z_xxxy = cbuffer.data(pg_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_z_z_xxxz = cbuffer.data(pg_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_z_z_xxyy = cbuffer.data(pg_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_z_z_xxyz = cbuffer.data(pg_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_z_z_xxzz = cbuffer.data(pg_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_z_z_xyyy = cbuffer.data(pg_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_z_xyyz = cbuffer.data(pg_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_z_xyzz = cbuffer.data(pg_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_z_z_xzzz = cbuffer.data(pg_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_z_z_yyyy = cbuffer.data(pg_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_z_yyyz = cbuffer.data(pg_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_z_z_yyzz = cbuffer.data(pg_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_z_yzzz = cbuffer.data(pg_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_z_zzzz = cbuffer.data(pg_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_x_x_xxxx = cbuffer.data(pg_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_x_x_xxxy = cbuffer.data(pg_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_x_x_xxxz = cbuffer.data(pg_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_x_x_xxyy = cbuffer.data(pg_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_x_x_xxyz = cbuffer.data(pg_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_x_x_xxzz = cbuffer.data(pg_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_x_x_xyyy = cbuffer.data(pg_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_x_x_xyyz = cbuffer.data(pg_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_x_x_xyzz = cbuffer.data(pg_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_x_x_xzzz = cbuffer.data(pg_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_x_x_yyyy = cbuffer.data(pg_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_x_x_yyyz = cbuffer.data(pg_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_x_x_yyzz = cbuffer.data(pg_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_x_x_yzzz = cbuffer.data(pg_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_x_x_zzzz = cbuffer.data(pg_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_x_y_xxxx = cbuffer.data(pg_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_x_y_xxxy = cbuffer.data(pg_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_x_y_xxxz = cbuffer.data(pg_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_x_y_xxyy = cbuffer.data(pg_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_x_y_xxyz = cbuffer.data(pg_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_x_y_xxzz = cbuffer.data(pg_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_x_y_xyyy = cbuffer.data(pg_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_x_y_xyyz = cbuffer.data(pg_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_x_y_xyzz = cbuffer.data(pg_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_x_y_xzzz = cbuffer.data(pg_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_x_y_yyyy = cbuffer.data(pg_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_x_y_yyyz = cbuffer.data(pg_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_x_y_yyzz = cbuffer.data(pg_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_x_y_yzzz = cbuffer.data(pg_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_x_y_zzzz = cbuffer.data(pg_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_x_z_xxxx = cbuffer.data(pg_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_x_z_xxxy = cbuffer.data(pg_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_x_z_xxxz = cbuffer.data(pg_geom_11_off + 167 * ccomps * dcomps);

            auto g_y_x_z_xxyy = cbuffer.data(pg_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_x_z_xxyz = cbuffer.data(pg_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_x_z_xxzz = cbuffer.data(pg_geom_11_off + 170 * ccomps * dcomps);

            auto g_y_x_z_xyyy = cbuffer.data(pg_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_x_z_xyyz = cbuffer.data(pg_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_x_z_xyzz = cbuffer.data(pg_geom_11_off + 173 * ccomps * dcomps);

            auto g_y_x_z_xzzz = cbuffer.data(pg_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_x_z_yyyy = cbuffer.data(pg_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_x_z_yyyz = cbuffer.data(pg_geom_11_off + 176 * ccomps * dcomps);

            auto g_y_x_z_yyzz = cbuffer.data(pg_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_x_z_yzzz = cbuffer.data(pg_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_x_z_zzzz = cbuffer.data(pg_geom_11_off + 179 * ccomps * dcomps);

            auto g_y_y_x_xxxx = cbuffer.data(pg_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_y_x_xxxy = cbuffer.data(pg_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_y_x_xxxz = cbuffer.data(pg_geom_11_off + 182 * ccomps * dcomps);

            auto g_y_y_x_xxyy = cbuffer.data(pg_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_y_x_xxyz = cbuffer.data(pg_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_y_x_xxzz = cbuffer.data(pg_geom_11_off + 185 * ccomps * dcomps);

            auto g_y_y_x_xyyy = cbuffer.data(pg_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_y_x_xyyz = cbuffer.data(pg_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_y_x_xyzz = cbuffer.data(pg_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_y_x_xzzz = cbuffer.data(pg_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_y_x_yyyy = cbuffer.data(pg_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_y_x_yyyz = cbuffer.data(pg_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_y_x_yyzz = cbuffer.data(pg_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_y_x_yzzz = cbuffer.data(pg_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_y_x_zzzz = cbuffer.data(pg_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_y_y_xxxx = cbuffer.data(pg_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_y_y_xxxy = cbuffer.data(pg_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_y_y_xxxz = cbuffer.data(pg_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_y_y_xxyy = cbuffer.data(pg_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_y_y_xxyz = cbuffer.data(pg_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_y_y_xxzz = cbuffer.data(pg_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_y_y_xyyy = cbuffer.data(pg_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_y_y_xyyz = cbuffer.data(pg_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_y_y_xyzz = cbuffer.data(pg_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_y_y_xzzz = cbuffer.data(pg_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_y_y_yyyy = cbuffer.data(pg_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_y_y_yyyz = cbuffer.data(pg_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_y_y_yyzz = cbuffer.data(pg_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_y_y_yzzz = cbuffer.data(pg_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_y_y_zzzz = cbuffer.data(pg_geom_11_off + 209 * ccomps * dcomps);

            auto g_y_y_z_xxxx = cbuffer.data(pg_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_y_z_xxxy = cbuffer.data(pg_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_y_z_xxxz = cbuffer.data(pg_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_y_z_xxyy = cbuffer.data(pg_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_y_z_xxyz = cbuffer.data(pg_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_y_z_xxzz = cbuffer.data(pg_geom_11_off + 215 * ccomps * dcomps);

            auto g_y_y_z_xyyy = cbuffer.data(pg_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_y_z_xyyz = cbuffer.data(pg_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_y_z_xyzz = cbuffer.data(pg_geom_11_off + 218 * ccomps * dcomps);

            auto g_y_y_z_xzzz = cbuffer.data(pg_geom_11_off + 219 * ccomps * dcomps);

            auto g_y_y_z_yyyy = cbuffer.data(pg_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_y_z_yyyz = cbuffer.data(pg_geom_11_off + 221 * ccomps * dcomps);

            auto g_y_y_z_yyzz = cbuffer.data(pg_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_y_z_yzzz = cbuffer.data(pg_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_y_z_zzzz = cbuffer.data(pg_geom_11_off + 224 * ccomps * dcomps);

            auto g_y_z_x_xxxx = cbuffer.data(pg_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_z_x_xxxy = cbuffer.data(pg_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_z_x_xxxz = cbuffer.data(pg_geom_11_off + 227 * ccomps * dcomps);

            auto g_y_z_x_xxyy = cbuffer.data(pg_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_z_x_xxyz = cbuffer.data(pg_geom_11_off + 229 * ccomps * dcomps);

            auto g_y_z_x_xxzz = cbuffer.data(pg_geom_11_off + 230 * ccomps * dcomps);

            auto g_y_z_x_xyyy = cbuffer.data(pg_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_z_x_xyyz = cbuffer.data(pg_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_z_x_xyzz = cbuffer.data(pg_geom_11_off + 233 * ccomps * dcomps);

            auto g_y_z_x_xzzz = cbuffer.data(pg_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_z_x_yyyy = cbuffer.data(pg_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_z_x_yyyz = cbuffer.data(pg_geom_11_off + 236 * ccomps * dcomps);

            auto g_y_z_x_yyzz = cbuffer.data(pg_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_z_x_yzzz = cbuffer.data(pg_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_z_x_zzzz = cbuffer.data(pg_geom_11_off + 239 * ccomps * dcomps);

            auto g_y_z_y_xxxx = cbuffer.data(pg_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_z_y_xxxy = cbuffer.data(pg_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_z_y_xxxz = cbuffer.data(pg_geom_11_off + 242 * ccomps * dcomps);

            auto g_y_z_y_xxyy = cbuffer.data(pg_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_z_y_xxyz = cbuffer.data(pg_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_z_y_xxzz = cbuffer.data(pg_geom_11_off + 245 * ccomps * dcomps);

            auto g_y_z_y_xyyy = cbuffer.data(pg_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_z_y_xyyz = cbuffer.data(pg_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_z_y_xyzz = cbuffer.data(pg_geom_11_off + 248 * ccomps * dcomps);

            auto g_y_z_y_xzzz = cbuffer.data(pg_geom_11_off + 249 * ccomps * dcomps);

            auto g_y_z_y_yyyy = cbuffer.data(pg_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_z_y_yyyz = cbuffer.data(pg_geom_11_off + 251 * ccomps * dcomps);

            auto g_y_z_y_yyzz = cbuffer.data(pg_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_z_y_yzzz = cbuffer.data(pg_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_z_y_zzzz = cbuffer.data(pg_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_z_z_xxxx = cbuffer.data(pg_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_z_z_xxxy = cbuffer.data(pg_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_z_z_xxxz = cbuffer.data(pg_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_z_z_xxyy = cbuffer.data(pg_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_z_z_xxyz = cbuffer.data(pg_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_z_z_xxzz = cbuffer.data(pg_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_z_z_xyyy = cbuffer.data(pg_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_z_z_xyyz = cbuffer.data(pg_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_z_z_xyzz = cbuffer.data(pg_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_z_z_xzzz = cbuffer.data(pg_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_z_z_yyyy = cbuffer.data(pg_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_z_z_yyyz = cbuffer.data(pg_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_z_z_yyzz = cbuffer.data(pg_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_z_z_yzzz = cbuffer.data(pg_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_z_z_zzzz = cbuffer.data(pg_geom_11_off + 269 * ccomps * dcomps);

            auto g_z_x_x_xxxx = cbuffer.data(pg_geom_11_off + 270 * ccomps * dcomps);

            auto g_z_x_x_xxxy = cbuffer.data(pg_geom_11_off + 271 * ccomps * dcomps);

            auto g_z_x_x_xxxz = cbuffer.data(pg_geom_11_off + 272 * ccomps * dcomps);

            auto g_z_x_x_xxyy = cbuffer.data(pg_geom_11_off + 273 * ccomps * dcomps);

            auto g_z_x_x_xxyz = cbuffer.data(pg_geom_11_off + 274 * ccomps * dcomps);

            auto g_z_x_x_xxzz = cbuffer.data(pg_geom_11_off + 275 * ccomps * dcomps);

            auto g_z_x_x_xyyy = cbuffer.data(pg_geom_11_off + 276 * ccomps * dcomps);

            auto g_z_x_x_xyyz = cbuffer.data(pg_geom_11_off + 277 * ccomps * dcomps);

            auto g_z_x_x_xyzz = cbuffer.data(pg_geom_11_off + 278 * ccomps * dcomps);

            auto g_z_x_x_xzzz = cbuffer.data(pg_geom_11_off + 279 * ccomps * dcomps);

            auto g_z_x_x_yyyy = cbuffer.data(pg_geom_11_off + 280 * ccomps * dcomps);

            auto g_z_x_x_yyyz = cbuffer.data(pg_geom_11_off + 281 * ccomps * dcomps);

            auto g_z_x_x_yyzz = cbuffer.data(pg_geom_11_off + 282 * ccomps * dcomps);

            auto g_z_x_x_yzzz = cbuffer.data(pg_geom_11_off + 283 * ccomps * dcomps);

            auto g_z_x_x_zzzz = cbuffer.data(pg_geom_11_off + 284 * ccomps * dcomps);

            auto g_z_x_y_xxxx = cbuffer.data(pg_geom_11_off + 285 * ccomps * dcomps);

            auto g_z_x_y_xxxy = cbuffer.data(pg_geom_11_off + 286 * ccomps * dcomps);

            auto g_z_x_y_xxxz = cbuffer.data(pg_geom_11_off + 287 * ccomps * dcomps);

            auto g_z_x_y_xxyy = cbuffer.data(pg_geom_11_off + 288 * ccomps * dcomps);

            auto g_z_x_y_xxyz = cbuffer.data(pg_geom_11_off + 289 * ccomps * dcomps);

            auto g_z_x_y_xxzz = cbuffer.data(pg_geom_11_off + 290 * ccomps * dcomps);

            auto g_z_x_y_xyyy = cbuffer.data(pg_geom_11_off + 291 * ccomps * dcomps);

            auto g_z_x_y_xyyz = cbuffer.data(pg_geom_11_off + 292 * ccomps * dcomps);

            auto g_z_x_y_xyzz = cbuffer.data(pg_geom_11_off + 293 * ccomps * dcomps);

            auto g_z_x_y_xzzz = cbuffer.data(pg_geom_11_off + 294 * ccomps * dcomps);

            auto g_z_x_y_yyyy = cbuffer.data(pg_geom_11_off + 295 * ccomps * dcomps);

            auto g_z_x_y_yyyz = cbuffer.data(pg_geom_11_off + 296 * ccomps * dcomps);

            auto g_z_x_y_yyzz = cbuffer.data(pg_geom_11_off + 297 * ccomps * dcomps);

            auto g_z_x_y_yzzz = cbuffer.data(pg_geom_11_off + 298 * ccomps * dcomps);

            auto g_z_x_y_zzzz = cbuffer.data(pg_geom_11_off + 299 * ccomps * dcomps);

            auto g_z_x_z_xxxx = cbuffer.data(pg_geom_11_off + 300 * ccomps * dcomps);

            auto g_z_x_z_xxxy = cbuffer.data(pg_geom_11_off + 301 * ccomps * dcomps);

            auto g_z_x_z_xxxz = cbuffer.data(pg_geom_11_off + 302 * ccomps * dcomps);

            auto g_z_x_z_xxyy = cbuffer.data(pg_geom_11_off + 303 * ccomps * dcomps);

            auto g_z_x_z_xxyz = cbuffer.data(pg_geom_11_off + 304 * ccomps * dcomps);

            auto g_z_x_z_xxzz = cbuffer.data(pg_geom_11_off + 305 * ccomps * dcomps);

            auto g_z_x_z_xyyy = cbuffer.data(pg_geom_11_off + 306 * ccomps * dcomps);

            auto g_z_x_z_xyyz = cbuffer.data(pg_geom_11_off + 307 * ccomps * dcomps);

            auto g_z_x_z_xyzz = cbuffer.data(pg_geom_11_off + 308 * ccomps * dcomps);

            auto g_z_x_z_xzzz = cbuffer.data(pg_geom_11_off + 309 * ccomps * dcomps);

            auto g_z_x_z_yyyy = cbuffer.data(pg_geom_11_off + 310 * ccomps * dcomps);

            auto g_z_x_z_yyyz = cbuffer.data(pg_geom_11_off + 311 * ccomps * dcomps);

            auto g_z_x_z_yyzz = cbuffer.data(pg_geom_11_off + 312 * ccomps * dcomps);

            auto g_z_x_z_yzzz = cbuffer.data(pg_geom_11_off + 313 * ccomps * dcomps);

            auto g_z_x_z_zzzz = cbuffer.data(pg_geom_11_off + 314 * ccomps * dcomps);

            auto g_z_y_x_xxxx = cbuffer.data(pg_geom_11_off + 315 * ccomps * dcomps);

            auto g_z_y_x_xxxy = cbuffer.data(pg_geom_11_off + 316 * ccomps * dcomps);

            auto g_z_y_x_xxxz = cbuffer.data(pg_geom_11_off + 317 * ccomps * dcomps);

            auto g_z_y_x_xxyy = cbuffer.data(pg_geom_11_off + 318 * ccomps * dcomps);

            auto g_z_y_x_xxyz = cbuffer.data(pg_geom_11_off + 319 * ccomps * dcomps);

            auto g_z_y_x_xxzz = cbuffer.data(pg_geom_11_off + 320 * ccomps * dcomps);

            auto g_z_y_x_xyyy = cbuffer.data(pg_geom_11_off + 321 * ccomps * dcomps);

            auto g_z_y_x_xyyz = cbuffer.data(pg_geom_11_off + 322 * ccomps * dcomps);

            auto g_z_y_x_xyzz = cbuffer.data(pg_geom_11_off + 323 * ccomps * dcomps);

            auto g_z_y_x_xzzz = cbuffer.data(pg_geom_11_off + 324 * ccomps * dcomps);

            auto g_z_y_x_yyyy = cbuffer.data(pg_geom_11_off + 325 * ccomps * dcomps);

            auto g_z_y_x_yyyz = cbuffer.data(pg_geom_11_off + 326 * ccomps * dcomps);

            auto g_z_y_x_yyzz = cbuffer.data(pg_geom_11_off + 327 * ccomps * dcomps);

            auto g_z_y_x_yzzz = cbuffer.data(pg_geom_11_off + 328 * ccomps * dcomps);

            auto g_z_y_x_zzzz = cbuffer.data(pg_geom_11_off + 329 * ccomps * dcomps);

            auto g_z_y_y_xxxx = cbuffer.data(pg_geom_11_off + 330 * ccomps * dcomps);

            auto g_z_y_y_xxxy = cbuffer.data(pg_geom_11_off + 331 * ccomps * dcomps);

            auto g_z_y_y_xxxz = cbuffer.data(pg_geom_11_off + 332 * ccomps * dcomps);

            auto g_z_y_y_xxyy = cbuffer.data(pg_geom_11_off + 333 * ccomps * dcomps);

            auto g_z_y_y_xxyz = cbuffer.data(pg_geom_11_off + 334 * ccomps * dcomps);

            auto g_z_y_y_xxzz = cbuffer.data(pg_geom_11_off + 335 * ccomps * dcomps);

            auto g_z_y_y_xyyy = cbuffer.data(pg_geom_11_off + 336 * ccomps * dcomps);

            auto g_z_y_y_xyyz = cbuffer.data(pg_geom_11_off + 337 * ccomps * dcomps);

            auto g_z_y_y_xyzz = cbuffer.data(pg_geom_11_off + 338 * ccomps * dcomps);

            auto g_z_y_y_xzzz = cbuffer.data(pg_geom_11_off + 339 * ccomps * dcomps);

            auto g_z_y_y_yyyy = cbuffer.data(pg_geom_11_off + 340 * ccomps * dcomps);

            auto g_z_y_y_yyyz = cbuffer.data(pg_geom_11_off + 341 * ccomps * dcomps);

            auto g_z_y_y_yyzz = cbuffer.data(pg_geom_11_off + 342 * ccomps * dcomps);

            auto g_z_y_y_yzzz = cbuffer.data(pg_geom_11_off + 343 * ccomps * dcomps);

            auto g_z_y_y_zzzz = cbuffer.data(pg_geom_11_off + 344 * ccomps * dcomps);

            auto g_z_y_z_xxxx = cbuffer.data(pg_geom_11_off + 345 * ccomps * dcomps);

            auto g_z_y_z_xxxy = cbuffer.data(pg_geom_11_off + 346 * ccomps * dcomps);

            auto g_z_y_z_xxxz = cbuffer.data(pg_geom_11_off + 347 * ccomps * dcomps);

            auto g_z_y_z_xxyy = cbuffer.data(pg_geom_11_off + 348 * ccomps * dcomps);

            auto g_z_y_z_xxyz = cbuffer.data(pg_geom_11_off + 349 * ccomps * dcomps);

            auto g_z_y_z_xxzz = cbuffer.data(pg_geom_11_off + 350 * ccomps * dcomps);

            auto g_z_y_z_xyyy = cbuffer.data(pg_geom_11_off + 351 * ccomps * dcomps);

            auto g_z_y_z_xyyz = cbuffer.data(pg_geom_11_off + 352 * ccomps * dcomps);

            auto g_z_y_z_xyzz = cbuffer.data(pg_geom_11_off + 353 * ccomps * dcomps);

            auto g_z_y_z_xzzz = cbuffer.data(pg_geom_11_off + 354 * ccomps * dcomps);

            auto g_z_y_z_yyyy = cbuffer.data(pg_geom_11_off + 355 * ccomps * dcomps);

            auto g_z_y_z_yyyz = cbuffer.data(pg_geom_11_off + 356 * ccomps * dcomps);

            auto g_z_y_z_yyzz = cbuffer.data(pg_geom_11_off + 357 * ccomps * dcomps);

            auto g_z_y_z_yzzz = cbuffer.data(pg_geom_11_off + 358 * ccomps * dcomps);

            auto g_z_y_z_zzzz = cbuffer.data(pg_geom_11_off + 359 * ccomps * dcomps);

            auto g_z_z_x_xxxx = cbuffer.data(pg_geom_11_off + 360 * ccomps * dcomps);

            auto g_z_z_x_xxxy = cbuffer.data(pg_geom_11_off + 361 * ccomps * dcomps);

            auto g_z_z_x_xxxz = cbuffer.data(pg_geom_11_off + 362 * ccomps * dcomps);

            auto g_z_z_x_xxyy = cbuffer.data(pg_geom_11_off + 363 * ccomps * dcomps);

            auto g_z_z_x_xxyz = cbuffer.data(pg_geom_11_off + 364 * ccomps * dcomps);

            auto g_z_z_x_xxzz = cbuffer.data(pg_geom_11_off + 365 * ccomps * dcomps);

            auto g_z_z_x_xyyy = cbuffer.data(pg_geom_11_off + 366 * ccomps * dcomps);

            auto g_z_z_x_xyyz = cbuffer.data(pg_geom_11_off + 367 * ccomps * dcomps);

            auto g_z_z_x_xyzz = cbuffer.data(pg_geom_11_off + 368 * ccomps * dcomps);

            auto g_z_z_x_xzzz = cbuffer.data(pg_geom_11_off + 369 * ccomps * dcomps);

            auto g_z_z_x_yyyy = cbuffer.data(pg_geom_11_off + 370 * ccomps * dcomps);

            auto g_z_z_x_yyyz = cbuffer.data(pg_geom_11_off + 371 * ccomps * dcomps);

            auto g_z_z_x_yyzz = cbuffer.data(pg_geom_11_off + 372 * ccomps * dcomps);

            auto g_z_z_x_yzzz = cbuffer.data(pg_geom_11_off + 373 * ccomps * dcomps);

            auto g_z_z_x_zzzz = cbuffer.data(pg_geom_11_off + 374 * ccomps * dcomps);

            auto g_z_z_y_xxxx = cbuffer.data(pg_geom_11_off + 375 * ccomps * dcomps);

            auto g_z_z_y_xxxy = cbuffer.data(pg_geom_11_off + 376 * ccomps * dcomps);

            auto g_z_z_y_xxxz = cbuffer.data(pg_geom_11_off + 377 * ccomps * dcomps);

            auto g_z_z_y_xxyy = cbuffer.data(pg_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_z_y_xxyz = cbuffer.data(pg_geom_11_off + 379 * ccomps * dcomps);

            auto g_z_z_y_xxzz = cbuffer.data(pg_geom_11_off + 380 * ccomps * dcomps);

            auto g_z_z_y_xyyy = cbuffer.data(pg_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_z_y_xyyz = cbuffer.data(pg_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_z_y_xyzz = cbuffer.data(pg_geom_11_off + 383 * ccomps * dcomps);

            auto g_z_z_y_xzzz = cbuffer.data(pg_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_z_y_yyyy = cbuffer.data(pg_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_z_y_yyyz = cbuffer.data(pg_geom_11_off + 386 * ccomps * dcomps);

            auto g_z_z_y_yyzz = cbuffer.data(pg_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_z_y_yzzz = cbuffer.data(pg_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_z_y_zzzz = cbuffer.data(pg_geom_11_off + 389 * ccomps * dcomps);

            auto g_z_z_z_xxxx = cbuffer.data(pg_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_z_z_xxxy = cbuffer.data(pg_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_z_z_xxxz = cbuffer.data(pg_geom_11_off + 392 * ccomps * dcomps);

            auto g_z_z_z_xxyy = cbuffer.data(pg_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_z_z_xxyz = cbuffer.data(pg_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_z_z_xxzz = cbuffer.data(pg_geom_11_off + 395 * ccomps * dcomps);

            auto g_z_z_z_xyyy = cbuffer.data(pg_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_z_z_xyyz = cbuffer.data(pg_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_z_z_xyzz = cbuffer.data(pg_geom_11_off + 398 * ccomps * dcomps);

            auto g_z_z_z_xzzz = cbuffer.data(pg_geom_11_off + 399 * ccomps * dcomps);

            auto g_z_z_z_yyyy = cbuffer.data(pg_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_z_z_yyyz = cbuffer.data(pg_geom_11_off + 401 * ccomps * dcomps);

            auto g_z_z_z_yyzz = cbuffer.data(pg_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_z_z_yzzz = cbuffer.data(pg_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_z_z_zzzz = cbuffer.data(pg_geom_11_off + 404 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dfxx

            const auto df_geom_11_off = idx_geom_11_dfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_x_xx_xxx = cbuffer.data(df_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_xxy = cbuffer.data(df_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_xxz = cbuffer.data(df_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xx_xyy = cbuffer.data(df_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xx_xyz = cbuffer.data(df_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xx_xzz = cbuffer.data(df_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xx_yyy = cbuffer.data(df_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xx_yyz = cbuffer.data(df_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xx_yzz = cbuffer.data(df_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xx_zzz = cbuffer.data(df_geom_11_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxx, g_0_x_x_xxy, g_0_x_x_xxz, g_0_x_x_xyy, g_0_x_x_xyz, g_0_x_x_xzz, g_0_x_x_yyy, g_0_x_x_yyz, g_0_x_x_yzz, g_0_x_x_zzz, g_x_0_x_xxx, g_x_0_x_xxy, g_x_0_x_xxz, g_x_0_x_xyy, g_x_0_x_xyz, g_x_0_x_xzz, g_x_0_x_yyy, g_x_0_x_yyz, g_x_0_x_yzz, g_x_0_x_zzz, g_x_x_x_xxx, g_x_x_x_xxxx, g_x_x_x_xxxy, g_x_x_x_xxxz, g_x_x_x_xxy, g_x_x_x_xxyy, g_x_x_x_xxyz, g_x_x_x_xxz, g_x_x_x_xxzz, g_x_x_x_xyy, g_x_x_x_xyyy, g_x_x_x_xyyz, g_x_x_x_xyz, g_x_x_x_xyzz, g_x_x_x_xzz, g_x_x_x_xzzz, g_x_x_x_yyy, g_x_x_x_yyz, g_x_x_x_yzz, g_x_x_x_zzz, g_x_x_xx_xxx, g_x_x_xx_xxy, g_x_x_xx_xxz, g_x_x_xx_xyy, g_x_x_xx_xyz, g_x_x_xx_xzz, g_x_x_xx_yyy, g_x_x_xx_yyz, g_x_x_xx_yzz, g_x_x_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xx_xxx[k] = -g_0_x_x_xxx[k] + g_x_0_x_xxx[k] - g_x_x_x_xxx[k] * ab_x + g_x_x_x_xxxx[k];

                g_x_x_xx_xxy[k] = -g_0_x_x_xxy[k] + g_x_0_x_xxy[k] - g_x_x_x_xxy[k] * ab_x + g_x_x_x_xxxy[k];

                g_x_x_xx_xxz[k] = -g_0_x_x_xxz[k] + g_x_0_x_xxz[k] - g_x_x_x_xxz[k] * ab_x + g_x_x_x_xxxz[k];

                g_x_x_xx_xyy[k] = -g_0_x_x_xyy[k] + g_x_0_x_xyy[k] - g_x_x_x_xyy[k] * ab_x + g_x_x_x_xxyy[k];

                g_x_x_xx_xyz[k] = -g_0_x_x_xyz[k] + g_x_0_x_xyz[k] - g_x_x_x_xyz[k] * ab_x + g_x_x_x_xxyz[k];

                g_x_x_xx_xzz[k] = -g_0_x_x_xzz[k] + g_x_0_x_xzz[k] - g_x_x_x_xzz[k] * ab_x + g_x_x_x_xxzz[k];

                g_x_x_xx_yyy[k] = -g_0_x_x_yyy[k] + g_x_0_x_yyy[k] - g_x_x_x_yyy[k] * ab_x + g_x_x_x_xyyy[k];

                g_x_x_xx_yyz[k] = -g_0_x_x_yyz[k] + g_x_0_x_yyz[k] - g_x_x_x_yyz[k] * ab_x + g_x_x_x_xyyz[k];

                g_x_x_xx_yzz[k] = -g_0_x_x_yzz[k] + g_x_0_x_yzz[k] - g_x_x_x_yzz[k] * ab_x + g_x_x_x_xyzz[k];

                g_x_x_xx_zzz[k] = -g_0_x_x_zzz[k] + g_x_0_x_zzz[k] - g_x_x_x_zzz[k] * ab_x + g_x_x_x_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_x_xy_xxx = cbuffer.data(df_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xy_xxy = cbuffer.data(df_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xy_xxz = cbuffer.data(df_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xy_xyy = cbuffer.data(df_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xy_xyz = cbuffer.data(df_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xy_xzz = cbuffer.data(df_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xy_yyy = cbuffer.data(df_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xy_yyz = cbuffer.data(df_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xy_yzz = cbuffer.data(df_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xy_zzz = cbuffer.data(df_geom_11_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xxx, g_x_x_x_xxxy, g_x_x_x_xxy, g_x_x_x_xxyy, g_x_x_x_xxyz, g_x_x_x_xxz, g_x_x_x_xyy, g_x_x_x_xyyy, g_x_x_x_xyyz, g_x_x_x_xyz, g_x_x_x_xyzz, g_x_x_x_xzz, g_x_x_x_yyy, g_x_x_x_yyyy, g_x_x_x_yyyz, g_x_x_x_yyz, g_x_x_x_yyzz, g_x_x_x_yzz, g_x_x_x_yzzz, g_x_x_x_zzz, g_x_x_xy_xxx, g_x_x_xy_xxy, g_x_x_xy_xxz, g_x_x_xy_xyy, g_x_x_xy_xyz, g_x_x_xy_xzz, g_x_x_xy_yyy, g_x_x_xy_yyz, g_x_x_xy_yzz, g_x_x_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xy_xxx[k] = -g_x_x_x_xxx[k] * ab_y + g_x_x_x_xxxy[k];

                g_x_x_xy_xxy[k] = -g_x_x_x_xxy[k] * ab_y + g_x_x_x_xxyy[k];

                g_x_x_xy_xxz[k] = -g_x_x_x_xxz[k] * ab_y + g_x_x_x_xxyz[k];

                g_x_x_xy_xyy[k] = -g_x_x_x_xyy[k] * ab_y + g_x_x_x_xyyy[k];

                g_x_x_xy_xyz[k] = -g_x_x_x_xyz[k] * ab_y + g_x_x_x_xyyz[k];

                g_x_x_xy_xzz[k] = -g_x_x_x_xzz[k] * ab_y + g_x_x_x_xyzz[k];

                g_x_x_xy_yyy[k] = -g_x_x_x_yyy[k] * ab_y + g_x_x_x_yyyy[k];

                g_x_x_xy_yyz[k] = -g_x_x_x_yyz[k] * ab_y + g_x_x_x_yyyz[k];

                g_x_x_xy_yzz[k] = -g_x_x_x_yzz[k] * ab_y + g_x_x_x_yyzz[k];

                g_x_x_xy_zzz[k] = -g_x_x_x_zzz[k] * ab_y + g_x_x_x_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_xz_xxx = cbuffer.data(df_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xz_xxy = cbuffer.data(df_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xz_xxz = cbuffer.data(df_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xz_xyy = cbuffer.data(df_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xz_xyz = cbuffer.data(df_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xz_xzz = cbuffer.data(df_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xz_yyy = cbuffer.data(df_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xz_yyz = cbuffer.data(df_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xz_yzz = cbuffer.data(df_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xz_zzz = cbuffer.data(df_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xxx, g_x_x_x_xxxz, g_x_x_x_xxy, g_x_x_x_xxyz, g_x_x_x_xxz, g_x_x_x_xxzz, g_x_x_x_xyy, g_x_x_x_xyyz, g_x_x_x_xyz, g_x_x_x_xyzz, g_x_x_x_xzz, g_x_x_x_xzzz, g_x_x_x_yyy, g_x_x_x_yyyz, g_x_x_x_yyz, g_x_x_x_yyzz, g_x_x_x_yzz, g_x_x_x_yzzz, g_x_x_x_zzz, g_x_x_x_zzzz, g_x_x_xz_xxx, g_x_x_xz_xxy, g_x_x_xz_xxz, g_x_x_xz_xyy, g_x_x_xz_xyz, g_x_x_xz_xzz, g_x_x_xz_yyy, g_x_x_xz_yyz, g_x_x_xz_yzz, g_x_x_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xz_xxx[k] = -g_x_x_x_xxx[k] * ab_z + g_x_x_x_xxxz[k];

                g_x_x_xz_xxy[k] = -g_x_x_x_xxy[k] * ab_z + g_x_x_x_xxyz[k];

                g_x_x_xz_xxz[k] = -g_x_x_x_xxz[k] * ab_z + g_x_x_x_xxzz[k];

                g_x_x_xz_xyy[k] = -g_x_x_x_xyy[k] * ab_z + g_x_x_x_xyyz[k];

                g_x_x_xz_xyz[k] = -g_x_x_x_xyz[k] * ab_z + g_x_x_x_xyzz[k];

                g_x_x_xz_xzz[k] = -g_x_x_x_xzz[k] * ab_z + g_x_x_x_xzzz[k];

                g_x_x_xz_yyy[k] = -g_x_x_x_yyy[k] * ab_z + g_x_x_x_yyyz[k];

                g_x_x_xz_yyz[k] = -g_x_x_x_yyz[k] * ab_z + g_x_x_x_yyzz[k];

                g_x_x_xz_yzz[k] = -g_x_x_x_yzz[k] * ab_z + g_x_x_x_yzzz[k];

                g_x_x_xz_zzz[k] = -g_x_x_x_zzz[k] * ab_z + g_x_x_x_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_x_x_yy_xxx = cbuffer.data(df_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_yy_xxy = cbuffer.data(df_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_yy_xxz = cbuffer.data(df_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_yy_xyy = cbuffer.data(df_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_yy_xyz = cbuffer.data(df_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_yy_xzz = cbuffer.data(df_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_yy_yyy = cbuffer.data(df_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_yy_yyz = cbuffer.data(df_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_yy_yzz = cbuffer.data(df_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_yy_zzz = cbuffer.data(df_geom_11_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_y_xxx, g_x_x_y_xxxy, g_x_x_y_xxy, g_x_x_y_xxyy, g_x_x_y_xxyz, g_x_x_y_xxz, g_x_x_y_xyy, g_x_x_y_xyyy, g_x_x_y_xyyz, g_x_x_y_xyz, g_x_x_y_xyzz, g_x_x_y_xzz, g_x_x_y_yyy, g_x_x_y_yyyy, g_x_x_y_yyyz, g_x_x_y_yyz, g_x_x_y_yyzz, g_x_x_y_yzz, g_x_x_y_yzzz, g_x_x_y_zzz, g_x_x_yy_xxx, g_x_x_yy_xxy, g_x_x_yy_xxz, g_x_x_yy_xyy, g_x_x_yy_xyz, g_x_x_yy_xzz, g_x_x_yy_yyy, g_x_x_yy_yyz, g_x_x_yy_yzz, g_x_x_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yy_xxx[k] = -g_x_x_y_xxx[k] * ab_y + g_x_x_y_xxxy[k];

                g_x_x_yy_xxy[k] = -g_x_x_y_xxy[k] * ab_y + g_x_x_y_xxyy[k];

                g_x_x_yy_xxz[k] = -g_x_x_y_xxz[k] * ab_y + g_x_x_y_xxyz[k];

                g_x_x_yy_xyy[k] = -g_x_x_y_xyy[k] * ab_y + g_x_x_y_xyyy[k];

                g_x_x_yy_xyz[k] = -g_x_x_y_xyz[k] * ab_y + g_x_x_y_xyyz[k];

                g_x_x_yy_xzz[k] = -g_x_x_y_xzz[k] * ab_y + g_x_x_y_xyzz[k];

                g_x_x_yy_yyy[k] = -g_x_x_y_yyy[k] * ab_y + g_x_x_y_yyyy[k];

                g_x_x_yy_yyz[k] = -g_x_x_y_yyz[k] * ab_y + g_x_x_y_yyyz[k];

                g_x_x_yy_yzz[k] = -g_x_x_y_yzz[k] * ab_y + g_x_x_y_yyzz[k];

                g_x_x_yy_zzz[k] = -g_x_x_y_zzz[k] * ab_y + g_x_x_y_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_x_x_yz_xxx = cbuffer.data(df_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_yz_xxy = cbuffer.data(df_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_yz_xxz = cbuffer.data(df_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_yz_xyy = cbuffer.data(df_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_yz_xyz = cbuffer.data(df_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_yz_xzz = cbuffer.data(df_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_yz_yyy = cbuffer.data(df_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_yz_yyz = cbuffer.data(df_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_yz_yzz = cbuffer.data(df_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_yz_zzz = cbuffer.data(df_geom_11_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yz_xxx, g_x_x_yz_xxy, g_x_x_yz_xxz, g_x_x_yz_xyy, g_x_x_yz_xyz, g_x_x_yz_xzz, g_x_x_yz_yyy, g_x_x_yz_yyz, g_x_x_yz_yzz, g_x_x_yz_zzz, g_x_x_z_xxx, g_x_x_z_xxxy, g_x_x_z_xxy, g_x_x_z_xxyy, g_x_x_z_xxyz, g_x_x_z_xxz, g_x_x_z_xyy, g_x_x_z_xyyy, g_x_x_z_xyyz, g_x_x_z_xyz, g_x_x_z_xyzz, g_x_x_z_xzz, g_x_x_z_yyy, g_x_x_z_yyyy, g_x_x_z_yyyz, g_x_x_z_yyz, g_x_x_z_yyzz, g_x_x_z_yzz, g_x_x_z_yzzz, g_x_x_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yz_xxx[k] = -g_x_x_z_xxx[k] * ab_y + g_x_x_z_xxxy[k];

                g_x_x_yz_xxy[k] = -g_x_x_z_xxy[k] * ab_y + g_x_x_z_xxyy[k];

                g_x_x_yz_xxz[k] = -g_x_x_z_xxz[k] * ab_y + g_x_x_z_xxyz[k];

                g_x_x_yz_xyy[k] = -g_x_x_z_xyy[k] * ab_y + g_x_x_z_xyyy[k];

                g_x_x_yz_xyz[k] = -g_x_x_z_xyz[k] * ab_y + g_x_x_z_xyyz[k];

                g_x_x_yz_xzz[k] = -g_x_x_z_xzz[k] * ab_y + g_x_x_z_xyzz[k];

                g_x_x_yz_yyy[k] = -g_x_x_z_yyy[k] * ab_y + g_x_x_z_yyyy[k];

                g_x_x_yz_yyz[k] = -g_x_x_z_yyz[k] * ab_y + g_x_x_z_yyyz[k];

                g_x_x_yz_yzz[k] = -g_x_x_z_yzz[k] * ab_y + g_x_x_z_yyzz[k];

                g_x_x_yz_zzz[k] = -g_x_x_z_zzz[k] * ab_y + g_x_x_z_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_x_x_zz_xxx = cbuffer.data(df_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_zz_xxy = cbuffer.data(df_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_zz_xxz = cbuffer.data(df_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_zz_xyy = cbuffer.data(df_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_zz_xyz = cbuffer.data(df_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_zz_xzz = cbuffer.data(df_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_zz_yyy = cbuffer.data(df_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_zz_yyz = cbuffer.data(df_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_zz_yzz = cbuffer.data(df_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_zz_zzz = cbuffer.data(df_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_z_xxx, g_x_x_z_xxxz, g_x_x_z_xxy, g_x_x_z_xxyz, g_x_x_z_xxz, g_x_x_z_xxzz, g_x_x_z_xyy, g_x_x_z_xyyz, g_x_x_z_xyz, g_x_x_z_xyzz, g_x_x_z_xzz, g_x_x_z_xzzz, g_x_x_z_yyy, g_x_x_z_yyyz, g_x_x_z_yyz, g_x_x_z_yyzz, g_x_x_z_yzz, g_x_x_z_yzzz, g_x_x_z_zzz, g_x_x_z_zzzz, g_x_x_zz_xxx, g_x_x_zz_xxy, g_x_x_zz_xxz, g_x_x_zz_xyy, g_x_x_zz_xyz, g_x_x_zz_xzz, g_x_x_zz_yyy, g_x_x_zz_yyz, g_x_x_zz_yzz, g_x_x_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zz_xxx[k] = -g_x_x_z_xxx[k] * ab_z + g_x_x_z_xxxz[k];

                g_x_x_zz_xxy[k] = -g_x_x_z_xxy[k] * ab_z + g_x_x_z_xxyz[k];

                g_x_x_zz_xxz[k] = -g_x_x_z_xxz[k] * ab_z + g_x_x_z_xxzz[k];

                g_x_x_zz_xyy[k] = -g_x_x_z_xyy[k] * ab_z + g_x_x_z_xyyz[k];

                g_x_x_zz_xyz[k] = -g_x_x_z_xyz[k] * ab_z + g_x_x_z_xyzz[k];

                g_x_x_zz_xzz[k] = -g_x_x_z_xzz[k] * ab_z + g_x_x_z_xzzz[k];

                g_x_x_zz_yyy[k] = -g_x_x_z_yyy[k] * ab_z + g_x_x_z_yyyz[k];

                g_x_x_zz_yyz[k] = -g_x_x_z_yyz[k] * ab_z + g_x_x_z_yyzz[k];

                g_x_x_zz_yzz[k] = -g_x_x_z_yzz[k] * ab_z + g_x_x_z_yzzz[k];

                g_x_x_zz_zzz[k] = -g_x_x_z_zzz[k] * ab_z + g_x_x_z_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_x_y_xx_xxx = cbuffer.data(df_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_xx_xxy = cbuffer.data(df_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_xx_xxz = cbuffer.data(df_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_xx_xyy = cbuffer.data(df_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_xx_xyz = cbuffer.data(df_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_xx_xzz = cbuffer.data(df_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_xx_yyy = cbuffer.data(df_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_xx_yyz = cbuffer.data(df_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_xx_yzz = cbuffer.data(df_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_xx_zzz = cbuffer.data(df_geom_11_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_x_xxx, g_0_y_x_xxy, g_0_y_x_xxz, g_0_y_x_xyy, g_0_y_x_xyz, g_0_y_x_xzz, g_0_y_x_yyy, g_0_y_x_yyz, g_0_y_x_yzz, g_0_y_x_zzz, g_x_y_x_xxx, g_x_y_x_xxxx, g_x_y_x_xxxy, g_x_y_x_xxxz, g_x_y_x_xxy, g_x_y_x_xxyy, g_x_y_x_xxyz, g_x_y_x_xxz, g_x_y_x_xxzz, g_x_y_x_xyy, g_x_y_x_xyyy, g_x_y_x_xyyz, g_x_y_x_xyz, g_x_y_x_xyzz, g_x_y_x_xzz, g_x_y_x_xzzz, g_x_y_x_yyy, g_x_y_x_yyz, g_x_y_x_yzz, g_x_y_x_zzz, g_x_y_xx_xxx, g_x_y_xx_xxy, g_x_y_xx_xxz, g_x_y_xx_xyy, g_x_y_xx_xyz, g_x_y_xx_xzz, g_x_y_xx_yyy, g_x_y_xx_yyz, g_x_y_xx_yzz, g_x_y_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xx_xxx[k] = -g_0_y_x_xxx[k] - g_x_y_x_xxx[k] * ab_x + g_x_y_x_xxxx[k];

                g_x_y_xx_xxy[k] = -g_0_y_x_xxy[k] - g_x_y_x_xxy[k] * ab_x + g_x_y_x_xxxy[k];

                g_x_y_xx_xxz[k] = -g_0_y_x_xxz[k] - g_x_y_x_xxz[k] * ab_x + g_x_y_x_xxxz[k];

                g_x_y_xx_xyy[k] = -g_0_y_x_xyy[k] - g_x_y_x_xyy[k] * ab_x + g_x_y_x_xxyy[k];

                g_x_y_xx_xyz[k] = -g_0_y_x_xyz[k] - g_x_y_x_xyz[k] * ab_x + g_x_y_x_xxyz[k];

                g_x_y_xx_xzz[k] = -g_0_y_x_xzz[k] - g_x_y_x_xzz[k] * ab_x + g_x_y_x_xxzz[k];

                g_x_y_xx_yyy[k] = -g_0_y_x_yyy[k] - g_x_y_x_yyy[k] * ab_x + g_x_y_x_xyyy[k];

                g_x_y_xx_yyz[k] = -g_0_y_x_yyz[k] - g_x_y_x_yyz[k] * ab_x + g_x_y_x_xyyz[k];

                g_x_y_xx_yzz[k] = -g_0_y_x_yzz[k] - g_x_y_x_yzz[k] * ab_x + g_x_y_x_xyzz[k];

                g_x_y_xx_zzz[k] = -g_0_y_x_zzz[k] - g_x_y_x_zzz[k] * ab_x + g_x_y_x_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_x_y_xy_xxx = cbuffer.data(df_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_xy_xxy = cbuffer.data(df_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_y_xy_xxz = cbuffer.data(df_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_xy_xyy = cbuffer.data(df_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_xy_xyz = cbuffer.data(df_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_y_xy_xzz = cbuffer.data(df_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_xy_yyy = cbuffer.data(df_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_xy_yyz = cbuffer.data(df_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_y_xy_yzz = cbuffer.data(df_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_xy_zzz = cbuffer.data(df_geom_11_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxx, g_0_y_y_xxy, g_0_y_y_xxz, g_0_y_y_xyy, g_0_y_y_xyz, g_0_y_y_xzz, g_0_y_y_yyy, g_0_y_y_yyz, g_0_y_y_yzz, g_0_y_y_zzz, g_x_y_xy_xxx, g_x_y_xy_xxy, g_x_y_xy_xxz, g_x_y_xy_xyy, g_x_y_xy_xyz, g_x_y_xy_xzz, g_x_y_xy_yyy, g_x_y_xy_yyz, g_x_y_xy_yzz, g_x_y_xy_zzz, g_x_y_y_xxx, g_x_y_y_xxxx, g_x_y_y_xxxy, g_x_y_y_xxxz, g_x_y_y_xxy, g_x_y_y_xxyy, g_x_y_y_xxyz, g_x_y_y_xxz, g_x_y_y_xxzz, g_x_y_y_xyy, g_x_y_y_xyyy, g_x_y_y_xyyz, g_x_y_y_xyz, g_x_y_y_xyzz, g_x_y_y_xzz, g_x_y_y_xzzz, g_x_y_y_yyy, g_x_y_y_yyz, g_x_y_y_yzz, g_x_y_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xy_xxx[k] = -g_0_y_y_xxx[k] - g_x_y_y_xxx[k] * ab_x + g_x_y_y_xxxx[k];

                g_x_y_xy_xxy[k] = -g_0_y_y_xxy[k] - g_x_y_y_xxy[k] * ab_x + g_x_y_y_xxxy[k];

                g_x_y_xy_xxz[k] = -g_0_y_y_xxz[k] - g_x_y_y_xxz[k] * ab_x + g_x_y_y_xxxz[k];

                g_x_y_xy_xyy[k] = -g_0_y_y_xyy[k] - g_x_y_y_xyy[k] * ab_x + g_x_y_y_xxyy[k];

                g_x_y_xy_xyz[k] = -g_0_y_y_xyz[k] - g_x_y_y_xyz[k] * ab_x + g_x_y_y_xxyz[k];

                g_x_y_xy_xzz[k] = -g_0_y_y_xzz[k] - g_x_y_y_xzz[k] * ab_x + g_x_y_y_xxzz[k];

                g_x_y_xy_yyy[k] = -g_0_y_y_yyy[k] - g_x_y_y_yyy[k] * ab_x + g_x_y_y_xyyy[k];

                g_x_y_xy_yyz[k] = -g_0_y_y_yyz[k] - g_x_y_y_yyz[k] * ab_x + g_x_y_y_xyyz[k];

                g_x_y_xy_yzz[k] = -g_0_y_y_yzz[k] - g_x_y_y_yzz[k] * ab_x + g_x_y_y_xyzz[k];

                g_x_y_xy_zzz[k] = -g_0_y_y_zzz[k] - g_x_y_y_zzz[k] * ab_x + g_x_y_y_xzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_x_y_xz_xxx = cbuffer.data(df_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_y_xz_xxy = cbuffer.data(df_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_xz_xxz = cbuffer.data(df_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_xz_xyy = cbuffer.data(df_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_y_xz_xyz = cbuffer.data(df_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_xz_xzz = cbuffer.data(df_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_xz_yyy = cbuffer.data(df_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_xz_yyz = cbuffer.data(df_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_xz_yzz = cbuffer.data(df_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_xz_zzz = cbuffer.data(df_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_x_xxx, g_x_y_x_xxxz, g_x_y_x_xxy, g_x_y_x_xxyz, g_x_y_x_xxz, g_x_y_x_xxzz, g_x_y_x_xyy, g_x_y_x_xyyz, g_x_y_x_xyz, g_x_y_x_xyzz, g_x_y_x_xzz, g_x_y_x_xzzz, g_x_y_x_yyy, g_x_y_x_yyyz, g_x_y_x_yyz, g_x_y_x_yyzz, g_x_y_x_yzz, g_x_y_x_yzzz, g_x_y_x_zzz, g_x_y_x_zzzz, g_x_y_xz_xxx, g_x_y_xz_xxy, g_x_y_xz_xxz, g_x_y_xz_xyy, g_x_y_xz_xyz, g_x_y_xz_xzz, g_x_y_xz_yyy, g_x_y_xz_yyz, g_x_y_xz_yzz, g_x_y_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xz_xxx[k] = -g_x_y_x_xxx[k] * ab_z + g_x_y_x_xxxz[k];

                g_x_y_xz_xxy[k] = -g_x_y_x_xxy[k] * ab_z + g_x_y_x_xxyz[k];

                g_x_y_xz_xxz[k] = -g_x_y_x_xxz[k] * ab_z + g_x_y_x_xxzz[k];

                g_x_y_xz_xyy[k] = -g_x_y_x_xyy[k] * ab_z + g_x_y_x_xyyz[k];

                g_x_y_xz_xyz[k] = -g_x_y_x_xyz[k] * ab_z + g_x_y_x_xyzz[k];

                g_x_y_xz_xzz[k] = -g_x_y_x_xzz[k] * ab_z + g_x_y_x_xzzz[k];

                g_x_y_xz_yyy[k] = -g_x_y_x_yyy[k] * ab_z + g_x_y_x_yyyz[k];

                g_x_y_xz_yyz[k] = -g_x_y_x_yyz[k] * ab_z + g_x_y_x_yyzz[k];

                g_x_y_xz_yzz[k] = -g_x_y_x_yzz[k] * ab_z + g_x_y_x_yzzz[k];

                g_x_y_xz_zzz[k] = -g_x_y_x_zzz[k] * ab_z + g_x_y_x_zzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_x_y_yy_xxx = cbuffer.data(df_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_yy_xxy = cbuffer.data(df_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_yy_xxz = cbuffer.data(df_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_yy_xyy = cbuffer.data(df_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_yy_xyz = cbuffer.data(df_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_yy_xzz = cbuffer.data(df_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_yy_yyy = cbuffer.data(df_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_yy_yyz = cbuffer.data(df_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_yy_yzz = cbuffer.data(df_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_yy_zzz = cbuffer.data(df_geom_11_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxx, g_x_0_y_xxy, g_x_0_y_xxz, g_x_0_y_xyy, g_x_0_y_xyz, g_x_0_y_xzz, g_x_0_y_yyy, g_x_0_y_yyz, g_x_0_y_yzz, g_x_0_y_zzz, g_x_y_y_xxx, g_x_y_y_xxxy, g_x_y_y_xxy, g_x_y_y_xxyy, g_x_y_y_xxyz, g_x_y_y_xxz, g_x_y_y_xyy, g_x_y_y_xyyy, g_x_y_y_xyyz, g_x_y_y_xyz, g_x_y_y_xyzz, g_x_y_y_xzz, g_x_y_y_yyy, g_x_y_y_yyyy, g_x_y_y_yyyz, g_x_y_y_yyz, g_x_y_y_yyzz, g_x_y_y_yzz, g_x_y_y_yzzz, g_x_y_y_zzz, g_x_y_yy_xxx, g_x_y_yy_xxy, g_x_y_yy_xxz, g_x_y_yy_xyy, g_x_y_yy_xyz, g_x_y_yy_xzz, g_x_y_yy_yyy, g_x_y_yy_yyz, g_x_y_yy_yzz, g_x_y_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yy_xxx[k] = g_x_0_y_xxx[k] - g_x_y_y_xxx[k] * ab_y + g_x_y_y_xxxy[k];

                g_x_y_yy_xxy[k] = g_x_0_y_xxy[k] - g_x_y_y_xxy[k] * ab_y + g_x_y_y_xxyy[k];

                g_x_y_yy_xxz[k] = g_x_0_y_xxz[k] - g_x_y_y_xxz[k] * ab_y + g_x_y_y_xxyz[k];

                g_x_y_yy_xyy[k] = g_x_0_y_xyy[k] - g_x_y_y_xyy[k] * ab_y + g_x_y_y_xyyy[k];

                g_x_y_yy_xyz[k] = g_x_0_y_xyz[k] - g_x_y_y_xyz[k] * ab_y + g_x_y_y_xyyz[k];

                g_x_y_yy_xzz[k] = g_x_0_y_xzz[k] - g_x_y_y_xzz[k] * ab_y + g_x_y_y_xyzz[k];

                g_x_y_yy_yyy[k] = g_x_0_y_yyy[k] - g_x_y_y_yyy[k] * ab_y + g_x_y_y_yyyy[k];

                g_x_y_yy_yyz[k] = g_x_0_y_yyz[k] - g_x_y_y_yyz[k] * ab_y + g_x_y_y_yyyz[k];

                g_x_y_yy_yzz[k] = g_x_0_y_yzz[k] - g_x_y_y_yzz[k] * ab_y + g_x_y_y_yyzz[k];

                g_x_y_yy_zzz[k] = g_x_0_y_zzz[k] - g_x_y_y_zzz[k] * ab_y + g_x_y_y_yzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_x_y_yz_xxx = cbuffer.data(df_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_yz_xxy = cbuffer.data(df_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_yz_xxz = cbuffer.data(df_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_yz_xyy = cbuffer.data(df_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_yz_xyz = cbuffer.data(df_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_yz_xzz = cbuffer.data(df_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_yz_yyy = cbuffer.data(df_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_yz_yyz = cbuffer.data(df_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_yz_yzz = cbuffer.data(df_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_yz_zzz = cbuffer.data(df_geom_11_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_y_xxx, g_x_y_y_xxxz, g_x_y_y_xxy, g_x_y_y_xxyz, g_x_y_y_xxz, g_x_y_y_xxzz, g_x_y_y_xyy, g_x_y_y_xyyz, g_x_y_y_xyz, g_x_y_y_xyzz, g_x_y_y_xzz, g_x_y_y_xzzz, g_x_y_y_yyy, g_x_y_y_yyyz, g_x_y_y_yyz, g_x_y_y_yyzz, g_x_y_y_yzz, g_x_y_y_yzzz, g_x_y_y_zzz, g_x_y_y_zzzz, g_x_y_yz_xxx, g_x_y_yz_xxy, g_x_y_yz_xxz, g_x_y_yz_xyy, g_x_y_yz_xyz, g_x_y_yz_xzz, g_x_y_yz_yyy, g_x_y_yz_yyz, g_x_y_yz_yzz, g_x_y_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yz_xxx[k] = -g_x_y_y_xxx[k] * ab_z + g_x_y_y_xxxz[k];

                g_x_y_yz_xxy[k] = -g_x_y_y_xxy[k] * ab_z + g_x_y_y_xxyz[k];

                g_x_y_yz_xxz[k] = -g_x_y_y_xxz[k] * ab_z + g_x_y_y_xxzz[k];

                g_x_y_yz_xyy[k] = -g_x_y_y_xyy[k] * ab_z + g_x_y_y_xyyz[k];

                g_x_y_yz_xyz[k] = -g_x_y_y_xyz[k] * ab_z + g_x_y_y_xyzz[k];

                g_x_y_yz_xzz[k] = -g_x_y_y_xzz[k] * ab_z + g_x_y_y_xzzz[k];

                g_x_y_yz_yyy[k] = -g_x_y_y_yyy[k] * ab_z + g_x_y_y_yyyz[k];

                g_x_y_yz_yyz[k] = -g_x_y_y_yyz[k] * ab_z + g_x_y_y_yyzz[k];

                g_x_y_yz_yzz[k] = -g_x_y_y_yzz[k] * ab_z + g_x_y_y_yzzz[k];

                g_x_y_yz_zzz[k] = -g_x_y_y_zzz[k] * ab_z + g_x_y_y_zzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_x_y_zz_xxx = cbuffer.data(df_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_zz_xxy = cbuffer.data(df_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_zz_xxz = cbuffer.data(df_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_zz_xyy = cbuffer.data(df_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_zz_xyz = cbuffer.data(df_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_zz_xzz = cbuffer.data(df_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_zz_yyy = cbuffer.data(df_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_zz_yyz = cbuffer.data(df_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_zz_yzz = cbuffer.data(df_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_zz_zzz = cbuffer.data(df_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_z_xxx, g_x_y_z_xxxz, g_x_y_z_xxy, g_x_y_z_xxyz, g_x_y_z_xxz, g_x_y_z_xxzz, g_x_y_z_xyy, g_x_y_z_xyyz, g_x_y_z_xyz, g_x_y_z_xyzz, g_x_y_z_xzz, g_x_y_z_xzzz, g_x_y_z_yyy, g_x_y_z_yyyz, g_x_y_z_yyz, g_x_y_z_yyzz, g_x_y_z_yzz, g_x_y_z_yzzz, g_x_y_z_zzz, g_x_y_z_zzzz, g_x_y_zz_xxx, g_x_y_zz_xxy, g_x_y_zz_xxz, g_x_y_zz_xyy, g_x_y_zz_xyz, g_x_y_zz_xzz, g_x_y_zz_yyy, g_x_y_zz_yyz, g_x_y_zz_yzz, g_x_y_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zz_xxx[k] = -g_x_y_z_xxx[k] * ab_z + g_x_y_z_xxxz[k];

                g_x_y_zz_xxy[k] = -g_x_y_z_xxy[k] * ab_z + g_x_y_z_xxyz[k];

                g_x_y_zz_xxz[k] = -g_x_y_z_xxz[k] * ab_z + g_x_y_z_xxzz[k];

                g_x_y_zz_xyy[k] = -g_x_y_z_xyy[k] * ab_z + g_x_y_z_xyyz[k];

                g_x_y_zz_xyz[k] = -g_x_y_z_xyz[k] * ab_z + g_x_y_z_xyzz[k];

                g_x_y_zz_xzz[k] = -g_x_y_z_xzz[k] * ab_z + g_x_y_z_xzzz[k];

                g_x_y_zz_yyy[k] = -g_x_y_z_yyy[k] * ab_z + g_x_y_z_yyyz[k];

                g_x_y_zz_yyz[k] = -g_x_y_z_yyz[k] * ab_z + g_x_y_z_yyzz[k];

                g_x_y_zz_yzz[k] = -g_x_y_z_yzz[k] * ab_z + g_x_y_z_yzzz[k];

                g_x_y_zz_zzz[k] = -g_x_y_z_zzz[k] * ab_z + g_x_y_z_zzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_x_z_xx_xxx = cbuffer.data(df_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_z_xx_xxy = cbuffer.data(df_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_z_xx_xxz = cbuffer.data(df_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_z_xx_xyy = cbuffer.data(df_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_z_xx_xyz = cbuffer.data(df_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_z_xx_xzz = cbuffer.data(df_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_z_xx_yyy = cbuffer.data(df_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_xx_yyz = cbuffer.data(df_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_xx_yzz = cbuffer.data(df_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_z_xx_zzz = cbuffer.data(df_geom_11_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_x_xxx, g_0_z_x_xxy, g_0_z_x_xxz, g_0_z_x_xyy, g_0_z_x_xyz, g_0_z_x_xzz, g_0_z_x_yyy, g_0_z_x_yyz, g_0_z_x_yzz, g_0_z_x_zzz, g_x_z_x_xxx, g_x_z_x_xxxx, g_x_z_x_xxxy, g_x_z_x_xxxz, g_x_z_x_xxy, g_x_z_x_xxyy, g_x_z_x_xxyz, g_x_z_x_xxz, g_x_z_x_xxzz, g_x_z_x_xyy, g_x_z_x_xyyy, g_x_z_x_xyyz, g_x_z_x_xyz, g_x_z_x_xyzz, g_x_z_x_xzz, g_x_z_x_xzzz, g_x_z_x_yyy, g_x_z_x_yyz, g_x_z_x_yzz, g_x_z_x_zzz, g_x_z_xx_xxx, g_x_z_xx_xxy, g_x_z_xx_xxz, g_x_z_xx_xyy, g_x_z_xx_xyz, g_x_z_xx_xzz, g_x_z_xx_yyy, g_x_z_xx_yyz, g_x_z_xx_yzz, g_x_z_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xx_xxx[k] = -g_0_z_x_xxx[k] - g_x_z_x_xxx[k] * ab_x + g_x_z_x_xxxx[k];

                g_x_z_xx_xxy[k] = -g_0_z_x_xxy[k] - g_x_z_x_xxy[k] * ab_x + g_x_z_x_xxxy[k];

                g_x_z_xx_xxz[k] = -g_0_z_x_xxz[k] - g_x_z_x_xxz[k] * ab_x + g_x_z_x_xxxz[k];

                g_x_z_xx_xyy[k] = -g_0_z_x_xyy[k] - g_x_z_x_xyy[k] * ab_x + g_x_z_x_xxyy[k];

                g_x_z_xx_xyz[k] = -g_0_z_x_xyz[k] - g_x_z_x_xyz[k] * ab_x + g_x_z_x_xxyz[k];

                g_x_z_xx_xzz[k] = -g_0_z_x_xzz[k] - g_x_z_x_xzz[k] * ab_x + g_x_z_x_xxzz[k];

                g_x_z_xx_yyy[k] = -g_0_z_x_yyy[k] - g_x_z_x_yyy[k] * ab_x + g_x_z_x_xyyy[k];

                g_x_z_xx_yyz[k] = -g_0_z_x_yyz[k] - g_x_z_x_yyz[k] * ab_x + g_x_z_x_xyyz[k];

                g_x_z_xx_yzz[k] = -g_0_z_x_yzz[k] - g_x_z_x_yzz[k] * ab_x + g_x_z_x_xyzz[k];

                g_x_z_xx_zzz[k] = -g_0_z_x_zzz[k] - g_x_z_x_zzz[k] * ab_x + g_x_z_x_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_x_z_xy_xxx = cbuffer.data(df_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_xy_xxy = cbuffer.data(df_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_z_xy_xxz = cbuffer.data(df_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_xy_xyy = cbuffer.data(df_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_xy_xyz = cbuffer.data(df_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_z_xy_xzz = cbuffer.data(df_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_z_xy_yyy = cbuffer.data(df_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_z_xy_yyz = cbuffer.data(df_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_z_xy_yzz = cbuffer.data(df_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_z_xy_zzz = cbuffer.data(df_geom_11_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_x_xxx, g_x_z_x_xxxy, g_x_z_x_xxy, g_x_z_x_xxyy, g_x_z_x_xxyz, g_x_z_x_xxz, g_x_z_x_xyy, g_x_z_x_xyyy, g_x_z_x_xyyz, g_x_z_x_xyz, g_x_z_x_xyzz, g_x_z_x_xzz, g_x_z_x_yyy, g_x_z_x_yyyy, g_x_z_x_yyyz, g_x_z_x_yyz, g_x_z_x_yyzz, g_x_z_x_yzz, g_x_z_x_yzzz, g_x_z_x_zzz, g_x_z_xy_xxx, g_x_z_xy_xxy, g_x_z_xy_xxz, g_x_z_xy_xyy, g_x_z_xy_xyz, g_x_z_xy_xzz, g_x_z_xy_yyy, g_x_z_xy_yyz, g_x_z_xy_yzz, g_x_z_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xy_xxx[k] = -g_x_z_x_xxx[k] * ab_y + g_x_z_x_xxxy[k];

                g_x_z_xy_xxy[k] = -g_x_z_x_xxy[k] * ab_y + g_x_z_x_xxyy[k];

                g_x_z_xy_xxz[k] = -g_x_z_x_xxz[k] * ab_y + g_x_z_x_xxyz[k];

                g_x_z_xy_xyy[k] = -g_x_z_x_xyy[k] * ab_y + g_x_z_x_xyyy[k];

                g_x_z_xy_xyz[k] = -g_x_z_x_xyz[k] * ab_y + g_x_z_x_xyyz[k];

                g_x_z_xy_xzz[k] = -g_x_z_x_xzz[k] * ab_y + g_x_z_x_xyzz[k];

                g_x_z_xy_yyy[k] = -g_x_z_x_yyy[k] * ab_y + g_x_z_x_yyyy[k];

                g_x_z_xy_yyz[k] = -g_x_z_x_yyz[k] * ab_y + g_x_z_x_yyyz[k];

                g_x_z_xy_yzz[k] = -g_x_z_x_yzz[k] * ab_y + g_x_z_x_yyzz[k];

                g_x_z_xy_zzz[k] = -g_x_z_x_zzz[k] * ab_y + g_x_z_x_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_x_z_xz_xxx = cbuffer.data(df_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_z_xz_xxy = cbuffer.data(df_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_z_xz_xxz = cbuffer.data(df_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_z_xz_xyy = cbuffer.data(df_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_z_xz_xyz = cbuffer.data(df_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_z_xz_xzz = cbuffer.data(df_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_z_xz_yyy = cbuffer.data(df_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_z_xz_yyz = cbuffer.data(df_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_z_xz_yzz = cbuffer.data(df_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_z_xz_zzz = cbuffer.data(df_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxx, g_0_z_z_xxy, g_0_z_z_xxz, g_0_z_z_xyy, g_0_z_z_xyz, g_0_z_z_xzz, g_0_z_z_yyy, g_0_z_z_yyz, g_0_z_z_yzz, g_0_z_z_zzz, g_x_z_xz_xxx, g_x_z_xz_xxy, g_x_z_xz_xxz, g_x_z_xz_xyy, g_x_z_xz_xyz, g_x_z_xz_xzz, g_x_z_xz_yyy, g_x_z_xz_yyz, g_x_z_xz_yzz, g_x_z_xz_zzz, g_x_z_z_xxx, g_x_z_z_xxxx, g_x_z_z_xxxy, g_x_z_z_xxxz, g_x_z_z_xxy, g_x_z_z_xxyy, g_x_z_z_xxyz, g_x_z_z_xxz, g_x_z_z_xxzz, g_x_z_z_xyy, g_x_z_z_xyyy, g_x_z_z_xyyz, g_x_z_z_xyz, g_x_z_z_xyzz, g_x_z_z_xzz, g_x_z_z_xzzz, g_x_z_z_yyy, g_x_z_z_yyz, g_x_z_z_yzz, g_x_z_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xz_xxx[k] = -g_0_z_z_xxx[k] - g_x_z_z_xxx[k] * ab_x + g_x_z_z_xxxx[k];

                g_x_z_xz_xxy[k] = -g_0_z_z_xxy[k] - g_x_z_z_xxy[k] * ab_x + g_x_z_z_xxxy[k];

                g_x_z_xz_xxz[k] = -g_0_z_z_xxz[k] - g_x_z_z_xxz[k] * ab_x + g_x_z_z_xxxz[k];

                g_x_z_xz_xyy[k] = -g_0_z_z_xyy[k] - g_x_z_z_xyy[k] * ab_x + g_x_z_z_xxyy[k];

                g_x_z_xz_xyz[k] = -g_0_z_z_xyz[k] - g_x_z_z_xyz[k] * ab_x + g_x_z_z_xxyz[k];

                g_x_z_xz_xzz[k] = -g_0_z_z_xzz[k] - g_x_z_z_xzz[k] * ab_x + g_x_z_z_xxzz[k];

                g_x_z_xz_yyy[k] = -g_0_z_z_yyy[k] - g_x_z_z_yyy[k] * ab_x + g_x_z_z_xyyy[k];

                g_x_z_xz_yyz[k] = -g_0_z_z_yyz[k] - g_x_z_z_yyz[k] * ab_x + g_x_z_z_xyyz[k];

                g_x_z_xz_yzz[k] = -g_0_z_z_yzz[k] - g_x_z_z_yzz[k] * ab_x + g_x_z_z_xyzz[k];

                g_x_z_xz_zzz[k] = -g_0_z_z_zzz[k] - g_x_z_z_zzz[k] * ab_x + g_x_z_z_xzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_x_z_yy_xxx = cbuffer.data(df_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_z_yy_xxy = cbuffer.data(df_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_z_yy_xxz = cbuffer.data(df_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_z_yy_xyy = cbuffer.data(df_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_z_yy_xyz = cbuffer.data(df_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_z_yy_xzz = cbuffer.data(df_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_z_yy_yyy = cbuffer.data(df_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_z_yy_yyz = cbuffer.data(df_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_z_yy_yzz = cbuffer.data(df_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_z_yy_zzz = cbuffer.data(df_geom_11_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_y_xxx, g_x_z_y_xxxy, g_x_z_y_xxy, g_x_z_y_xxyy, g_x_z_y_xxyz, g_x_z_y_xxz, g_x_z_y_xyy, g_x_z_y_xyyy, g_x_z_y_xyyz, g_x_z_y_xyz, g_x_z_y_xyzz, g_x_z_y_xzz, g_x_z_y_yyy, g_x_z_y_yyyy, g_x_z_y_yyyz, g_x_z_y_yyz, g_x_z_y_yyzz, g_x_z_y_yzz, g_x_z_y_yzzz, g_x_z_y_zzz, g_x_z_yy_xxx, g_x_z_yy_xxy, g_x_z_yy_xxz, g_x_z_yy_xyy, g_x_z_yy_xyz, g_x_z_yy_xzz, g_x_z_yy_yyy, g_x_z_yy_yyz, g_x_z_yy_yzz, g_x_z_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yy_xxx[k] = -g_x_z_y_xxx[k] * ab_y + g_x_z_y_xxxy[k];

                g_x_z_yy_xxy[k] = -g_x_z_y_xxy[k] * ab_y + g_x_z_y_xxyy[k];

                g_x_z_yy_xxz[k] = -g_x_z_y_xxz[k] * ab_y + g_x_z_y_xxyz[k];

                g_x_z_yy_xyy[k] = -g_x_z_y_xyy[k] * ab_y + g_x_z_y_xyyy[k];

                g_x_z_yy_xyz[k] = -g_x_z_y_xyz[k] * ab_y + g_x_z_y_xyyz[k];

                g_x_z_yy_xzz[k] = -g_x_z_y_xzz[k] * ab_y + g_x_z_y_xyzz[k];

                g_x_z_yy_yyy[k] = -g_x_z_y_yyy[k] * ab_y + g_x_z_y_yyyy[k];

                g_x_z_yy_yyz[k] = -g_x_z_y_yyz[k] * ab_y + g_x_z_y_yyyz[k];

                g_x_z_yy_yzz[k] = -g_x_z_y_yzz[k] * ab_y + g_x_z_y_yyzz[k];

                g_x_z_yy_zzz[k] = -g_x_z_y_zzz[k] * ab_y + g_x_z_y_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_x_z_yz_xxx = cbuffer.data(df_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_z_yz_xxy = cbuffer.data(df_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_z_yz_xxz = cbuffer.data(df_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_z_yz_xyy = cbuffer.data(df_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_z_yz_xyz = cbuffer.data(df_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_z_yz_xzz = cbuffer.data(df_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_z_yz_yyy = cbuffer.data(df_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_z_yz_yyz = cbuffer.data(df_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_z_yz_yzz = cbuffer.data(df_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_z_yz_zzz = cbuffer.data(df_geom_11_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yz_xxx, g_x_z_yz_xxy, g_x_z_yz_xxz, g_x_z_yz_xyy, g_x_z_yz_xyz, g_x_z_yz_xzz, g_x_z_yz_yyy, g_x_z_yz_yyz, g_x_z_yz_yzz, g_x_z_yz_zzz, g_x_z_z_xxx, g_x_z_z_xxxy, g_x_z_z_xxy, g_x_z_z_xxyy, g_x_z_z_xxyz, g_x_z_z_xxz, g_x_z_z_xyy, g_x_z_z_xyyy, g_x_z_z_xyyz, g_x_z_z_xyz, g_x_z_z_xyzz, g_x_z_z_xzz, g_x_z_z_yyy, g_x_z_z_yyyy, g_x_z_z_yyyz, g_x_z_z_yyz, g_x_z_z_yyzz, g_x_z_z_yzz, g_x_z_z_yzzz, g_x_z_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yz_xxx[k] = -g_x_z_z_xxx[k] * ab_y + g_x_z_z_xxxy[k];

                g_x_z_yz_xxy[k] = -g_x_z_z_xxy[k] * ab_y + g_x_z_z_xxyy[k];

                g_x_z_yz_xxz[k] = -g_x_z_z_xxz[k] * ab_y + g_x_z_z_xxyz[k];

                g_x_z_yz_xyy[k] = -g_x_z_z_xyy[k] * ab_y + g_x_z_z_xyyy[k];

                g_x_z_yz_xyz[k] = -g_x_z_z_xyz[k] * ab_y + g_x_z_z_xyyz[k];

                g_x_z_yz_xzz[k] = -g_x_z_z_xzz[k] * ab_y + g_x_z_z_xyzz[k];

                g_x_z_yz_yyy[k] = -g_x_z_z_yyy[k] * ab_y + g_x_z_z_yyyy[k];

                g_x_z_yz_yyz[k] = -g_x_z_z_yyz[k] * ab_y + g_x_z_z_yyyz[k];

                g_x_z_yz_yzz[k] = -g_x_z_z_yzz[k] * ab_y + g_x_z_z_yyzz[k];

                g_x_z_yz_zzz[k] = -g_x_z_z_zzz[k] * ab_y + g_x_z_z_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_x_z_zz_xxx = cbuffer.data(df_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_z_zz_xxy = cbuffer.data(df_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_z_zz_xxz = cbuffer.data(df_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_z_zz_xyy = cbuffer.data(df_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_z_zz_xyz = cbuffer.data(df_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_z_zz_xzz = cbuffer.data(df_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_z_zz_yyy = cbuffer.data(df_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_z_zz_yyz = cbuffer.data(df_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_z_zz_yzz = cbuffer.data(df_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_z_zz_zzz = cbuffer.data(df_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxx, g_x_0_z_xxy, g_x_0_z_xxz, g_x_0_z_xyy, g_x_0_z_xyz, g_x_0_z_xzz, g_x_0_z_yyy, g_x_0_z_yyz, g_x_0_z_yzz, g_x_0_z_zzz, g_x_z_z_xxx, g_x_z_z_xxxz, g_x_z_z_xxy, g_x_z_z_xxyz, g_x_z_z_xxz, g_x_z_z_xxzz, g_x_z_z_xyy, g_x_z_z_xyyz, g_x_z_z_xyz, g_x_z_z_xyzz, g_x_z_z_xzz, g_x_z_z_xzzz, g_x_z_z_yyy, g_x_z_z_yyyz, g_x_z_z_yyz, g_x_z_z_yyzz, g_x_z_z_yzz, g_x_z_z_yzzz, g_x_z_z_zzz, g_x_z_z_zzzz, g_x_z_zz_xxx, g_x_z_zz_xxy, g_x_z_zz_xxz, g_x_z_zz_xyy, g_x_z_zz_xyz, g_x_z_zz_xzz, g_x_z_zz_yyy, g_x_z_zz_yyz, g_x_z_zz_yzz, g_x_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zz_xxx[k] = g_x_0_z_xxx[k] - g_x_z_z_xxx[k] * ab_z + g_x_z_z_xxxz[k];

                g_x_z_zz_xxy[k] = g_x_0_z_xxy[k] - g_x_z_z_xxy[k] * ab_z + g_x_z_z_xxyz[k];

                g_x_z_zz_xxz[k] = g_x_0_z_xxz[k] - g_x_z_z_xxz[k] * ab_z + g_x_z_z_xxzz[k];

                g_x_z_zz_xyy[k] = g_x_0_z_xyy[k] - g_x_z_z_xyy[k] * ab_z + g_x_z_z_xyyz[k];

                g_x_z_zz_xyz[k] = g_x_0_z_xyz[k] - g_x_z_z_xyz[k] * ab_z + g_x_z_z_xyzz[k];

                g_x_z_zz_xzz[k] = g_x_0_z_xzz[k] - g_x_z_z_xzz[k] * ab_z + g_x_z_z_xzzz[k];

                g_x_z_zz_yyy[k] = g_x_0_z_yyy[k] - g_x_z_z_yyy[k] * ab_z + g_x_z_z_yyyz[k];

                g_x_z_zz_yyz[k] = g_x_0_z_yyz[k] - g_x_z_z_yyz[k] * ab_z + g_x_z_z_yyzz[k];

                g_x_z_zz_yzz[k] = g_x_0_z_yzz[k] - g_x_z_z_yzz[k] * ab_z + g_x_z_z_yzzz[k];

                g_x_z_zz_zzz[k] = g_x_0_z_zzz[k] - g_x_z_z_zzz[k] * ab_z + g_x_z_z_zzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_y_x_xx_xxx = cbuffer.data(df_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_x_xx_xxy = cbuffer.data(df_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_x_xx_xxz = cbuffer.data(df_geom_11_off + 182 * ccomps * dcomps);

            auto g_y_x_xx_xyy = cbuffer.data(df_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_x_xx_xyz = cbuffer.data(df_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_x_xx_xzz = cbuffer.data(df_geom_11_off + 185 * ccomps * dcomps);

            auto g_y_x_xx_yyy = cbuffer.data(df_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_x_xx_yyz = cbuffer.data(df_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_x_xx_yzz = cbuffer.data(df_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_x_xx_zzz = cbuffer.data(df_geom_11_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_xxx, g_y_0_x_xxy, g_y_0_x_xxz, g_y_0_x_xyy, g_y_0_x_xyz, g_y_0_x_xzz, g_y_0_x_yyy, g_y_0_x_yyz, g_y_0_x_yzz, g_y_0_x_zzz, g_y_x_x_xxx, g_y_x_x_xxxx, g_y_x_x_xxxy, g_y_x_x_xxxz, g_y_x_x_xxy, g_y_x_x_xxyy, g_y_x_x_xxyz, g_y_x_x_xxz, g_y_x_x_xxzz, g_y_x_x_xyy, g_y_x_x_xyyy, g_y_x_x_xyyz, g_y_x_x_xyz, g_y_x_x_xyzz, g_y_x_x_xzz, g_y_x_x_xzzz, g_y_x_x_yyy, g_y_x_x_yyz, g_y_x_x_yzz, g_y_x_x_zzz, g_y_x_xx_xxx, g_y_x_xx_xxy, g_y_x_xx_xxz, g_y_x_xx_xyy, g_y_x_xx_xyz, g_y_x_xx_xzz, g_y_x_xx_yyy, g_y_x_xx_yyz, g_y_x_xx_yzz, g_y_x_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xx_xxx[k] = g_y_0_x_xxx[k] - g_y_x_x_xxx[k] * ab_x + g_y_x_x_xxxx[k];

                g_y_x_xx_xxy[k] = g_y_0_x_xxy[k] - g_y_x_x_xxy[k] * ab_x + g_y_x_x_xxxy[k];

                g_y_x_xx_xxz[k] = g_y_0_x_xxz[k] - g_y_x_x_xxz[k] * ab_x + g_y_x_x_xxxz[k];

                g_y_x_xx_xyy[k] = g_y_0_x_xyy[k] - g_y_x_x_xyy[k] * ab_x + g_y_x_x_xxyy[k];

                g_y_x_xx_xyz[k] = g_y_0_x_xyz[k] - g_y_x_x_xyz[k] * ab_x + g_y_x_x_xxyz[k];

                g_y_x_xx_xzz[k] = g_y_0_x_xzz[k] - g_y_x_x_xzz[k] * ab_x + g_y_x_x_xxzz[k];

                g_y_x_xx_yyy[k] = g_y_0_x_yyy[k] - g_y_x_x_yyy[k] * ab_x + g_y_x_x_xyyy[k];

                g_y_x_xx_yyz[k] = g_y_0_x_yyz[k] - g_y_x_x_yyz[k] * ab_x + g_y_x_x_xyyz[k];

                g_y_x_xx_yzz[k] = g_y_0_x_yzz[k] - g_y_x_x_yzz[k] * ab_x + g_y_x_x_xyzz[k];

                g_y_x_xx_zzz[k] = g_y_0_x_zzz[k] - g_y_x_x_zzz[k] * ab_x + g_y_x_x_xzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_y_x_xy_xxx = cbuffer.data(df_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_x_xy_xxy = cbuffer.data(df_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_x_xy_xxz = cbuffer.data(df_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_x_xy_xyy = cbuffer.data(df_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_x_xy_xyz = cbuffer.data(df_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_x_xy_xzz = cbuffer.data(df_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_x_xy_yyy = cbuffer.data(df_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_x_xy_yyz = cbuffer.data(df_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_x_xy_yzz = cbuffer.data(df_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_x_xy_zzz = cbuffer.data(df_geom_11_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxx, g_y_0_y_xxy, g_y_0_y_xxz, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xzz, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yzz, g_y_0_y_zzz, g_y_x_xy_xxx, g_y_x_xy_xxy, g_y_x_xy_xxz, g_y_x_xy_xyy, g_y_x_xy_xyz, g_y_x_xy_xzz, g_y_x_xy_yyy, g_y_x_xy_yyz, g_y_x_xy_yzz, g_y_x_xy_zzz, g_y_x_y_xxx, g_y_x_y_xxxx, g_y_x_y_xxxy, g_y_x_y_xxxz, g_y_x_y_xxy, g_y_x_y_xxyy, g_y_x_y_xxyz, g_y_x_y_xxz, g_y_x_y_xxzz, g_y_x_y_xyy, g_y_x_y_xyyy, g_y_x_y_xyyz, g_y_x_y_xyz, g_y_x_y_xyzz, g_y_x_y_xzz, g_y_x_y_xzzz, g_y_x_y_yyy, g_y_x_y_yyz, g_y_x_y_yzz, g_y_x_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xy_xxx[k] = g_y_0_y_xxx[k] - g_y_x_y_xxx[k] * ab_x + g_y_x_y_xxxx[k];

                g_y_x_xy_xxy[k] = g_y_0_y_xxy[k] - g_y_x_y_xxy[k] * ab_x + g_y_x_y_xxxy[k];

                g_y_x_xy_xxz[k] = g_y_0_y_xxz[k] - g_y_x_y_xxz[k] * ab_x + g_y_x_y_xxxz[k];

                g_y_x_xy_xyy[k] = g_y_0_y_xyy[k] - g_y_x_y_xyy[k] * ab_x + g_y_x_y_xxyy[k];

                g_y_x_xy_xyz[k] = g_y_0_y_xyz[k] - g_y_x_y_xyz[k] * ab_x + g_y_x_y_xxyz[k];

                g_y_x_xy_xzz[k] = g_y_0_y_xzz[k] - g_y_x_y_xzz[k] * ab_x + g_y_x_y_xxzz[k];

                g_y_x_xy_yyy[k] = g_y_0_y_yyy[k] - g_y_x_y_yyy[k] * ab_x + g_y_x_y_xyyy[k];

                g_y_x_xy_yyz[k] = g_y_0_y_yyz[k] - g_y_x_y_yyz[k] * ab_x + g_y_x_y_xyyz[k];

                g_y_x_xy_yzz[k] = g_y_0_y_yzz[k] - g_y_x_y_yzz[k] * ab_x + g_y_x_y_xyzz[k];

                g_y_x_xy_zzz[k] = g_y_0_y_zzz[k] - g_y_x_y_zzz[k] * ab_x + g_y_x_y_xzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_y_x_xz_xxx = cbuffer.data(df_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_x_xz_xxy = cbuffer.data(df_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_x_xz_xxz = cbuffer.data(df_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_x_xz_xyy = cbuffer.data(df_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_x_xz_xyz = cbuffer.data(df_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_x_xz_xzz = cbuffer.data(df_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_x_xz_yyy = cbuffer.data(df_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_x_xz_yyz = cbuffer.data(df_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_x_xz_yzz = cbuffer.data(df_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_x_xz_zzz = cbuffer.data(df_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_x_xxx, g_y_x_x_xxxz, g_y_x_x_xxy, g_y_x_x_xxyz, g_y_x_x_xxz, g_y_x_x_xxzz, g_y_x_x_xyy, g_y_x_x_xyyz, g_y_x_x_xyz, g_y_x_x_xyzz, g_y_x_x_xzz, g_y_x_x_xzzz, g_y_x_x_yyy, g_y_x_x_yyyz, g_y_x_x_yyz, g_y_x_x_yyzz, g_y_x_x_yzz, g_y_x_x_yzzz, g_y_x_x_zzz, g_y_x_x_zzzz, g_y_x_xz_xxx, g_y_x_xz_xxy, g_y_x_xz_xxz, g_y_x_xz_xyy, g_y_x_xz_xyz, g_y_x_xz_xzz, g_y_x_xz_yyy, g_y_x_xz_yyz, g_y_x_xz_yzz, g_y_x_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xz_xxx[k] = -g_y_x_x_xxx[k] * ab_z + g_y_x_x_xxxz[k];

                g_y_x_xz_xxy[k] = -g_y_x_x_xxy[k] * ab_z + g_y_x_x_xxyz[k];

                g_y_x_xz_xxz[k] = -g_y_x_x_xxz[k] * ab_z + g_y_x_x_xxzz[k];

                g_y_x_xz_xyy[k] = -g_y_x_x_xyy[k] * ab_z + g_y_x_x_xyyz[k];

                g_y_x_xz_xyz[k] = -g_y_x_x_xyz[k] * ab_z + g_y_x_x_xyzz[k];

                g_y_x_xz_xzz[k] = -g_y_x_x_xzz[k] * ab_z + g_y_x_x_xzzz[k];

                g_y_x_xz_yyy[k] = -g_y_x_x_yyy[k] * ab_z + g_y_x_x_yyyz[k];

                g_y_x_xz_yyz[k] = -g_y_x_x_yyz[k] * ab_z + g_y_x_x_yyzz[k];

                g_y_x_xz_yzz[k] = -g_y_x_x_yzz[k] * ab_z + g_y_x_x_yzzz[k];

                g_y_x_xz_zzz[k] = -g_y_x_x_zzz[k] * ab_z + g_y_x_x_zzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_y_x_yy_xxx = cbuffer.data(df_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_x_yy_xxy = cbuffer.data(df_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_x_yy_xxz = cbuffer.data(df_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_x_yy_xyy = cbuffer.data(df_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_x_yy_xyz = cbuffer.data(df_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_x_yy_xzz = cbuffer.data(df_geom_11_off + 215 * ccomps * dcomps);

            auto g_y_x_yy_yyy = cbuffer.data(df_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_x_yy_yyz = cbuffer.data(df_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_x_yy_yzz = cbuffer.data(df_geom_11_off + 218 * ccomps * dcomps);

            auto g_y_x_yy_zzz = cbuffer.data(df_geom_11_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_y_xxx, g_0_x_y_xxy, g_0_x_y_xxz, g_0_x_y_xyy, g_0_x_y_xyz, g_0_x_y_xzz, g_0_x_y_yyy, g_0_x_y_yyz, g_0_x_y_yzz, g_0_x_y_zzz, g_y_x_y_xxx, g_y_x_y_xxxy, g_y_x_y_xxy, g_y_x_y_xxyy, g_y_x_y_xxyz, g_y_x_y_xxz, g_y_x_y_xyy, g_y_x_y_xyyy, g_y_x_y_xyyz, g_y_x_y_xyz, g_y_x_y_xyzz, g_y_x_y_xzz, g_y_x_y_yyy, g_y_x_y_yyyy, g_y_x_y_yyyz, g_y_x_y_yyz, g_y_x_y_yyzz, g_y_x_y_yzz, g_y_x_y_yzzz, g_y_x_y_zzz, g_y_x_yy_xxx, g_y_x_yy_xxy, g_y_x_yy_xxz, g_y_x_yy_xyy, g_y_x_yy_xyz, g_y_x_yy_xzz, g_y_x_yy_yyy, g_y_x_yy_yyz, g_y_x_yy_yzz, g_y_x_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yy_xxx[k] = -g_0_x_y_xxx[k] - g_y_x_y_xxx[k] * ab_y + g_y_x_y_xxxy[k];

                g_y_x_yy_xxy[k] = -g_0_x_y_xxy[k] - g_y_x_y_xxy[k] * ab_y + g_y_x_y_xxyy[k];

                g_y_x_yy_xxz[k] = -g_0_x_y_xxz[k] - g_y_x_y_xxz[k] * ab_y + g_y_x_y_xxyz[k];

                g_y_x_yy_xyy[k] = -g_0_x_y_xyy[k] - g_y_x_y_xyy[k] * ab_y + g_y_x_y_xyyy[k];

                g_y_x_yy_xyz[k] = -g_0_x_y_xyz[k] - g_y_x_y_xyz[k] * ab_y + g_y_x_y_xyyz[k];

                g_y_x_yy_xzz[k] = -g_0_x_y_xzz[k] - g_y_x_y_xzz[k] * ab_y + g_y_x_y_xyzz[k];

                g_y_x_yy_yyy[k] = -g_0_x_y_yyy[k] - g_y_x_y_yyy[k] * ab_y + g_y_x_y_yyyy[k];

                g_y_x_yy_yyz[k] = -g_0_x_y_yyz[k] - g_y_x_y_yyz[k] * ab_y + g_y_x_y_yyyz[k];

                g_y_x_yy_yzz[k] = -g_0_x_y_yzz[k] - g_y_x_y_yzz[k] * ab_y + g_y_x_y_yyzz[k];

                g_y_x_yy_zzz[k] = -g_0_x_y_zzz[k] - g_y_x_y_zzz[k] * ab_y + g_y_x_y_yzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_y_x_yz_xxx = cbuffer.data(df_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_x_yz_xxy = cbuffer.data(df_geom_11_off + 221 * ccomps * dcomps);

            auto g_y_x_yz_xxz = cbuffer.data(df_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_x_yz_xyy = cbuffer.data(df_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_x_yz_xyz = cbuffer.data(df_geom_11_off + 224 * ccomps * dcomps);

            auto g_y_x_yz_xzz = cbuffer.data(df_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_x_yz_yyy = cbuffer.data(df_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_x_yz_yyz = cbuffer.data(df_geom_11_off + 227 * ccomps * dcomps);

            auto g_y_x_yz_yzz = cbuffer.data(df_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_x_yz_zzz = cbuffer.data(df_geom_11_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_y_xxx, g_y_x_y_xxxz, g_y_x_y_xxy, g_y_x_y_xxyz, g_y_x_y_xxz, g_y_x_y_xxzz, g_y_x_y_xyy, g_y_x_y_xyyz, g_y_x_y_xyz, g_y_x_y_xyzz, g_y_x_y_xzz, g_y_x_y_xzzz, g_y_x_y_yyy, g_y_x_y_yyyz, g_y_x_y_yyz, g_y_x_y_yyzz, g_y_x_y_yzz, g_y_x_y_yzzz, g_y_x_y_zzz, g_y_x_y_zzzz, g_y_x_yz_xxx, g_y_x_yz_xxy, g_y_x_yz_xxz, g_y_x_yz_xyy, g_y_x_yz_xyz, g_y_x_yz_xzz, g_y_x_yz_yyy, g_y_x_yz_yyz, g_y_x_yz_yzz, g_y_x_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yz_xxx[k] = -g_y_x_y_xxx[k] * ab_z + g_y_x_y_xxxz[k];

                g_y_x_yz_xxy[k] = -g_y_x_y_xxy[k] * ab_z + g_y_x_y_xxyz[k];

                g_y_x_yz_xxz[k] = -g_y_x_y_xxz[k] * ab_z + g_y_x_y_xxzz[k];

                g_y_x_yz_xyy[k] = -g_y_x_y_xyy[k] * ab_z + g_y_x_y_xyyz[k];

                g_y_x_yz_xyz[k] = -g_y_x_y_xyz[k] * ab_z + g_y_x_y_xyzz[k];

                g_y_x_yz_xzz[k] = -g_y_x_y_xzz[k] * ab_z + g_y_x_y_xzzz[k];

                g_y_x_yz_yyy[k] = -g_y_x_y_yyy[k] * ab_z + g_y_x_y_yyyz[k];

                g_y_x_yz_yyz[k] = -g_y_x_y_yyz[k] * ab_z + g_y_x_y_yyzz[k];

                g_y_x_yz_yzz[k] = -g_y_x_y_yzz[k] * ab_z + g_y_x_y_yzzz[k];

                g_y_x_yz_zzz[k] = -g_y_x_y_zzz[k] * ab_z + g_y_x_y_zzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_y_x_zz_xxx = cbuffer.data(df_geom_11_off + 230 * ccomps * dcomps);

            auto g_y_x_zz_xxy = cbuffer.data(df_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_x_zz_xxz = cbuffer.data(df_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_x_zz_xyy = cbuffer.data(df_geom_11_off + 233 * ccomps * dcomps);

            auto g_y_x_zz_xyz = cbuffer.data(df_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_x_zz_xzz = cbuffer.data(df_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_x_zz_yyy = cbuffer.data(df_geom_11_off + 236 * ccomps * dcomps);

            auto g_y_x_zz_yyz = cbuffer.data(df_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_x_zz_yzz = cbuffer.data(df_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_x_zz_zzz = cbuffer.data(df_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_z_xxx, g_y_x_z_xxxz, g_y_x_z_xxy, g_y_x_z_xxyz, g_y_x_z_xxz, g_y_x_z_xxzz, g_y_x_z_xyy, g_y_x_z_xyyz, g_y_x_z_xyz, g_y_x_z_xyzz, g_y_x_z_xzz, g_y_x_z_xzzz, g_y_x_z_yyy, g_y_x_z_yyyz, g_y_x_z_yyz, g_y_x_z_yyzz, g_y_x_z_yzz, g_y_x_z_yzzz, g_y_x_z_zzz, g_y_x_z_zzzz, g_y_x_zz_xxx, g_y_x_zz_xxy, g_y_x_zz_xxz, g_y_x_zz_xyy, g_y_x_zz_xyz, g_y_x_zz_xzz, g_y_x_zz_yyy, g_y_x_zz_yyz, g_y_x_zz_yzz, g_y_x_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zz_xxx[k] = -g_y_x_z_xxx[k] * ab_z + g_y_x_z_xxxz[k];

                g_y_x_zz_xxy[k] = -g_y_x_z_xxy[k] * ab_z + g_y_x_z_xxyz[k];

                g_y_x_zz_xxz[k] = -g_y_x_z_xxz[k] * ab_z + g_y_x_z_xxzz[k];

                g_y_x_zz_xyy[k] = -g_y_x_z_xyy[k] * ab_z + g_y_x_z_xyyz[k];

                g_y_x_zz_xyz[k] = -g_y_x_z_xyz[k] * ab_z + g_y_x_z_xyzz[k];

                g_y_x_zz_xzz[k] = -g_y_x_z_xzz[k] * ab_z + g_y_x_z_xzzz[k];

                g_y_x_zz_yyy[k] = -g_y_x_z_yyy[k] * ab_z + g_y_x_z_yyyz[k];

                g_y_x_zz_yyz[k] = -g_y_x_z_yyz[k] * ab_z + g_y_x_z_yyzz[k];

                g_y_x_zz_yzz[k] = -g_y_x_z_yzz[k] * ab_z + g_y_x_z_yzzz[k];

                g_y_x_zz_zzz[k] = -g_y_x_z_zzz[k] * ab_z + g_y_x_z_zzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_y_y_xx_xxx = cbuffer.data(df_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_y_xx_xxy = cbuffer.data(df_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_y_xx_xxz = cbuffer.data(df_geom_11_off + 242 * ccomps * dcomps);

            auto g_y_y_xx_xyy = cbuffer.data(df_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_y_xx_xyz = cbuffer.data(df_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_y_xx_xzz = cbuffer.data(df_geom_11_off + 245 * ccomps * dcomps);

            auto g_y_y_xx_yyy = cbuffer.data(df_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_y_xx_yyz = cbuffer.data(df_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_y_xx_yzz = cbuffer.data(df_geom_11_off + 248 * ccomps * dcomps);

            auto g_y_y_xx_zzz = cbuffer.data(df_geom_11_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_x_xxx, g_y_y_x_xxxx, g_y_y_x_xxxy, g_y_y_x_xxxz, g_y_y_x_xxy, g_y_y_x_xxyy, g_y_y_x_xxyz, g_y_y_x_xxz, g_y_y_x_xxzz, g_y_y_x_xyy, g_y_y_x_xyyy, g_y_y_x_xyyz, g_y_y_x_xyz, g_y_y_x_xyzz, g_y_y_x_xzz, g_y_y_x_xzzz, g_y_y_x_yyy, g_y_y_x_yyz, g_y_y_x_yzz, g_y_y_x_zzz, g_y_y_xx_xxx, g_y_y_xx_xxy, g_y_y_xx_xxz, g_y_y_xx_xyy, g_y_y_xx_xyz, g_y_y_xx_xzz, g_y_y_xx_yyy, g_y_y_xx_yyz, g_y_y_xx_yzz, g_y_y_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xx_xxx[k] = -g_y_y_x_xxx[k] * ab_x + g_y_y_x_xxxx[k];

                g_y_y_xx_xxy[k] = -g_y_y_x_xxy[k] * ab_x + g_y_y_x_xxxy[k];

                g_y_y_xx_xxz[k] = -g_y_y_x_xxz[k] * ab_x + g_y_y_x_xxxz[k];

                g_y_y_xx_xyy[k] = -g_y_y_x_xyy[k] * ab_x + g_y_y_x_xxyy[k];

                g_y_y_xx_xyz[k] = -g_y_y_x_xyz[k] * ab_x + g_y_y_x_xxyz[k];

                g_y_y_xx_xzz[k] = -g_y_y_x_xzz[k] * ab_x + g_y_y_x_xxzz[k];

                g_y_y_xx_yyy[k] = -g_y_y_x_yyy[k] * ab_x + g_y_y_x_xyyy[k];

                g_y_y_xx_yyz[k] = -g_y_y_x_yyz[k] * ab_x + g_y_y_x_xyyz[k];

                g_y_y_xx_yzz[k] = -g_y_y_x_yzz[k] * ab_x + g_y_y_x_xyzz[k];

                g_y_y_xx_zzz[k] = -g_y_y_x_zzz[k] * ab_x + g_y_y_x_xzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_y_y_xy_xxx = cbuffer.data(df_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_y_xy_xxy = cbuffer.data(df_geom_11_off + 251 * ccomps * dcomps);

            auto g_y_y_xy_xxz = cbuffer.data(df_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_y_xy_xyy = cbuffer.data(df_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_y_xy_xyz = cbuffer.data(df_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_y_xy_xzz = cbuffer.data(df_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_y_xy_yyy = cbuffer.data(df_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_y_xy_yyz = cbuffer.data(df_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_y_xy_yzz = cbuffer.data(df_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_y_xy_zzz = cbuffer.data(df_geom_11_off + 259 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xy_xxx, g_y_y_xy_xxy, g_y_y_xy_xxz, g_y_y_xy_xyy, g_y_y_xy_xyz, g_y_y_xy_xzz, g_y_y_xy_yyy, g_y_y_xy_yyz, g_y_y_xy_yzz, g_y_y_xy_zzz, g_y_y_y_xxx, g_y_y_y_xxxx, g_y_y_y_xxxy, g_y_y_y_xxxz, g_y_y_y_xxy, g_y_y_y_xxyy, g_y_y_y_xxyz, g_y_y_y_xxz, g_y_y_y_xxzz, g_y_y_y_xyy, g_y_y_y_xyyy, g_y_y_y_xyyz, g_y_y_y_xyz, g_y_y_y_xyzz, g_y_y_y_xzz, g_y_y_y_xzzz, g_y_y_y_yyy, g_y_y_y_yyz, g_y_y_y_yzz, g_y_y_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xy_xxx[k] = -g_y_y_y_xxx[k] * ab_x + g_y_y_y_xxxx[k];

                g_y_y_xy_xxy[k] = -g_y_y_y_xxy[k] * ab_x + g_y_y_y_xxxy[k];

                g_y_y_xy_xxz[k] = -g_y_y_y_xxz[k] * ab_x + g_y_y_y_xxxz[k];

                g_y_y_xy_xyy[k] = -g_y_y_y_xyy[k] * ab_x + g_y_y_y_xxyy[k];

                g_y_y_xy_xyz[k] = -g_y_y_y_xyz[k] * ab_x + g_y_y_y_xxyz[k];

                g_y_y_xy_xzz[k] = -g_y_y_y_xzz[k] * ab_x + g_y_y_y_xxzz[k];

                g_y_y_xy_yyy[k] = -g_y_y_y_yyy[k] * ab_x + g_y_y_y_xyyy[k];

                g_y_y_xy_yyz[k] = -g_y_y_y_yyz[k] * ab_x + g_y_y_y_xyyz[k];

                g_y_y_xy_yzz[k] = -g_y_y_y_yzz[k] * ab_x + g_y_y_y_xyzz[k];

                g_y_y_xy_zzz[k] = -g_y_y_y_zzz[k] * ab_x + g_y_y_y_xzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_y_y_xz_xxx = cbuffer.data(df_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_y_xz_xxy = cbuffer.data(df_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_y_xz_xxz = cbuffer.data(df_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_y_xz_xyy = cbuffer.data(df_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_y_xz_xyz = cbuffer.data(df_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_y_xz_xzz = cbuffer.data(df_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_y_xz_yyy = cbuffer.data(df_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_y_xz_yyz = cbuffer.data(df_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_y_xz_yzz = cbuffer.data(df_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_y_xz_zzz = cbuffer.data(df_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xz_xxx, g_y_y_xz_xxy, g_y_y_xz_xxz, g_y_y_xz_xyy, g_y_y_xz_xyz, g_y_y_xz_xzz, g_y_y_xz_yyy, g_y_y_xz_yyz, g_y_y_xz_yzz, g_y_y_xz_zzz, g_y_y_z_xxx, g_y_y_z_xxxx, g_y_y_z_xxxy, g_y_y_z_xxxz, g_y_y_z_xxy, g_y_y_z_xxyy, g_y_y_z_xxyz, g_y_y_z_xxz, g_y_y_z_xxzz, g_y_y_z_xyy, g_y_y_z_xyyy, g_y_y_z_xyyz, g_y_y_z_xyz, g_y_y_z_xyzz, g_y_y_z_xzz, g_y_y_z_xzzz, g_y_y_z_yyy, g_y_y_z_yyz, g_y_y_z_yzz, g_y_y_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xz_xxx[k] = -g_y_y_z_xxx[k] * ab_x + g_y_y_z_xxxx[k];

                g_y_y_xz_xxy[k] = -g_y_y_z_xxy[k] * ab_x + g_y_y_z_xxxy[k];

                g_y_y_xz_xxz[k] = -g_y_y_z_xxz[k] * ab_x + g_y_y_z_xxxz[k];

                g_y_y_xz_xyy[k] = -g_y_y_z_xyy[k] * ab_x + g_y_y_z_xxyy[k];

                g_y_y_xz_xyz[k] = -g_y_y_z_xyz[k] * ab_x + g_y_y_z_xxyz[k];

                g_y_y_xz_xzz[k] = -g_y_y_z_xzz[k] * ab_x + g_y_y_z_xxzz[k];

                g_y_y_xz_yyy[k] = -g_y_y_z_yyy[k] * ab_x + g_y_y_z_xyyy[k];

                g_y_y_xz_yyz[k] = -g_y_y_z_yyz[k] * ab_x + g_y_y_z_xyyz[k];

                g_y_y_xz_yzz[k] = -g_y_y_z_yzz[k] * ab_x + g_y_y_z_xyzz[k];

                g_y_y_xz_zzz[k] = -g_y_y_z_zzz[k] * ab_x + g_y_y_z_xzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_y_y_yy_xxx = cbuffer.data(df_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_y_yy_xxy = cbuffer.data(df_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_y_yy_xxz = cbuffer.data(df_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_y_yy_xyy = cbuffer.data(df_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_y_yy_xyz = cbuffer.data(df_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_y_yy_xzz = cbuffer.data(df_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_y_yy_yyy = cbuffer.data(df_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_y_yy_yyz = cbuffer.data(df_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_y_yy_yzz = cbuffer.data(df_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_y_yy_zzz = cbuffer.data(df_geom_11_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxx, g_0_y_y_xxy, g_0_y_y_xxz, g_0_y_y_xyy, g_0_y_y_xyz, g_0_y_y_xzz, g_0_y_y_yyy, g_0_y_y_yyz, g_0_y_y_yzz, g_0_y_y_zzz, g_y_0_y_xxx, g_y_0_y_xxy, g_y_0_y_xxz, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xzz, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yzz, g_y_0_y_zzz, g_y_y_y_xxx, g_y_y_y_xxxy, g_y_y_y_xxy, g_y_y_y_xxyy, g_y_y_y_xxyz, g_y_y_y_xxz, g_y_y_y_xyy, g_y_y_y_xyyy, g_y_y_y_xyyz, g_y_y_y_xyz, g_y_y_y_xyzz, g_y_y_y_xzz, g_y_y_y_yyy, g_y_y_y_yyyy, g_y_y_y_yyyz, g_y_y_y_yyz, g_y_y_y_yyzz, g_y_y_y_yzz, g_y_y_y_yzzz, g_y_y_y_zzz, g_y_y_yy_xxx, g_y_y_yy_xxy, g_y_y_yy_xxz, g_y_y_yy_xyy, g_y_y_yy_xyz, g_y_y_yy_xzz, g_y_y_yy_yyy, g_y_y_yy_yyz, g_y_y_yy_yzz, g_y_y_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yy_xxx[k] = -g_0_y_y_xxx[k] + g_y_0_y_xxx[k] - g_y_y_y_xxx[k] * ab_y + g_y_y_y_xxxy[k];

                g_y_y_yy_xxy[k] = -g_0_y_y_xxy[k] + g_y_0_y_xxy[k] - g_y_y_y_xxy[k] * ab_y + g_y_y_y_xxyy[k];

                g_y_y_yy_xxz[k] = -g_0_y_y_xxz[k] + g_y_0_y_xxz[k] - g_y_y_y_xxz[k] * ab_y + g_y_y_y_xxyz[k];

                g_y_y_yy_xyy[k] = -g_0_y_y_xyy[k] + g_y_0_y_xyy[k] - g_y_y_y_xyy[k] * ab_y + g_y_y_y_xyyy[k];

                g_y_y_yy_xyz[k] = -g_0_y_y_xyz[k] + g_y_0_y_xyz[k] - g_y_y_y_xyz[k] * ab_y + g_y_y_y_xyyz[k];

                g_y_y_yy_xzz[k] = -g_0_y_y_xzz[k] + g_y_0_y_xzz[k] - g_y_y_y_xzz[k] * ab_y + g_y_y_y_xyzz[k];

                g_y_y_yy_yyy[k] = -g_0_y_y_yyy[k] + g_y_0_y_yyy[k] - g_y_y_y_yyy[k] * ab_y + g_y_y_y_yyyy[k];

                g_y_y_yy_yyz[k] = -g_0_y_y_yyz[k] + g_y_0_y_yyz[k] - g_y_y_y_yyz[k] * ab_y + g_y_y_y_yyyz[k];

                g_y_y_yy_yzz[k] = -g_0_y_y_yzz[k] + g_y_0_y_yzz[k] - g_y_y_y_yzz[k] * ab_y + g_y_y_y_yyzz[k];

                g_y_y_yy_zzz[k] = -g_0_y_y_zzz[k] + g_y_0_y_zzz[k] - g_y_y_y_zzz[k] * ab_y + g_y_y_y_yzzz[k];
            }

            /// Set up 280-290 components of targeted buffer : cbuffer.data(

            auto g_y_y_yz_xxx = cbuffer.data(df_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_y_yz_xxy = cbuffer.data(df_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_y_yz_xxz = cbuffer.data(df_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_y_yz_xyy = cbuffer.data(df_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_y_yz_xyz = cbuffer.data(df_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_y_yz_xzz = cbuffer.data(df_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_y_yz_yyy = cbuffer.data(df_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_y_yz_yyz = cbuffer.data(df_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_y_yz_yzz = cbuffer.data(df_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_y_yz_zzz = cbuffer.data(df_geom_11_off + 289 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_y_xxx, g_y_y_y_xxxz, g_y_y_y_xxy, g_y_y_y_xxyz, g_y_y_y_xxz, g_y_y_y_xxzz, g_y_y_y_xyy, g_y_y_y_xyyz, g_y_y_y_xyz, g_y_y_y_xyzz, g_y_y_y_xzz, g_y_y_y_xzzz, g_y_y_y_yyy, g_y_y_y_yyyz, g_y_y_y_yyz, g_y_y_y_yyzz, g_y_y_y_yzz, g_y_y_y_yzzz, g_y_y_y_zzz, g_y_y_y_zzzz, g_y_y_yz_xxx, g_y_y_yz_xxy, g_y_y_yz_xxz, g_y_y_yz_xyy, g_y_y_yz_xyz, g_y_y_yz_xzz, g_y_y_yz_yyy, g_y_y_yz_yyz, g_y_y_yz_yzz, g_y_y_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yz_xxx[k] = -g_y_y_y_xxx[k] * ab_z + g_y_y_y_xxxz[k];

                g_y_y_yz_xxy[k] = -g_y_y_y_xxy[k] * ab_z + g_y_y_y_xxyz[k];

                g_y_y_yz_xxz[k] = -g_y_y_y_xxz[k] * ab_z + g_y_y_y_xxzz[k];

                g_y_y_yz_xyy[k] = -g_y_y_y_xyy[k] * ab_z + g_y_y_y_xyyz[k];

                g_y_y_yz_xyz[k] = -g_y_y_y_xyz[k] * ab_z + g_y_y_y_xyzz[k];

                g_y_y_yz_xzz[k] = -g_y_y_y_xzz[k] * ab_z + g_y_y_y_xzzz[k];

                g_y_y_yz_yyy[k] = -g_y_y_y_yyy[k] * ab_z + g_y_y_y_yyyz[k];

                g_y_y_yz_yyz[k] = -g_y_y_y_yyz[k] * ab_z + g_y_y_y_yyzz[k];

                g_y_y_yz_yzz[k] = -g_y_y_y_yzz[k] * ab_z + g_y_y_y_yzzz[k];

                g_y_y_yz_zzz[k] = -g_y_y_y_zzz[k] * ab_z + g_y_y_y_zzzz[k];
            }

            /// Set up 290-300 components of targeted buffer : cbuffer.data(

            auto g_y_y_zz_xxx = cbuffer.data(df_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_y_zz_xxy = cbuffer.data(df_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_y_zz_xxz = cbuffer.data(df_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_y_zz_xyy = cbuffer.data(df_geom_11_off + 293 * ccomps * dcomps);

            auto g_y_y_zz_xyz = cbuffer.data(df_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_y_zz_xzz = cbuffer.data(df_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_y_zz_yyy = cbuffer.data(df_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_y_zz_yyz = cbuffer.data(df_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_y_zz_yzz = cbuffer.data(df_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_y_zz_zzz = cbuffer.data(df_geom_11_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_z_xxx, g_y_y_z_xxxz, g_y_y_z_xxy, g_y_y_z_xxyz, g_y_y_z_xxz, g_y_y_z_xxzz, g_y_y_z_xyy, g_y_y_z_xyyz, g_y_y_z_xyz, g_y_y_z_xyzz, g_y_y_z_xzz, g_y_y_z_xzzz, g_y_y_z_yyy, g_y_y_z_yyyz, g_y_y_z_yyz, g_y_y_z_yyzz, g_y_y_z_yzz, g_y_y_z_yzzz, g_y_y_z_zzz, g_y_y_z_zzzz, g_y_y_zz_xxx, g_y_y_zz_xxy, g_y_y_zz_xxz, g_y_y_zz_xyy, g_y_y_zz_xyz, g_y_y_zz_xzz, g_y_y_zz_yyy, g_y_y_zz_yyz, g_y_y_zz_yzz, g_y_y_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zz_xxx[k] = -g_y_y_z_xxx[k] * ab_z + g_y_y_z_xxxz[k];

                g_y_y_zz_xxy[k] = -g_y_y_z_xxy[k] * ab_z + g_y_y_z_xxyz[k];

                g_y_y_zz_xxz[k] = -g_y_y_z_xxz[k] * ab_z + g_y_y_z_xxzz[k];

                g_y_y_zz_xyy[k] = -g_y_y_z_xyy[k] * ab_z + g_y_y_z_xyyz[k];

                g_y_y_zz_xyz[k] = -g_y_y_z_xyz[k] * ab_z + g_y_y_z_xyzz[k];

                g_y_y_zz_xzz[k] = -g_y_y_z_xzz[k] * ab_z + g_y_y_z_xzzz[k];

                g_y_y_zz_yyy[k] = -g_y_y_z_yyy[k] * ab_z + g_y_y_z_yyyz[k];

                g_y_y_zz_yyz[k] = -g_y_y_z_yyz[k] * ab_z + g_y_y_z_yyzz[k];

                g_y_y_zz_yzz[k] = -g_y_y_z_yzz[k] * ab_z + g_y_y_z_yzzz[k];

                g_y_y_zz_zzz[k] = -g_y_y_z_zzz[k] * ab_z + g_y_y_z_zzzz[k];
            }

            /// Set up 300-310 components of targeted buffer : cbuffer.data(

            auto g_y_z_xx_xxx = cbuffer.data(df_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_z_xx_xxy = cbuffer.data(df_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_z_xx_xxz = cbuffer.data(df_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_z_xx_xyy = cbuffer.data(df_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_z_xx_xyz = cbuffer.data(df_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_z_xx_xzz = cbuffer.data(df_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_z_xx_yyy = cbuffer.data(df_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_z_xx_yyz = cbuffer.data(df_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_z_xx_yzz = cbuffer.data(df_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_z_xx_zzz = cbuffer.data(df_geom_11_off + 309 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_x_xxx, g_y_z_x_xxxx, g_y_z_x_xxxy, g_y_z_x_xxxz, g_y_z_x_xxy, g_y_z_x_xxyy, g_y_z_x_xxyz, g_y_z_x_xxz, g_y_z_x_xxzz, g_y_z_x_xyy, g_y_z_x_xyyy, g_y_z_x_xyyz, g_y_z_x_xyz, g_y_z_x_xyzz, g_y_z_x_xzz, g_y_z_x_xzzz, g_y_z_x_yyy, g_y_z_x_yyz, g_y_z_x_yzz, g_y_z_x_zzz, g_y_z_xx_xxx, g_y_z_xx_xxy, g_y_z_xx_xxz, g_y_z_xx_xyy, g_y_z_xx_xyz, g_y_z_xx_xzz, g_y_z_xx_yyy, g_y_z_xx_yyz, g_y_z_xx_yzz, g_y_z_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xx_xxx[k] = -g_y_z_x_xxx[k] * ab_x + g_y_z_x_xxxx[k];

                g_y_z_xx_xxy[k] = -g_y_z_x_xxy[k] * ab_x + g_y_z_x_xxxy[k];

                g_y_z_xx_xxz[k] = -g_y_z_x_xxz[k] * ab_x + g_y_z_x_xxxz[k];

                g_y_z_xx_xyy[k] = -g_y_z_x_xyy[k] * ab_x + g_y_z_x_xxyy[k];

                g_y_z_xx_xyz[k] = -g_y_z_x_xyz[k] * ab_x + g_y_z_x_xxyz[k];

                g_y_z_xx_xzz[k] = -g_y_z_x_xzz[k] * ab_x + g_y_z_x_xxzz[k];

                g_y_z_xx_yyy[k] = -g_y_z_x_yyy[k] * ab_x + g_y_z_x_xyyy[k];

                g_y_z_xx_yyz[k] = -g_y_z_x_yyz[k] * ab_x + g_y_z_x_xyyz[k];

                g_y_z_xx_yzz[k] = -g_y_z_x_yzz[k] * ab_x + g_y_z_x_xyzz[k];

                g_y_z_xx_zzz[k] = -g_y_z_x_zzz[k] * ab_x + g_y_z_x_xzzz[k];
            }

            /// Set up 310-320 components of targeted buffer : cbuffer.data(

            auto g_y_z_xy_xxx = cbuffer.data(df_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_z_xy_xxy = cbuffer.data(df_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_z_xy_xxz = cbuffer.data(df_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_z_xy_xyy = cbuffer.data(df_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_z_xy_xyz = cbuffer.data(df_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_z_xy_xzz = cbuffer.data(df_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_z_xy_yyy = cbuffer.data(df_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_z_xy_yyz = cbuffer.data(df_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_z_xy_yzz = cbuffer.data(df_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_z_xy_zzz = cbuffer.data(df_geom_11_off + 319 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xy_xxx, g_y_z_xy_xxy, g_y_z_xy_xxz, g_y_z_xy_xyy, g_y_z_xy_xyz, g_y_z_xy_xzz, g_y_z_xy_yyy, g_y_z_xy_yyz, g_y_z_xy_yzz, g_y_z_xy_zzz, g_y_z_y_xxx, g_y_z_y_xxxx, g_y_z_y_xxxy, g_y_z_y_xxxz, g_y_z_y_xxy, g_y_z_y_xxyy, g_y_z_y_xxyz, g_y_z_y_xxz, g_y_z_y_xxzz, g_y_z_y_xyy, g_y_z_y_xyyy, g_y_z_y_xyyz, g_y_z_y_xyz, g_y_z_y_xyzz, g_y_z_y_xzz, g_y_z_y_xzzz, g_y_z_y_yyy, g_y_z_y_yyz, g_y_z_y_yzz, g_y_z_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xy_xxx[k] = -g_y_z_y_xxx[k] * ab_x + g_y_z_y_xxxx[k];

                g_y_z_xy_xxy[k] = -g_y_z_y_xxy[k] * ab_x + g_y_z_y_xxxy[k];

                g_y_z_xy_xxz[k] = -g_y_z_y_xxz[k] * ab_x + g_y_z_y_xxxz[k];

                g_y_z_xy_xyy[k] = -g_y_z_y_xyy[k] * ab_x + g_y_z_y_xxyy[k];

                g_y_z_xy_xyz[k] = -g_y_z_y_xyz[k] * ab_x + g_y_z_y_xxyz[k];

                g_y_z_xy_xzz[k] = -g_y_z_y_xzz[k] * ab_x + g_y_z_y_xxzz[k];

                g_y_z_xy_yyy[k] = -g_y_z_y_yyy[k] * ab_x + g_y_z_y_xyyy[k];

                g_y_z_xy_yyz[k] = -g_y_z_y_yyz[k] * ab_x + g_y_z_y_xyyz[k];

                g_y_z_xy_yzz[k] = -g_y_z_y_yzz[k] * ab_x + g_y_z_y_xyzz[k];

                g_y_z_xy_zzz[k] = -g_y_z_y_zzz[k] * ab_x + g_y_z_y_xzzz[k];
            }

            /// Set up 320-330 components of targeted buffer : cbuffer.data(

            auto g_y_z_xz_xxx = cbuffer.data(df_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_z_xz_xxy = cbuffer.data(df_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_z_xz_xxz = cbuffer.data(df_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_z_xz_xyy = cbuffer.data(df_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_z_xz_xyz = cbuffer.data(df_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_z_xz_xzz = cbuffer.data(df_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_z_xz_yyy = cbuffer.data(df_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_z_xz_yyz = cbuffer.data(df_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_z_xz_yzz = cbuffer.data(df_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_z_xz_zzz = cbuffer.data(df_geom_11_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xz_xxx, g_y_z_xz_xxy, g_y_z_xz_xxz, g_y_z_xz_xyy, g_y_z_xz_xyz, g_y_z_xz_xzz, g_y_z_xz_yyy, g_y_z_xz_yyz, g_y_z_xz_yzz, g_y_z_xz_zzz, g_y_z_z_xxx, g_y_z_z_xxxx, g_y_z_z_xxxy, g_y_z_z_xxxz, g_y_z_z_xxy, g_y_z_z_xxyy, g_y_z_z_xxyz, g_y_z_z_xxz, g_y_z_z_xxzz, g_y_z_z_xyy, g_y_z_z_xyyy, g_y_z_z_xyyz, g_y_z_z_xyz, g_y_z_z_xyzz, g_y_z_z_xzz, g_y_z_z_xzzz, g_y_z_z_yyy, g_y_z_z_yyz, g_y_z_z_yzz, g_y_z_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xz_xxx[k] = -g_y_z_z_xxx[k] * ab_x + g_y_z_z_xxxx[k];

                g_y_z_xz_xxy[k] = -g_y_z_z_xxy[k] * ab_x + g_y_z_z_xxxy[k];

                g_y_z_xz_xxz[k] = -g_y_z_z_xxz[k] * ab_x + g_y_z_z_xxxz[k];

                g_y_z_xz_xyy[k] = -g_y_z_z_xyy[k] * ab_x + g_y_z_z_xxyy[k];

                g_y_z_xz_xyz[k] = -g_y_z_z_xyz[k] * ab_x + g_y_z_z_xxyz[k];

                g_y_z_xz_xzz[k] = -g_y_z_z_xzz[k] * ab_x + g_y_z_z_xxzz[k];

                g_y_z_xz_yyy[k] = -g_y_z_z_yyy[k] * ab_x + g_y_z_z_xyyy[k];

                g_y_z_xz_yyz[k] = -g_y_z_z_yyz[k] * ab_x + g_y_z_z_xyyz[k];

                g_y_z_xz_yzz[k] = -g_y_z_z_yzz[k] * ab_x + g_y_z_z_xyzz[k];

                g_y_z_xz_zzz[k] = -g_y_z_z_zzz[k] * ab_x + g_y_z_z_xzzz[k];
            }

            /// Set up 330-340 components of targeted buffer : cbuffer.data(

            auto g_y_z_yy_xxx = cbuffer.data(df_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_z_yy_xxy = cbuffer.data(df_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_z_yy_xxz = cbuffer.data(df_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_z_yy_xyy = cbuffer.data(df_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_z_yy_xyz = cbuffer.data(df_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_z_yy_xzz = cbuffer.data(df_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_z_yy_yyy = cbuffer.data(df_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_z_yy_yyz = cbuffer.data(df_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_z_yy_yzz = cbuffer.data(df_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_z_yy_zzz = cbuffer.data(df_geom_11_off + 339 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_y_xxx, g_0_z_y_xxy, g_0_z_y_xxz, g_0_z_y_xyy, g_0_z_y_xyz, g_0_z_y_xzz, g_0_z_y_yyy, g_0_z_y_yyz, g_0_z_y_yzz, g_0_z_y_zzz, g_y_z_y_xxx, g_y_z_y_xxxy, g_y_z_y_xxy, g_y_z_y_xxyy, g_y_z_y_xxyz, g_y_z_y_xxz, g_y_z_y_xyy, g_y_z_y_xyyy, g_y_z_y_xyyz, g_y_z_y_xyz, g_y_z_y_xyzz, g_y_z_y_xzz, g_y_z_y_yyy, g_y_z_y_yyyy, g_y_z_y_yyyz, g_y_z_y_yyz, g_y_z_y_yyzz, g_y_z_y_yzz, g_y_z_y_yzzz, g_y_z_y_zzz, g_y_z_yy_xxx, g_y_z_yy_xxy, g_y_z_yy_xxz, g_y_z_yy_xyy, g_y_z_yy_xyz, g_y_z_yy_xzz, g_y_z_yy_yyy, g_y_z_yy_yyz, g_y_z_yy_yzz, g_y_z_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yy_xxx[k] = -g_0_z_y_xxx[k] - g_y_z_y_xxx[k] * ab_y + g_y_z_y_xxxy[k];

                g_y_z_yy_xxy[k] = -g_0_z_y_xxy[k] - g_y_z_y_xxy[k] * ab_y + g_y_z_y_xxyy[k];

                g_y_z_yy_xxz[k] = -g_0_z_y_xxz[k] - g_y_z_y_xxz[k] * ab_y + g_y_z_y_xxyz[k];

                g_y_z_yy_xyy[k] = -g_0_z_y_xyy[k] - g_y_z_y_xyy[k] * ab_y + g_y_z_y_xyyy[k];

                g_y_z_yy_xyz[k] = -g_0_z_y_xyz[k] - g_y_z_y_xyz[k] * ab_y + g_y_z_y_xyyz[k];

                g_y_z_yy_xzz[k] = -g_0_z_y_xzz[k] - g_y_z_y_xzz[k] * ab_y + g_y_z_y_xyzz[k];

                g_y_z_yy_yyy[k] = -g_0_z_y_yyy[k] - g_y_z_y_yyy[k] * ab_y + g_y_z_y_yyyy[k];

                g_y_z_yy_yyz[k] = -g_0_z_y_yyz[k] - g_y_z_y_yyz[k] * ab_y + g_y_z_y_yyyz[k];

                g_y_z_yy_yzz[k] = -g_0_z_y_yzz[k] - g_y_z_y_yzz[k] * ab_y + g_y_z_y_yyzz[k];

                g_y_z_yy_zzz[k] = -g_0_z_y_zzz[k] - g_y_z_y_zzz[k] * ab_y + g_y_z_y_yzzz[k];
            }

            /// Set up 340-350 components of targeted buffer : cbuffer.data(

            auto g_y_z_yz_xxx = cbuffer.data(df_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_z_yz_xxy = cbuffer.data(df_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_z_yz_xxz = cbuffer.data(df_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_z_yz_xyy = cbuffer.data(df_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_z_yz_xyz = cbuffer.data(df_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_z_yz_xzz = cbuffer.data(df_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_z_yz_yyy = cbuffer.data(df_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_z_yz_yyz = cbuffer.data(df_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_z_yz_yzz = cbuffer.data(df_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_z_yz_zzz = cbuffer.data(df_geom_11_off + 349 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxx, g_0_z_z_xxy, g_0_z_z_xxz, g_0_z_z_xyy, g_0_z_z_xyz, g_0_z_z_xzz, g_0_z_z_yyy, g_0_z_z_yyz, g_0_z_z_yzz, g_0_z_z_zzz, g_y_z_yz_xxx, g_y_z_yz_xxy, g_y_z_yz_xxz, g_y_z_yz_xyy, g_y_z_yz_xyz, g_y_z_yz_xzz, g_y_z_yz_yyy, g_y_z_yz_yyz, g_y_z_yz_yzz, g_y_z_yz_zzz, g_y_z_z_xxx, g_y_z_z_xxxy, g_y_z_z_xxy, g_y_z_z_xxyy, g_y_z_z_xxyz, g_y_z_z_xxz, g_y_z_z_xyy, g_y_z_z_xyyy, g_y_z_z_xyyz, g_y_z_z_xyz, g_y_z_z_xyzz, g_y_z_z_xzz, g_y_z_z_yyy, g_y_z_z_yyyy, g_y_z_z_yyyz, g_y_z_z_yyz, g_y_z_z_yyzz, g_y_z_z_yzz, g_y_z_z_yzzz, g_y_z_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yz_xxx[k] = -g_0_z_z_xxx[k] - g_y_z_z_xxx[k] * ab_y + g_y_z_z_xxxy[k];

                g_y_z_yz_xxy[k] = -g_0_z_z_xxy[k] - g_y_z_z_xxy[k] * ab_y + g_y_z_z_xxyy[k];

                g_y_z_yz_xxz[k] = -g_0_z_z_xxz[k] - g_y_z_z_xxz[k] * ab_y + g_y_z_z_xxyz[k];

                g_y_z_yz_xyy[k] = -g_0_z_z_xyy[k] - g_y_z_z_xyy[k] * ab_y + g_y_z_z_xyyy[k];

                g_y_z_yz_xyz[k] = -g_0_z_z_xyz[k] - g_y_z_z_xyz[k] * ab_y + g_y_z_z_xyyz[k];

                g_y_z_yz_xzz[k] = -g_0_z_z_xzz[k] - g_y_z_z_xzz[k] * ab_y + g_y_z_z_xyzz[k];

                g_y_z_yz_yyy[k] = -g_0_z_z_yyy[k] - g_y_z_z_yyy[k] * ab_y + g_y_z_z_yyyy[k];

                g_y_z_yz_yyz[k] = -g_0_z_z_yyz[k] - g_y_z_z_yyz[k] * ab_y + g_y_z_z_yyyz[k];

                g_y_z_yz_yzz[k] = -g_0_z_z_yzz[k] - g_y_z_z_yzz[k] * ab_y + g_y_z_z_yyzz[k];

                g_y_z_yz_zzz[k] = -g_0_z_z_zzz[k] - g_y_z_z_zzz[k] * ab_y + g_y_z_z_yzzz[k];
            }

            /// Set up 350-360 components of targeted buffer : cbuffer.data(

            auto g_y_z_zz_xxx = cbuffer.data(df_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_z_zz_xxy = cbuffer.data(df_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_z_zz_xxz = cbuffer.data(df_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_z_zz_xyy = cbuffer.data(df_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_z_zz_xyz = cbuffer.data(df_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_z_zz_xzz = cbuffer.data(df_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_z_zz_yyy = cbuffer.data(df_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_z_zz_yyz = cbuffer.data(df_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_z_zz_yzz = cbuffer.data(df_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_z_zz_zzz = cbuffer.data(df_geom_11_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxx, g_y_0_z_xxy, g_y_0_z_xxz, g_y_0_z_xyy, g_y_0_z_xyz, g_y_0_z_xzz, g_y_0_z_yyy, g_y_0_z_yyz, g_y_0_z_yzz, g_y_0_z_zzz, g_y_z_z_xxx, g_y_z_z_xxxz, g_y_z_z_xxy, g_y_z_z_xxyz, g_y_z_z_xxz, g_y_z_z_xxzz, g_y_z_z_xyy, g_y_z_z_xyyz, g_y_z_z_xyz, g_y_z_z_xyzz, g_y_z_z_xzz, g_y_z_z_xzzz, g_y_z_z_yyy, g_y_z_z_yyyz, g_y_z_z_yyz, g_y_z_z_yyzz, g_y_z_z_yzz, g_y_z_z_yzzz, g_y_z_z_zzz, g_y_z_z_zzzz, g_y_z_zz_xxx, g_y_z_zz_xxy, g_y_z_zz_xxz, g_y_z_zz_xyy, g_y_z_zz_xyz, g_y_z_zz_xzz, g_y_z_zz_yyy, g_y_z_zz_yyz, g_y_z_zz_yzz, g_y_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zz_xxx[k] = g_y_0_z_xxx[k] - g_y_z_z_xxx[k] * ab_z + g_y_z_z_xxxz[k];

                g_y_z_zz_xxy[k] = g_y_0_z_xxy[k] - g_y_z_z_xxy[k] * ab_z + g_y_z_z_xxyz[k];

                g_y_z_zz_xxz[k] = g_y_0_z_xxz[k] - g_y_z_z_xxz[k] * ab_z + g_y_z_z_xxzz[k];

                g_y_z_zz_xyy[k] = g_y_0_z_xyy[k] - g_y_z_z_xyy[k] * ab_z + g_y_z_z_xyyz[k];

                g_y_z_zz_xyz[k] = g_y_0_z_xyz[k] - g_y_z_z_xyz[k] * ab_z + g_y_z_z_xyzz[k];

                g_y_z_zz_xzz[k] = g_y_0_z_xzz[k] - g_y_z_z_xzz[k] * ab_z + g_y_z_z_xzzz[k];

                g_y_z_zz_yyy[k] = g_y_0_z_yyy[k] - g_y_z_z_yyy[k] * ab_z + g_y_z_z_yyyz[k];

                g_y_z_zz_yyz[k] = g_y_0_z_yyz[k] - g_y_z_z_yyz[k] * ab_z + g_y_z_z_yyzz[k];

                g_y_z_zz_yzz[k] = g_y_0_z_yzz[k] - g_y_z_z_yzz[k] * ab_z + g_y_z_z_yzzz[k];

                g_y_z_zz_zzz[k] = g_y_0_z_zzz[k] - g_y_z_z_zzz[k] * ab_z + g_y_z_z_zzzz[k];
            }

            /// Set up 360-370 components of targeted buffer : cbuffer.data(

            auto g_z_x_xx_xxx = cbuffer.data(df_geom_11_off + 360 * ccomps * dcomps);

            auto g_z_x_xx_xxy = cbuffer.data(df_geom_11_off + 361 * ccomps * dcomps);

            auto g_z_x_xx_xxz = cbuffer.data(df_geom_11_off + 362 * ccomps * dcomps);

            auto g_z_x_xx_xyy = cbuffer.data(df_geom_11_off + 363 * ccomps * dcomps);

            auto g_z_x_xx_xyz = cbuffer.data(df_geom_11_off + 364 * ccomps * dcomps);

            auto g_z_x_xx_xzz = cbuffer.data(df_geom_11_off + 365 * ccomps * dcomps);

            auto g_z_x_xx_yyy = cbuffer.data(df_geom_11_off + 366 * ccomps * dcomps);

            auto g_z_x_xx_yyz = cbuffer.data(df_geom_11_off + 367 * ccomps * dcomps);

            auto g_z_x_xx_yzz = cbuffer.data(df_geom_11_off + 368 * ccomps * dcomps);

            auto g_z_x_xx_zzz = cbuffer.data(df_geom_11_off + 369 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_xxx, g_z_0_x_xxy, g_z_0_x_xxz, g_z_0_x_xyy, g_z_0_x_xyz, g_z_0_x_xzz, g_z_0_x_yyy, g_z_0_x_yyz, g_z_0_x_yzz, g_z_0_x_zzz, g_z_x_x_xxx, g_z_x_x_xxxx, g_z_x_x_xxxy, g_z_x_x_xxxz, g_z_x_x_xxy, g_z_x_x_xxyy, g_z_x_x_xxyz, g_z_x_x_xxz, g_z_x_x_xxzz, g_z_x_x_xyy, g_z_x_x_xyyy, g_z_x_x_xyyz, g_z_x_x_xyz, g_z_x_x_xyzz, g_z_x_x_xzz, g_z_x_x_xzzz, g_z_x_x_yyy, g_z_x_x_yyz, g_z_x_x_yzz, g_z_x_x_zzz, g_z_x_xx_xxx, g_z_x_xx_xxy, g_z_x_xx_xxz, g_z_x_xx_xyy, g_z_x_xx_xyz, g_z_x_xx_xzz, g_z_x_xx_yyy, g_z_x_xx_yyz, g_z_x_xx_yzz, g_z_x_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xx_xxx[k] = g_z_0_x_xxx[k] - g_z_x_x_xxx[k] * ab_x + g_z_x_x_xxxx[k];

                g_z_x_xx_xxy[k] = g_z_0_x_xxy[k] - g_z_x_x_xxy[k] * ab_x + g_z_x_x_xxxy[k];

                g_z_x_xx_xxz[k] = g_z_0_x_xxz[k] - g_z_x_x_xxz[k] * ab_x + g_z_x_x_xxxz[k];

                g_z_x_xx_xyy[k] = g_z_0_x_xyy[k] - g_z_x_x_xyy[k] * ab_x + g_z_x_x_xxyy[k];

                g_z_x_xx_xyz[k] = g_z_0_x_xyz[k] - g_z_x_x_xyz[k] * ab_x + g_z_x_x_xxyz[k];

                g_z_x_xx_xzz[k] = g_z_0_x_xzz[k] - g_z_x_x_xzz[k] * ab_x + g_z_x_x_xxzz[k];

                g_z_x_xx_yyy[k] = g_z_0_x_yyy[k] - g_z_x_x_yyy[k] * ab_x + g_z_x_x_xyyy[k];

                g_z_x_xx_yyz[k] = g_z_0_x_yyz[k] - g_z_x_x_yyz[k] * ab_x + g_z_x_x_xyyz[k];

                g_z_x_xx_yzz[k] = g_z_0_x_yzz[k] - g_z_x_x_yzz[k] * ab_x + g_z_x_x_xyzz[k];

                g_z_x_xx_zzz[k] = g_z_0_x_zzz[k] - g_z_x_x_zzz[k] * ab_x + g_z_x_x_xzzz[k];
            }

            /// Set up 370-380 components of targeted buffer : cbuffer.data(

            auto g_z_x_xy_xxx = cbuffer.data(df_geom_11_off + 370 * ccomps * dcomps);

            auto g_z_x_xy_xxy = cbuffer.data(df_geom_11_off + 371 * ccomps * dcomps);

            auto g_z_x_xy_xxz = cbuffer.data(df_geom_11_off + 372 * ccomps * dcomps);

            auto g_z_x_xy_xyy = cbuffer.data(df_geom_11_off + 373 * ccomps * dcomps);

            auto g_z_x_xy_xyz = cbuffer.data(df_geom_11_off + 374 * ccomps * dcomps);

            auto g_z_x_xy_xzz = cbuffer.data(df_geom_11_off + 375 * ccomps * dcomps);

            auto g_z_x_xy_yyy = cbuffer.data(df_geom_11_off + 376 * ccomps * dcomps);

            auto g_z_x_xy_yyz = cbuffer.data(df_geom_11_off + 377 * ccomps * dcomps);

            auto g_z_x_xy_yzz = cbuffer.data(df_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_x_xy_zzz = cbuffer.data(df_geom_11_off + 379 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_x_xxx, g_z_x_x_xxxy, g_z_x_x_xxy, g_z_x_x_xxyy, g_z_x_x_xxyz, g_z_x_x_xxz, g_z_x_x_xyy, g_z_x_x_xyyy, g_z_x_x_xyyz, g_z_x_x_xyz, g_z_x_x_xyzz, g_z_x_x_xzz, g_z_x_x_yyy, g_z_x_x_yyyy, g_z_x_x_yyyz, g_z_x_x_yyz, g_z_x_x_yyzz, g_z_x_x_yzz, g_z_x_x_yzzz, g_z_x_x_zzz, g_z_x_xy_xxx, g_z_x_xy_xxy, g_z_x_xy_xxz, g_z_x_xy_xyy, g_z_x_xy_xyz, g_z_x_xy_xzz, g_z_x_xy_yyy, g_z_x_xy_yyz, g_z_x_xy_yzz, g_z_x_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xy_xxx[k] = -g_z_x_x_xxx[k] * ab_y + g_z_x_x_xxxy[k];

                g_z_x_xy_xxy[k] = -g_z_x_x_xxy[k] * ab_y + g_z_x_x_xxyy[k];

                g_z_x_xy_xxz[k] = -g_z_x_x_xxz[k] * ab_y + g_z_x_x_xxyz[k];

                g_z_x_xy_xyy[k] = -g_z_x_x_xyy[k] * ab_y + g_z_x_x_xyyy[k];

                g_z_x_xy_xyz[k] = -g_z_x_x_xyz[k] * ab_y + g_z_x_x_xyyz[k];

                g_z_x_xy_xzz[k] = -g_z_x_x_xzz[k] * ab_y + g_z_x_x_xyzz[k];

                g_z_x_xy_yyy[k] = -g_z_x_x_yyy[k] * ab_y + g_z_x_x_yyyy[k];

                g_z_x_xy_yyz[k] = -g_z_x_x_yyz[k] * ab_y + g_z_x_x_yyyz[k];

                g_z_x_xy_yzz[k] = -g_z_x_x_yzz[k] * ab_y + g_z_x_x_yyzz[k];

                g_z_x_xy_zzz[k] = -g_z_x_x_zzz[k] * ab_y + g_z_x_x_yzzz[k];
            }

            /// Set up 380-390 components of targeted buffer : cbuffer.data(

            auto g_z_x_xz_xxx = cbuffer.data(df_geom_11_off + 380 * ccomps * dcomps);

            auto g_z_x_xz_xxy = cbuffer.data(df_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_x_xz_xxz = cbuffer.data(df_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_x_xz_xyy = cbuffer.data(df_geom_11_off + 383 * ccomps * dcomps);

            auto g_z_x_xz_xyz = cbuffer.data(df_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_x_xz_xzz = cbuffer.data(df_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_x_xz_yyy = cbuffer.data(df_geom_11_off + 386 * ccomps * dcomps);

            auto g_z_x_xz_yyz = cbuffer.data(df_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_x_xz_yzz = cbuffer.data(df_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_x_xz_zzz = cbuffer.data(df_geom_11_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz, g_z_x_xz_xxx, g_z_x_xz_xxy, g_z_x_xz_xxz, g_z_x_xz_xyy, g_z_x_xz_xyz, g_z_x_xz_xzz, g_z_x_xz_yyy, g_z_x_xz_yyz, g_z_x_xz_yzz, g_z_x_xz_zzz, g_z_x_z_xxx, g_z_x_z_xxxx, g_z_x_z_xxxy, g_z_x_z_xxxz, g_z_x_z_xxy, g_z_x_z_xxyy, g_z_x_z_xxyz, g_z_x_z_xxz, g_z_x_z_xxzz, g_z_x_z_xyy, g_z_x_z_xyyy, g_z_x_z_xyyz, g_z_x_z_xyz, g_z_x_z_xyzz, g_z_x_z_xzz, g_z_x_z_xzzz, g_z_x_z_yyy, g_z_x_z_yyz, g_z_x_z_yzz, g_z_x_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xz_xxx[k] = g_z_0_z_xxx[k] - g_z_x_z_xxx[k] * ab_x + g_z_x_z_xxxx[k];

                g_z_x_xz_xxy[k] = g_z_0_z_xxy[k] - g_z_x_z_xxy[k] * ab_x + g_z_x_z_xxxy[k];

                g_z_x_xz_xxz[k] = g_z_0_z_xxz[k] - g_z_x_z_xxz[k] * ab_x + g_z_x_z_xxxz[k];

                g_z_x_xz_xyy[k] = g_z_0_z_xyy[k] - g_z_x_z_xyy[k] * ab_x + g_z_x_z_xxyy[k];

                g_z_x_xz_xyz[k] = g_z_0_z_xyz[k] - g_z_x_z_xyz[k] * ab_x + g_z_x_z_xxyz[k];

                g_z_x_xz_xzz[k] = g_z_0_z_xzz[k] - g_z_x_z_xzz[k] * ab_x + g_z_x_z_xxzz[k];

                g_z_x_xz_yyy[k] = g_z_0_z_yyy[k] - g_z_x_z_yyy[k] * ab_x + g_z_x_z_xyyy[k];

                g_z_x_xz_yyz[k] = g_z_0_z_yyz[k] - g_z_x_z_yyz[k] * ab_x + g_z_x_z_xyyz[k];

                g_z_x_xz_yzz[k] = g_z_0_z_yzz[k] - g_z_x_z_yzz[k] * ab_x + g_z_x_z_xyzz[k];

                g_z_x_xz_zzz[k] = g_z_0_z_zzz[k] - g_z_x_z_zzz[k] * ab_x + g_z_x_z_xzzz[k];
            }

            /// Set up 390-400 components of targeted buffer : cbuffer.data(

            auto g_z_x_yy_xxx = cbuffer.data(df_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_x_yy_xxy = cbuffer.data(df_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_x_yy_xxz = cbuffer.data(df_geom_11_off + 392 * ccomps * dcomps);

            auto g_z_x_yy_xyy = cbuffer.data(df_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_x_yy_xyz = cbuffer.data(df_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_x_yy_xzz = cbuffer.data(df_geom_11_off + 395 * ccomps * dcomps);

            auto g_z_x_yy_yyy = cbuffer.data(df_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_x_yy_yyz = cbuffer.data(df_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_x_yy_yzz = cbuffer.data(df_geom_11_off + 398 * ccomps * dcomps);

            auto g_z_x_yy_zzz = cbuffer.data(df_geom_11_off + 399 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_y_xxx, g_z_x_y_xxxy, g_z_x_y_xxy, g_z_x_y_xxyy, g_z_x_y_xxyz, g_z_x_y_xxz, g_z_x_y_xyy, g_z_x_y_xyyy, g_z_x_y_xyyz, g_z_x_y_xyz, g_z_x_y_xyzz, g_z_x_y_xzz, g_z_x_y_yyy, g_z_x_y_yyyy, g_z_x_y_yyyz, g_z_x_y_yyz, g_z_x_y_yyzz, g_z_x_y_yzz, g_z_x_y_yzzz, g_z_x_y_zzz, g_z_x_yy_xxx, g_z_x_yy_xxy, g_z_x_yy_xxz, g_z_x_yy_xyy, g_z_x_yy_xyz, g_z_x_yy_xzz, g_z_x_yy_yyy, g_z_x_yy_yyz, g_z_x_yy_yzz, g_z_x_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yy_xxx[k] = -g_z_x_y_xxx[k] * ab_y + g_z_x_y_xxxy[k];

                g_z_x_yy_xxy[k] = -g_z_x_y_xxy[k] * ab_y + g_z_x_y_xxyy[k];

                g_z_x_yy_xxz[k] = -g_z_x_y_xxz[k] * ab_y + g_z_x_y_xxyz[k];

                g_z_x_yy_xyy[k] = -g_z_x_y_xyy[k] * ab_y + g_z_x_y_xyyy[k];

                g_z_x_yy_xyz[k] = -g_z_x_y_xyz[k] * ab_y + g_z_x_y_xyyz[k];

                g_z_x_yy_xzz[k] = -g_z_x_y_xzz[k] * ab_y + g_z_x_y_xyzz[k];

                g_z_x_yy_yyy[k] = -g_z_x_y_yyy[k] * ab_y + g_z_x_y_yyyy[k];

                g_z_x_yy_yyz[k] = -g_z_x_y_yyz[k] * ab_y + g_z_x_y_yyyz[k];

                g_z_x_yy_yzz[k] = -g_z_x_y_yzz[k] * ab_y + g_z_x_y_yyzz[k];

                g_z_x_yy_zzz[k] = -g_z_x_y_zzz[k] * ab_y + g_z_x_y_yzzz[k];
            }

            /// Set up 400-410 components of targeted buffer : cbuffer.data(

            auto g_z_x_yz_xxx = cbuffer.data(df_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_x_yz_xxy = cbuffer.data(df_geom_11_off + 401 * ccomps * dcomps);

            auto g_z_x_yz_xxz = cbuffer.data(df_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_x_yz_xyy = cbuffer.data(df_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_x_yz_xyz = cbuffer.data(df_geom_11_off + 404 * ccomps * dcomps);

            auto g_z_x_yz_xzz = cbuffer.data(df_geom_11_off + 405 * ccomps * dcomps);

            auto g_z_x_yz_yyy = cbuffer.data(df_geom_11_off + 406 * ccomps * dcomps);

            auto g_z_x_yz_yyz = cbuffer.data(df_geom_11_off + 407 * ccomps * dcomps);

            auto g_z_x_yz_yzz = cbuffer.data(df_geom_11_off + 408 * ccomps * dcomps);

            auto g_z_x_yz_zzz = cbuffer.data(df_geom_11_off + 409 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yz_xxx, g_z_x_yz_xxy, g_z_x_yz_xxz, g_z_x_yz_xyy, g_z_x_yz_xyz, g_z_x_yz_xzz, g_z_x_yz_yyy, g_z_x_yz_yyz, g_z_x_yz_yzz, g_z_x_yz_zzz, g_z_x_z_xxx, g_z_x_z_xxxy, g_z_x_z_xxy, g_z_x_z_xxyy, g_z_x_z_xxyz, g_z_x_z_xxz, g_z_x_z_xyy, g_z_x_z_xyyy, g_z_x_z_xyyz, g_z_x_z_xyz, g_z_x_z_xyzz, g_z_x_z_xzz, g_z_x_z_yyy, g_z_x_z_yyyy, g_z_x_z_yyyz, g_z_x_z_yyz, g_z_x_z_yyzz, g_z_x_z_yzz, g_z_x_z_yzzz, g_z_x_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yz_xxx[k] = -g_z_x_z_xxx[k] * ab_y + g_z_x_z_xxxy[k];

                g_z_x_yz_xxy[k] = -g_z_x_z_xxy[k] * ab_y + g_z_x_z_xxyy[k];

                g_z_x_yz_xxz[k] = -g_z_x_z_xxz[k] * ab_y + g_z_x_z_xxyz[k];

                g_z_x_yz_xyy[k] = -g_z_x_z_xyy[k] * ab_y + g_z_x_z_xyyy[k];

                g_z_x_yz_xyz[k] = -g_z_x_z_xyz[k] * ab_y + g_z_x_z_xyyz[k];

                g_z_x_yz_xzz[k] = -g_z_x_z_xzz[k] * ab_y + g_z_x_z_xyzz[k];

                g_z_x_yz_yyy[k] = -g_z_x_z_yyy[k] * ab_y + g_z_x_z_yyyy[k];

                g_z_x_yz_yyz[k] = -g_z_x_z_yyz[k] * ab_y + g_z_x_z_yyyz[k];

                g_z_x_yz_yzz[k] = -g_z_x_z_yzz[k] * ab_y + g_z_x_z_yyzz[k];

                g_z_x_yz_zzz[k] = -g_z_x_z_zzz[k] * ab_y + g_z_x_z_yzzz[k];
            }

            /// Set up 410-420 components of targeted buffer : cbuffer.data(

            auto g_z_x_zz_xxx = cbuffer.data(df_geom_11_off + 410 * ccomps * dcomps);

            auto g_z_x_zz_xxy = cbuffer.data(df_geom_11_off + 411 * ccomps * dcomps);

            auto g_z_x_zz_xxz = cbuffer.data(df_geom_11_off + 412 * ccomps * dcomps);

            auto g_z_x_zz_xyy = cbuffer.data(df_geom_11_off + 413 * ccomps * dcomps);

            auto g_z_x_zz_xyz = cbuffer.data(df_geom_11_off + 414 * ccomps * dcomps);

            auto g_z_x_zz_xzz = cbuffer.data(df_geom_11_off + 415 * ccomps * dcomps);

            auto g_z_x_zz_yyy = cbuffer.data(df_geom_11_off + 416 * ccomps * dcomps);

            auto g_z_x_zz_yyz = cbuffer.data(df_geom_11_off + 417 * ccomps * dcomps);

            auto g_z_x_zz_yzz = cbuffer.data(df_geom_11_off + 418 * ccomps * dcomps);

            auto g_z_x_zz_zzz = cbuffer.data(df_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_z_xxx, g_0_x_z_xxy, g_0_x_z_xxz, g_0_x_z_xyy, g_0_x_z_xyz, g_0_x_z_xzz, g_0_x_z_yyy, g_0_x_z_yyz, g_0_x_z_yzz, g_0_x_z_zzz, g_z_x_z_xxx, g_z_x_z_xxxz, g_z_x_z_xxy, g_z_x_z_xxyz, g_z_x_z_xxz, g_z_x_z_xxzz, g_z_x_z_xyy, g_z_x_z_xyyz, g_z_x_z_xyz, g_z_x_z_xyzz, g_z_x_z_xzz, g_z_x_z_xzzz, g_z_x_z_yyy, g_z_x_z_yyyz, g_z_x_z_yyz, g_z_x_z_yyzz, g_z_x_z_yzz, g_z_x_z_yzzz, g_z_x_z_zzz, g_z_x_z_zzzz, g_z_x_zz_xxx, g_z_x_zz_xxy, g_z_x_zz_xxz, g_z_x_zz_xyy, g_z_x_zz_xyz, g_z_x_zz_xzz, g_z_x_zz_yyy, g_z_x_zz_yyz, g_z_x_zz_yzz, g_z_x_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zz_xxx[k] = -g_0_x_z_xxx[k] - g_z_x_z_xxx[k] * ab_z + g_z_x_z_xxxz[k];

                g_z_x_zz_xxy[k] = -g_0_x_z_xxy[k] - g_z_x_z_xxy[k] * ab_z + g_z_x_z_xxyz[k];

                g_z_x_zz_xxz[k] = -g_0_x_z_xxz[k] - g_z_x_z_xxz[k] * ab_z + g_z_x_z_xxzz[k];

                g_z_x_zz_xyy[k] = -g_0_x_z_xyy[k] - g_z_x_z_xyy[k] * ab_z + g_z_x_z_xyyz[k];

                g_z_x_zz_xyz[k] = -g_0_x_z_xyz[k] - g_z_x_z_xyz[k] * ab_z + g_z_x_z_xyzz[k];

                g_z_x_zz_xzz[k] = -g_0_x_z_xzz[k] - g_z_x_z_xzz[k] * ab_z + g_z_x_z_xzzz[k];

                g_z_x_zz_yyy[k] = -g_0_x_z_yyy[k] - g_z_x_z_yyy[k] * ab_z + g_z_x_z_yyyz[k];

                g_z_x_zz_yyz[k] = -g_0_x_z_yyz[k] - g_z_x_z_yyz[k] * ab_z + g_z_x_z_yyzz[k];

                g_z_x_zz_yzz[k] = -g_0_x_z_yzz[k] - g_z_x_z_yzz[k] * ab_z + g_z_x_z_yzzz[k];

                g_z_x_zz_zzz[k] = -g_0_x_z_zzz[k] - g_z_x_z_zzz[k] * ab_z + g_z_x_z_zzzz[k];
            }

            /// Set up 420-430 components of targeted buffer : cbuffer.data(

            auto g_z_y_xx_xxx = cbuffer.data(df_geom_11_off + 420 * ccomps * dcomps);

            auto g_z_y_xx_xxy = cbuffer.data(df_geom_11_off + 421 * ccomps * dcomps);

            auto g_z_y_xx_xxz = cbuffer.data(df_geom_11_off + 422 * ccomps * dcomps);

            auto g_z_y_xx_xyy = cbuffer.data(df_geom_11_off + 423 * ccomps * dcomps);

            auto g_z_y_xx_xyz = cbuffer.data(df_geom_11_off + 424 * ccomps * dcomps);

            auto g_z_y_xx_xzz = cbuffer.data(df_geom_11_off + 425 * ccomps * dcomps);

            auto g_z_y_xx_yyy = cbuffer.data(df_geom_11_off + 426 * ccomps * dcomps);

            auto g_z_y_xx_yyz = cbuffer.data(df_geom_11_off + 427 * ccomps * dcomps);

            auto g_z_y_xx_yzz = cbuffer.data(df_geom_11_off + 428 * ccomps * dcomps);

            auto g_z_y_xx_zzz = cbuffer.data(df_geom_11_off + 429 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_x_xxx, g_z_y_x_xxxx, g_z_y_x_xxxy, g_z_y_x_xxxz, g_z_y_x_xxy, g_z_y_x_xxyy, g_z_y_x_xxyz, g_z_y_x_xxz, g_z_y_x_xxzz, g_z_y_x_xyy, g_z_y_x_xyyy, g_z_y_x_xyyz, g_z_y_x_xyz, g_z_y_x_xyzz, g_z_y_x_xzz, g_z_y_x_xzzz, g_z_y_x_yyy, g_z_y_x_yyz, g_z_y_x_yzz, g_z_y_x_zzz, g_z_y_xx_xxx, g_z_y_xx_xxy, g_z_y_xx_xxz, g_z_y_xx_xyy, g_z_y_xx_xyz, g_z_y_xx_xzz, g_z_y_xx_yyy, g_z_y_xx_yyz, g_z_y_xx_yzz, g_z_y_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xx_xxx[k] = -g_z_y_x_xxx[k] * ab_x + g_z_y_x_xxxx[k];

                g_z_y_xx_xxy[k] = -g_z_y_x_xxy[k] * ab_x + g_z_y_x_xxxy[k];

                g_z_y_xx_xxz[k] = -g_z_y_x_xxz[k] * ab_x + g_z_y_x_xxxz[k];

                g_z_y_xx_xyy[k] = -g_z_y_x_xyy[k] * ab_x + g_z_y_x_xxyy[k];

                g_z_y_xx_xyz[k] = -g_z_y_x_xyz[k] * ab_x + g_z_y_x_xxyz[k];

                g_z_y_xx_xzz[k] = -g_z_y_x_xzz[k] * ab_x + g_z_y_x_xxzz[k];

                g_z_y_xx_yyy[k] = -g_z_y_x_yyy[k] * ab_x + g_z_y_x_xyyy[k];

                g_z_y_xx_yyz[k] = -g_z_y_x_yyz[k] * ab_x + g_z_y_x_xyyz[k];

                g_z_y_xx_yzz[k] = -g_z_y_x_yzz[k] * ab_x + g_z_y_x_xyzz[k];

                g_z_y_xx_zzz[k] = -g_z_y_x_zzz[k] * ab_x + g_z_y_x_xzzz[k];
            }

            /// Set up 430-440 components of targeted buffer : cbuffer.data(

            auto g_z_y_xy_xxx = cbuffer.data(df_geom_11_off + 430 * ccomps * dcomps);

            auto g_z_y_xy_xxy = cbuffer.data(df_geom_11_off + 431 * ccomps * dcomps);

            auto g_z_y_xy_xxz = cbuffer.data(df_geom_11_off + 432 * ccomps * dcomps);

            auto g_z_y_xy_xyy = cbuffer.data(df_geom_11_off + 433 * ccomps * dcomps);

            auto g_z_y_xy_xyz = cbuffer.data(df_geom_11_off + 434 * ccomps * dcomps);

            auto g_z_y_xy_xzz = cbuffer.data(df_geom_11_off + 435 * ccomps * dcomps);

            auto g_z_y_xy_yyy = cbuffer.data(df_geom_11_off + 436 * ccomps * dcomps);

            auto g_z_y_xy_yyz = cbuffer.data(df_geom_11_off + 437 * ccomps * dcomps);

            auto g_z_y_xy_yzz = cbuffer.data(df_geom_11_off + 438 * ccomps * dcomps);

            auto g_z_y_xy_zzz = cbuffer.data(df_geom_11_off + 439 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xy_xxx, g_z_y_xy_xxy, g_z_y_xy_xxz, g_z_y_xy_xyy, g_z_y_xy_xyz, g_z_y_xy_xzz, g_z_y_xy_yyy, g_z_y_xy_yyz, g_z_y_xy_yzz, g_z_y_xy_zzz, g_z_y_y_xxx, g_z_y_y_xxxx, g_z_y_y_xxxy, g_z_y_y_xxxz, g_z_y_y_xxy, g_z_y_y_xxyy, g_z_y_y_xxyz, g_z_y_y_xxz, g_z_y_y_xxzz, g_z_y_y_xyy, g_z_y_y_xyyy, g_z_y_y_xyyz, g_z_y_y_xyz, g_z_y_y_xyzz, g_z_y_y_xzz, g_z_y_y_xzzz, g_z_y_y_yyy, g_z_y_y_yyz, g_z_y_y_yzz, g_z_y_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xy_xxx[k] = -g_z_y_y_xxx[k] * ab_x + g_z_y_y_xxxx[k];

                g_z_y_xy_xxy[k] = -g_z_y_y_xxy[k] * ab_x + g_z_y_y_xxxy[k];

                g_z_y_xy_xxz[k] = -g_z_y_y_xxz[k] * ab_x + g_z_y_y_xxxz[k];

                g_z_y_xy_xyy[k] = -g_z_y_y_xyy[k] * ab_x + g_z_y_y_xxyy[k];

                g_z_y_xy_xyz[k] = -g_z_y_y_xyz[k] * ab_x + g_z_y_y_xxyz[k];

                g_z_y_xy_xzz[k] = -g_z_y_y_xzz[k] * ab_x + g_z_y_y_xxzz[k];

                g_z_y_xy_yyy[k] = -g_z_y_y_yyy[k] * ab_x + g_z_y_y_xyyy[k];

                g_z_y_xy_yyz[k] = -g_z_y_y_yyz[k] * ab_x + g_z_y_y_xyyz[k];

                g_z_y_xy_yzz[k] = -g_z_y_y_yzz[k] * ab_x + g_z_y_y_xyzz[k];

                g_z_y_xy_zzz[k] = -g_z_y_y_zzz[k] * ab_x + g_z_y_y_xzzz[k];
            }

            /// Set up 440-450 components of targeted buffer : cbuffer.data(

            auto g_z_y_xz_xxx = cbuffer.data(df_geom_11_off + 440 * ccomps * dcomps);

            auto g_z_y_xz_xxy = cbuffer.data(df_geom_11_off + 441 * ccomps * dcomps);

            auto g_z_y_xz_xxz = cbuffer.data(df_geom_11_off + 442 * ccomps * dcomps);

            auto g_z_y_xz_xyy = cbuffer.data(df_geom_11_off + 443 * ccomps * dcomps);

            auto g_z_y_xz_xyz = cbuffer.data(df_geom_11_off + 444 * ccomps * dcomps);

            auto g_z_y_xz_xzz = cbuffer.data(df_geom_11_off + 445 * ccomps * dcomps);

            auto g_z_y_xz_yyy = cbuffer.data(df_geom_11_off + 446 * ccomps * dcomps);

            auto g_z_y_xz_yyz = cbuffer.data(df_geom_11_off + 447 * ccomps * dcomps);

            auto g_z_y_xz_yzz = cbuffer.data(df_geom_11_off + 448 * ccomps * dcomps);

            auto g_z_y_xz_zzz = cbuffer.data(df_geom_11_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xz_xxx, g_z_y_xz_xxy, g_z_y_xz_xxz, g_z_y_xz_xyy, g_z_y_xz_xyz, g_z_y_xz_xzz, g_z_y_xz_yyy, g_z_y_xz_yyz, g_z_y_xz_yzz, g_z_y_xz_zzz, g_z_y_z_xxx, g_z_y_z_xxxx, g_z_y_z_xxxy, g_z_y_z_xxxz, g_z_y_z_xxy, g_z_y_z_xxyy, g_z_y_z_xxyz, g_z_y_z_xxz, g_z_y_z_xxzz, g_z_y_z_xyy, g_z_y_z_xyyy, g_z_y_z_xyyz, g_z_y_z_xyz, g_z_y_z_xyzz, g_z_y_z_xzz, g_z_y_z_xzzz, g_z_y_z_yyy, g_z_y_z_yyz, g_z_y_z_yzz, g_z_y_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xz_xxx[k] = -g_z_y_z_xxx[k] * ab_x + g_z_y_z_xxxx[k];

                g_z_y_xz_xxy[k] = -g_z_y_z_xxy[k] * ab_x + g_z_y_z_xxxy[k];

                g_z_y_xz_xxz[k] = -g_z_y_z_xxz[k] * ab_x + g_z_y_z_xxxz[k];

                g_z_y_xz_xyy[k] = -g_z_y_z_xyy[k] * ab_x + g_z_y_z_xxyy[k];

                g_z_y_xz_xyz[k] = -g_z_y_z_xyz[k] * ab_x + g_z_y_z_xxyz[k];

                g_z_y_xz_xzz[k] = -g_z_y_z_xzz[k] * ab_x + g_z_y_z_xxzz[k];

                g_z_y_xz_yyy[k] = -g_z_y_z_yyy[k] * ab_x + g_z_y_z_xyyy[k];

                g_z_y_xz_yyz[k] = -g_z_y_z_yyz[k] * ab_x + g_z_y_z_xyyz[k];

                g_z_y_xz_yzz[k] = -g_z_y_z_yzz[k] * ab_x + g_z_y_z_xyzz[k];

                g_z_y_xz_zzz[k] = -g_z_y_z_zzz[k] * ab_x + g_z_y_z_xzzz[k];
            }

            /// Set up 450-460 components of targeted buffer : cbuffer.data(

            auto g_z_y_yy_xxx = cbuffer.data(df_geom_11_off + 450 * ccomps * dcomps);

            auto g_z_y_yy_xxy = cbuffer.data(df_geom_11_off + 451 * ccomps * dcomps);

            auto g_z_y_yy_xxz = cbuffer.data(df_geom_11_off + 452 * ccomps * dcomps);

            auto g_z_y_yy_xyy = cbuffer.data(df_geom_11_off + 453 * ccomps * dcomps);

            auto g_z_y_yy_xyz = cbuffer.data(df_geom_11_off + 454 * ccomps * dcomps);

            auto g_z_y_yy_xzz = cbuffer.data(df_geom_11_off + 455 * ccomps * dcomps);

            auto g_z_y_yy_yyy = cbuffer.data(df_geom_11_off + 456 * ccomps * dcomps);

            auto g_z_y_yy_yyz = cbuffer.data(df_geom_11_off + 457 * ccomps * dcomps);

            auto g_z_y_yy_yzz = cbuffer.data(df_geom_11_off + 458 * ccomps * dcomps);

            auto g_z_y_yy_zzz = cbuffer.data(df_geom_11_off + 459 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_xxx, g_z_0_y_xxy, g_z_0_y_xxz, g_z_0_y_xyy, g_z_0_y_xyz, g_z_0_y_xzz, g_z_0_y_yyy, g_z_0_y_yyz, g_z_0_y_yzz, g_z_0_y_zzz, g_z_y_y_xxx, g_z_y_y_xxxy, g_z_y_y_xxy, g_z_y_y_xxyy, g_z_y_y_xxyz, g_z_y_y_xxz, g_z_y_y_xyy, g_z_y_y_xyyy, g_z_y_y_xyyz, g_z_y_y_xyz, g_z_y_y_xyzz, g_z_y_y_xzz, g_z_y_y_yyy, g_z_y_y_yyyy, g_z_y_y_yyyz, g_z_y_y_yyz, g_z_y_y_yyzz, g_z_y_y_yzz, g_z_y_y_yzzz, g_z_y_y_zzz, g_z_y_yy_xxx, g_z_y_yy_xxy, g_z_y_yy_xxz, g_z_y_yy_xyy, g_z_y_yy_xyz, g_z_y_yy_xzz, g_z_y_yy_yyy, g_z_y_yy_yyz, g_z_y_yy_yzz, g_z_y_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yy_xxx[k] = g_z_0_y_xxx[k] - g_z_y_y_xxx[k] * ab_y + g_z_y_y_xxxy[k];

                g_z_y_yy_xxy[k] = g_z_0_y_xxy[k] - g_z_y_y_xxy[k] * ab_y + g_z_y_y_xxyy[k];

                g_z_y_yy_xxz[k] = g_z_0_y_xxz[k] - g_z_y_y_xxz[k] * ab_y + g_z_y_y_xxyz[k];

                g_z_y_yy_xyy[k] = g_z_0_y_xyy[k] - g_z_y_y_xyy[k] * ab_y + g_z_y_y_xyyy[k];

                g_z_y_yy_xyz[k] = g_z_0_y_xyz[k] - g_z_y_y_xyz[k] * ab_y + g_z_y_y_xyyz[k];

                g_z_y_yy_xzz[k] = g_z_0_y_xzz[k] - g_z_y_y_xzz[k] * ab_y + g_z_y_y_xyzz[k];

                g_z_y_yy_yyy[k] = g_z_0_y_yyy[k] - g_z_y_y_yyy[k] * ab_y + g_z_y_y_yyyy[k];

                g_z_y_yy_yyz[k] = g_z_0_y_yyz[k] - g_z_y_y_yyz[k] * ab_y + g_z_y_y_yyyz[k];

                g_z_y_yy_yzz[k] = g_z_0_y_yzz[k] - g_z_y_y_yzz[k] * ab_y + g_z_y_y_yyzz[k];

                g_z_y_yy_zzz[k] = g_z_0_y_zzz[k] - g_z_y_y_zzz[k] * ab_y + g_z_y_y_yzzz[k];
            }

            /// Set up 460-470 components of targeted buffer : cbuffer.data(

            auto g_z_y_yz_xxx = cbuffer.data(df_geom_11_off + 460 * ccomps * dcomps);

            auto g_z_y_yz_xxy = cbuffer.data(df_geom_11_off + 461 * ccomps * dcomps);

            auto g_z_y_yz_xxz = cbuffer.data(df_geom_11_off + 462 * ccomps * dcomps);

            auto g_z_y_yz_xyy = cbuffer.data(df_geom_11_off + 463 * ccomps * dcomps);

            auto g_z_y_yz_xyz = cbuffer.data(df_geom_11_off + 464 * ccomps * dcomps);

            auto g_z_y_yz_xzz = cbuffer.data(df_geom_11_off + 465 * ccomps * dcomps);

            auto g_z_y_yz_yyy = cbuffer.data(df_geom_11_off + 466 * ccomps * dcomps);

            auto g_z_y_yz_yyz = cbuffer.data(df_geom_11_off + 467 * ccomps * dcomps);

            auto g_z_y_yz_yzz = cbuffer.data(df_geom_11_off + 468 * ccomps * dcomps);

            auto g_z_y_yz_zzz = cbuffer.data(df_geom_11_off + 469 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz, g_z_y_yz_xxx, g_z_y_yz_xxy, g_z_y_yz_xxz, g_z_y_yz_xyy, g_z_y_yz_xyz, g_z_y_yz_xzz, g_z_y_yz_yyy, g_z_y_yz_yyz, g_z_y_yz_yzz, g_z_y_yz_zzz, g_z_y_z_xxx, g_z_y_z_xxxy, g_z_y_z_xxy, g_z_y_z_xxyy, g_z_y_z_xxyz, g_z_y_z_xxz, g_z_y_z_xyy, g_z_y_z_xyyy, g_z_y_z_xyyz, g_z_y_z_xyz, g_z_y_z_xyzz, g_z_y_z_xzz, g_z_y_z_yyy, g_z_y_z_yyyy, g_z_y_z_yyyz, g_z_y_z_yyz, g_z_y_z_yyzz, g_z_y_z_yzz, g_z_y_z_yzzz, g_z_y_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yz_xxx[k] = g_z_0_z_xxx[k] - g_z_y_z_xxx[k] * ab_y + g_z_y_z_xxxy[k];

                g_z_y_yz_xxy[k] = g_z_0_z_xxy[k] - g_z_y_z_xxy[k] * ab_y + g_z_y_z_xxyy[k];

                g_z_y_yz_xxz[k] = g_z_0_z_xxz[k] - g_z_y_z_xxz[k] * ab_y + g_z_y_z_xxyz[k];

                g_z_y_yz_xyy[k] = g_z_0_z_xyy[k] - g_z_y_z_xyy[k] * ab_y + g_z_y_z_xyyy[k];

                g_z_y_yz_xyz[k] = g_z_0_z_xyz[k] - g_z_y_z_xyz[k] * ab_y + g_z_y_z_xyyz[k];

                g_z_y_yz_xzz[k] = g_z_0_z_xzz[k] - g_z_y_z_xzz[k] * ab_y + g_z_y_z_xyzz[k];

                g_z_y_yz_yyy[k] = g_z_0_z_yyy[k] - g_z_y_z_yyy[k] * ab_y + g_z_y_z_yyyy[k];

                g_z_y_yz_yyz[k] = g_z_0_z_yyz[k] - g_z_y_z_yyz[k] * ab_y + g_z_y_z_yyyz[k];

                g_z_y_yz_yzz[k] = g_z_0_z_yzz[k] - g_z_y_z_yzz[k] * ab_y + g_z_y_z_yyzz[k];

                g_z_y_yz_zzz[k] = g_z_0_z_zzz[k] - g_z_y_z_zzz[k] * ab_y + g_z_y_z_yzzz[k];
            }

            /// Set up 470-480 components of targeted buffer : cbuffer.data(

            auto g_z_y_zz_xxx = cbuffer.data(df_geom_11_off + 470 * ccomps * dcomps);

            auto g_z_y_zz_xxy = cbuffer.data(df_geom_11_off + 471 * ccomps * dcomps);

            auto g_z_y_zz_xxz = cbuffer.data(df_geom_11_off + 472 * ccomps * dcomps);

            auto g_z_y_zz_xyy = cbuffer.data(df_geom_11_off + 473 * ccomps * dcomps);

            auto g_z_y_zz_xyz = cbuffer.data(df_geom_11_off + 474 * ccomps * dcomps);

            auto g_z_y_zz_xzz = cbuffer.data(df_geom_11_off + 475 * ccomps * dcomps);

            auto g_z_y_zz_yyy = cbuffer.data(df_geom_11_off + 476 * ccomps * dcomps);

            auto g_z_y_zz_yyz = cbuffer.data(df_geom_11_off + 477 * ccomps * dcomps);

            auto g_z_y_zz_yzz = cbuffer.data(df_geom_11_off + 478 * ccomps * dcomps);

            auto g_z_y_zz_zzz = cbuffer.data(df_geom_11_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_z_xxx, g_0_y_z_xxy, g_0_y_z_xxz, g_0_y_z_xyy, g_0_y_z_xyz, g_0_y_z_xzz, g_0_y_z_yyy, g_0_y_z_yyz, g_0_y_z_yzz, g_0_y_z_zzz, g_z_y_z_xxx, g_z_y_z_xxxz, g_z_y_z_xxy, g_z_y_z_xxyz, g_z_y_z_xxz, g_z_y_z_xxzz, g_z_y_z_xyy, g_z_y_z_xyyz, g_z_y_z_xyz, g_z_y_z_xyzz, g_z_y_z_xzz, g_z_y_z_xzzz, g_z_y_z_yyy, g_z_y_z_yyyz, g_z_y_z_yyz, g_z_y_z_yyzz, g_z_y_z_yzz, g_z_y_z_yzzz, g_z_y_z_zzz, g_z_y_z_zzzz, g_z_y_zz_xxx, g_z_y_zz_xxy, g_z_y_zz_xxz, g_z_y_zz_xyy, g_z_y_zz_xyz, g_z_y_zz_xzz, g_z_y_zz_yyy, g_z_y_zz_yyz, g_z_y_zz_yzz, g_z_y_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zz_xxx[k] = -g_0_y_z_xxx[k] - g_z_y_z_xxx[k] * ab_z + g_z_y_z_xxxz[k];

                g_z_y_zz_xxy[k] = -g_0_y_z_xxy[k] - g_z_y_z_xxy[k] * ab_z + g_z_y_z_xxyz[k];

                g_z_y_zz_xxz[k] = -g_0_y_z_xxz[k] - g_z_y_z_xxz[k] * ab_z + g_z_y_z_xxzz[k];

                g_z_y_zz_xyy[k] = -g_0_y_z_xyy[k] - g_z_y_z_xyy[k] * ab_z + g_z_y_z_xyyz[k];

                g_z_y_zz_xyz[k] = -g_0_y_z_xyz[k] - g_z_y_z_xyz[k] * ab_z + g_z_y_z_xyzz[k];

                g_z_y_zz_xzz[k] = -g_0_y_z_xzz[k] - g_z_y_z_xzz[k] * ab_z + g_z_y_z_xzzz[k];

                g_z_y_zz_yyy[k] = -g_0_y_z_yyy[k] - g_z_y_z_yyy[k] * ab_z + g_z_y_z_yyyz[k];

                g_z_y_zz_yyz[k] = -g_0_y_z_yyz[k] - g_z_y_z_yyz[k] * ab_z + g_z_y_z_yyzz[k];

                g_z_y_zz_yzz[k] = -g_0_y_z_yzz[k] - g_z_y_z_yzz[k] * ab_z + g_z_y_z_yzzz[k];

                g_z_y_zz_zzz[k] = -g_0_y_z_zzz[k] - g_z_y_z_zzz[k] * ab_z + g_z_y_z_zzzz[k];
            }

            /// Set up 480-490 components of targeted buffer : cbuffer.data(

            auto g_z_z_xx_xxx = cbuffer.data(df_geom_11_off + 480 * ccomps * dcomps);

            auto g_z_z_xx_xxy = cbuffer.data(df_geom_11_off + 481 * ccomps * dcomps);

            auto g_z_z_xx_xxz = cbuffer.data(df_geom_11_off + 482 * ccomps * dcomps);

            auto g_z_z_xx_xyy = cbuffer.data(df_geom_11_off + 483 * ccomps * dcomps);

            auto g_z_z_xx_xyz = cbuffer.data(df_geom_11_off + 484 * ccomps * dcomps);

            auto g_z_z_xx_xzz = cbuffer.data(df_geom_11_off + 485 * ccomps * dcomps);

            auto g_z_z_xx_yyy = cbuffer.data(df_geom_11_off + 486 * ccomps * dcomps);

            auto g_z_z_xx_yyz = cbuffer.data(df_geom_11_off + 487 * ccomps * dcomps);

            auto g_z_z_xx_yzz = cbuffer.data(df_geom_11_off + 488 * ccomps * dcomps);

            auto g_z_z_xx_zzz = cbuffer.data(df_geom_11_off + 489 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_x_xxx, g_z_z_x_xxxx, g_z_z_x_xxxy, g_z_z_x_xxxz, g_z_z_x_xxy, g_z_z_x_xxyy, g_z_z_x_xxyz, g_z_z_x_xxz, g_z_z_x_xxzz, g_z_z_x_xyy, g_z_z_x_xyyy, g_z_z_x_xyyz, g_z_z_x_xyz, g_z_z_x_xyzz, g_z_z_x_xzz, g_z_z_x_xzzz, g_z_z_x_yyy, g_z_z_x_yyz, g_z_z_x_yzz, g_z_z_x_zzz, g_z_z_xx_xxx, g_z_z_xx_xxy, g_z_z_xx_xxz, g_z_z_xx_xyy, g_z_z_xx_xyz, g_z_z_xx_xzz, g_z_z_xx_yyy, g_z_z_xx_yyz, g_z_z_xx_yzz, g_z_z_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xx_xxx[k] = -g_z_z_x_xxx[k] * ab_x + g_z_z_x_xxxx[k];

                g_z_z_xx_xxy[k] = -g_z_z_x_xxy[k] * ab_x + g_z_z_x_xxxy[k];

                g_z_z_xx_xxz[k] = -g_z_z_x_xxz[k] * ab_x + g_z_z_x_xxxz[k];

                g_z_z_xx_xyy[k] = -g_z_z_x_xyy[k] * ab_x + g_z_z_x_xxyy[k];

                g_z_z_xx_xyz[k] = -g_z_z_x_xyz[k] * ab_x + g_z_z_x_xxyz[k];

                g_z_z_xx_xzz[k] = -g_z_z_x_xzz[k] * ab_x + g_z_z_x_xxzz[k];

                g_z_z_xx_yyy[k] = -g_z_z_x_yyy[k] * ab_x + g_z_z_x_xyyy[k];

                g_z_z_xx_yyz[k] = -g_z_z_x_yyz[k] * ab_x + g_z_z_x_xyyz[k];

                g_z_z_xx_yzz[k] = -g_z_z_x_yzz[k] * ab_x + g_z_z_x_xyzz[k];

                g_z_z_xx_zzz[k] = -g_z_z_x_zzz[k] * ab_x + g_z_z_x_xzzz[k];
            }

            /// Set up 490-500 components of targeted buffer : cbuffer.data(

            auto g_z_z_xy_xxx = cbuffer.data(df_geom_11_off + 490 * ccomps * dcomps);

            auto g_z_z_xy_xxy = cbuffer.data(df_geom_11_off + 491 * ccomps * dcomps);

            auto g_z_z_xy_xxz = cbuffer.data(df_geom_11_off + 492 * ccomps * dcomps);

            auto g_z_z_xy_xyy = cbuffer.data(df_geom_11_off + 493 * ccomps * dcomps);

            auto g_z_z_xy_xyz = cbuffer.data(df_geom_11_off + 494 * ccomps * dcomps);

            auto g_z_z_xy_xzz = cbuffer.data(df_geom_11_off + 495 * ccomps * dcomps);

            auto g_z_z_xy_yyy = cbuffer.data(df_geom_11_off + 496 * ccomps * dcomps);

            auto g_z_z_xy_yyz = cbuffer.data(df_geom_11_off + 497 * ccomps * dcomps);

            auto g_z_z_xy_yzz = cbuffer.data(df_geom_11_off + 498 * ccomps * dcomps);

            auto g_z_z_xy_zzz = cbuffer.data(df_geom_11_off + 499 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xy_xxx, g_z_z_xy_xxy, g_z_z_xy_xxz, g_z_z_xy_xyy, g_z_z_xy_xyz, g_z_z_xy_xzz, g_z_z_xy_yyy, g_z_z_xy_yyz, g_z_z_xy_yzz, g_z_z_xy_zzz, g_z_z_y_xxx, g_z_z_y_xxxx, g_z_z_y_xxxy, g_z_z_y_xxxz, g_z_z_y_xxy, g_z_z_y_xxyy, g_z_z_y_xxyz, g_z_z_y_xxz, g_z_z_y_xxzz, g_z_z_y_xyy, g_z_z_y_xyyy, g_z_z_y_xyyz, g_z_z_y_xyz, g_z_z_y_xyzz, g_z_z_y_xzz, g_z_z_y_xzzz, g_z_z_y_yyy, g_z_z_y_yyz, g_z_z_y_yzz, g_z_z_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xy_xxx[k] = -g_z_z_y_xxx[k] * ab_x + g_z_z_y_xxxx[k];

                g_z_z_xy_xxy[k] = -g_z_z_y_xxy[k] * ab_x + g_z_z_y_xxxy[k];

                g_z_z_xy_xxz[k] = -g_z_z_y_xxz[k] * ab_x + g_z_z_y_xxxz[k];

                g_z_z_xy_xyy[k] = -g_z_z_y_xyy[k] * ab_x + g_z_z_y_xxyy[k];

                g_z_z_xy_xyz[k] = -g_z_z_y_xyz[k] * ab_x + g_z_z_y_xxyz[k];

                g_z_z_xy_xzz[k] = -g_z_z_y_xzz[k] * ab_x + g_z_z_y_xxzz[k];

                g_z_z_xy_yyy[k] = -g_z_z_y_yyy[k] * ab_x + g_z_z_y_xyyy[k];

                g_z_z_xy_yyz[k] = -g_z_z_y_yyz[k] * ab_x + g_z_z_y_xyyz[k];

                g_z_z_xy_yzz[k] = -g_z_z_y_yzz[k] * ab_x + g_z_z_y_xyzz[k];

                g_z_z_xy_zzz[k] = -g_z_z_y_zzz[k] * ab_x + g_z_z_y_xzzz[k];
            }

            /// Set up 500-510 components of targeted buffer : cbuffer.data(

            auto g_z_z_xz_xxx = cbuffer.data(df_geom_11_off + 500 * ccomps * dcomps);

            auto g_z_z_xz_xxy = cbuffer.data(df_geom_11_off + 501 * ccomps * dcomps);

            auto g_z_z_xz_xxz = cbuffer.data(df_geom_11_off + 502 * ccomps * dcomps);

            auto g_z_z_xz_xyy = cbuffer.data(df_geom_11_off + 503 * ccomps * dcomps);

            auto g_z_z_xz_xyz = cbuffer.data(df_geom_11_off + 504 * ccomps * dcomps);

            auto g_z_z_xz_xzz = cbuffer.data(df_geom_11_off + 505 * ccomps * dcomps);

            auto g_z_z_xz_yyy = cbuffer.data(df_geom_11_off + 506 * ccomps * dcomps);

            auto g_z_z_xz_yyz = cbuffer.data(df_geom_11_off + 507 * ccomps * dcomps);

            auto g_z_z_xz_yzz = cbuffer.data(df_geom_11_off + 508 * ccomps * dcomps);

            auto g_z_z_xz_zzz = cbuffer.data(df_geom_11_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xz_xxx, g_z_z_xz_xxy, g_z_z_xz_xxz, g_z_z_xz_xyy, g_z_z_xz_xyz, g_z_z_xz_xzz, g_z_z_xz_yyy, g_z_z_xz_yyz, g_z_z_xz_yzz, g_z_z_xz_zzz, g_z_z_z_xxx, g_z_z_z_xxxx, g_z_z_z_xxxy, g_z_z_z_xxxz, g_z_z_z_xxy, g_z_z_z_xxyy, g_z_z_z_xxyz, g_z_z_z_xxz, g_z_z_z_xxzz, g_z_z_z_xyy, g_z_z_z_xyyy, g_z_z_z_xyyz, g_z_z_z_xyz, g_z_z_z_xyzz, g_z_z_z_xzz, g_z_z_z_xzzz, g_z_z_z_yyy, g_z_z_z_yyz, g_z_z_z_yzz, g_z_z_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xz_xxx[k] = -g_z_z_z_xxx[k] * ab_x + g_z_z_z_xxxx[k];

                g_z_z_xz_xxy[k] = -g_z_z_z_xxy[k] * ab_x + g_z_z_z_xxxy[k];

                g_z_z_xz_xxz[k] = -g_z_z_z_xxz[k] * ab_x + g_z_z_z_xxxz[k];

                g_z_z_xz_xyy[k] = -g_z_z_z_xyy[k] * ab_x + g_z_z_z_xxyy[k];

                g_z_z_xz_xyz[k] = -g_z_z_z_xyz[k] * ab_x + g_z_z_z_xxyz[k];

                g_z_z_xz_xzz[k] = -g_z_z_z_xzz[k] * ab_x + g_z_z_z_xxzz[k];

                g_z_z_xz_yyy[k] = -g_z_z_z_yyy[k] * ab_x + g_z_z_z_xyyy[k];

                g_z_z_xz_yyz[k] = -g_z_z_z_yyz[k] * ab_x + g_z_z_z_xyyz[k];

                g_z_z_xz_yzz[k] = -g_z_z_z_yzz[k] * ab_x + g_z_z_z_xyzz[k];

                g_z_z_xz_zzz[k] = -g_z_z_z_zzz[k] * ab_x + g_z_z_z_xzzz[k];
            }

            /// Set up 510-520 components of targeted buffer : cbuffer.data(

            auto g_z_z_yy_xxx = cbuffer.data(df_geom_11_off + 510 * ccomps * dcomps);

            auto g_z_z_yy_xxy = cbuffer.data(df_geom_11_off + 511 * ccomps * dcomps);

            auto g_z_z_yy_xxz = cbuffer.data(df_geom_11_off + 512 * ccomps * dcomps);

            auto g_z_z_yy_xyy = cbuffer.data(df_geom_11_off + 513 * ccomps * dcomps);

            auto g_z_z_yy_xyz = cbuffer.data(df_geom_11_off + 514 * ccomps * dcomps);

            auto g_z_z_yy_xzz = cbuffer.data(df_geom_11_off + 515 * ccomps * dcomps);

            auto g_z_z_yy_yyy = cbuffer.data(df_geom_11_off + 516 * ccomps * dcomps);

            auto g_z_z_yy_yyz = cbuffer.data(df_geom_11_off + 517 * ccomps * dcomps);

            auto g_z_z_yy_yzz = cbuffer.data(df_geom_11_off + 518 * ccomps * dcomps);

            auto g_z_z_yy_zzz = cbuffer.data(df_geom_11_off + 519 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_y_xxx, g_z_z_y_xxxy, g_z_z_y_xxy, g_z_z_y_xxyy, g_z_z_y_xxyz, g_z_z_y_xxz, g_z_z_y_xyy, g_z_z_y_xyyy, g_z_z_y_xyyz, g_z_z_y_xyz, g_z_z_y_xyzz, g_z_z_y_xzz, g_z_z_y_yyy, g_z_z_y_yyyy, g_z_z_y_yyyz, g_z_z_y_yyz, g_z_z_y_yyzz, g_z_z_y_yzz, g_z_z_y_yzzz, g_z_z_y_zzz, g_z_z_yy_xxx, g_z_z_yy_xxy, g_z_z_yy_xxz, g_z_z_yy_xyy, g_z_z_yy_xyz, g_z_z_yy_xzz, g_z_z_yy_yyy, g_z_z_yy_yyz, g_z_z_yy_yzz, g_z_z_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yy_xxx[k] = -g_z_z_y_xxx[k] * ab_y + g_z_z_y_xxxy[k];

                g_z_z_yy_xxy[k] = -g_z_z_y_xxy[k] * ab_y + g_z_z_y_xxyy[k];

                g_z_z_yy_xxz[k] = -g_z_z_y_xxz[k] * ab_y + g_z_z_y_xxyz[k];

                g_z_z_yy_xyy[k] = -g_z_z_y_xyy[k] * ab_y + g_z_z_y_xyyy[k];

                g_z_z_yy_xyz[k] = -g_z_z_y_xyz[k] * ab_y + g_z_z_y_xyyz[k];

                g_z_z_yy_xzz[k] = -g_z_z_y_xzz[k] * ab_y + g_z_z_y_xyzz[k];

                g_z_z_yy_yyy[k] = -g_z_z_y_yyy[k] * ab_y + g_z_z_y_yyyy[k];

                g_z_z_yy_yyz[k] = -g_z_z_y_yyz[k] * ab_y + g_z_z_y_yyyz[k];

                g_z_z_yy_yzz[k] = -g_z_z_y_yzz[k] * ab_y + g_z_z_y_yyzz[k];

                g_z_z_yy_zzz[k] = -g_z_z_y_zzz[k] * ab_y + g_z_z_y_yzzz[k];
            }

            /// Set up 520-530 components of targeted buffer : cbuffer.data(

            auto g_z_z_yz_xxx = cbuffer.data(df_geom_11_off + 520 * ccomps * dcomps);

            auto g_z_z_yz_xxy = cbuffer.data(df_geom_11_off + 521 * ccomps * dcomps);

            auto g_z_z_yz_xxz = cbuffer.data(df_geom_11_off + 522 * ccomps * dcomps);

            auto g_z_z_yz_xyy = cbuffer.data(df_geom_11_off + 523 * ccomps * dcomps);

            auto g_z_z_yz_xyz = cbuffer.data(df_geom_11_off + 524 * ccomps * dcomps);

            auto g_z_z_yz_xzz = cbuffer.data(df_geom_11_off + 525 * ccomps * dcomps);

            auto g_z_z_yz_yyy = cbuffer.data(df_geom_11_off + 526 * ccomps * dcomps);

            auto g_z_z_yz_yyz = cbuffer.data(df_geom_11_off + 527 * ccomps * dcomps);

            auto g_z_z_yz_yzz = cbuffer.data(df_geom_11_off + 528 * ccomps * dcomps);

            auto g_z_z_yz_zzz = cbuffer.data(df_geom_11_off + 529 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yz_xxx, g_z_z_yz_xxy, g_z_z_yz_xxz, g_z_z_yz_xyy, g_z_z_yz_xyz, g_z_z_yz_xzz, g_z_z_yz_yyy, g_z_z_yz_yyz, g_z_z_yz_yzz, g_z_z_yz_zzz, g_z_z_z_xxx, g_z_z_z_xxxy, g_z_z_z_xxy, g_z_z_z_xxyy, g_z_z_z_xxyz, g_z_z_z_xxz, g_z_z_z_xyy, g_z_z_z_xyyy, g_z_z_z_xyyz, g_z_z_z_xyz, g_z_z_z_xyzz, g_z_z_z_xzz, g_z_z_z_yyy, g_z_z_z_yyyy, g_z_z_z_yyyz, g_z_z_z_yyz, g_z_z_z_yyzz, g_z_z_z_yzz, g_z_z_z_yzzz, g_z_z_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yz_xxx[k] = -g_z_z_z_xxx[k] * ab_y + g_z_z_z_xxxy[k];

                g_z_z_yz_xxy[k] = -g_z_z_z_xxy[k] * ab_y + g_z_z_z_xxyy[k];

                g_z_z_yz_xxz[k] = -g_z_z_z_xxz[k] * ab_y + g_z_z_z_xxyz[k];

                g_z_z_yz_xyy[k] = -g_z_z_z_xyy[k] * ab_y + g_z_z_z_xyyy[k];

                g_z_z_yz_xyz[k] = -g_z_z_z_xyz[k] * ab_y + g_z_z_z_xyyz[k];

                g_z_z_yz_xzz[k] = -g_z_z_z_xzz[k] * ab_y + g_z_z_z_xyzz[k];

                g_z_z_yz_yyy[k] = -g_z_z_z_yyy[k] * ab_y + g_z_z_z_yyyy[k];

                g_z_z_yz_yyz[k] = -g_z_z_z_yyz[k] * ab_y + g_z_z_z_yyyz[k];

                g_z_z_yz_yzz[k] = -g_z_z_z_yzz[k] * ab_y + g_z_z_z_yyzz[k];

                g_z_z_yz_zzz[k] = -g_z_z_z_zzz[k] * ab_y + g_z_z_z_yzzz[k];
            }

            /// Set up 530-540 components of targeted buffer : cbuffer.data(

            auto g_z_z_zz_xxx = cbuffer.data(df_geom_11_off + 530 * ccomps * dcomps);

            auto g_z_z_zz_xxy = cbuffer.data(df_geom_11_off + 531 * ccomps * dcomps);

            auto g_z_z_zz_xxz = cbuffer.data(df_geom_11_off + 532 * ccomps * dcomps);

            auto g_z_z_zz_xyy = cbuffer.data(df_geom_11_off + 533 * ccomps * dcomps);

            auto g_z_z_zz_xyz = cbuffer.data(df_geom_11_off + 534 * ccomps * dcomps);

            auto g_z_z_zz_xzz = cbuffer.data(df_geom_11_off + 535 * ccomps * dcomps);

            auto g_z_z_zz_yyy = cbuffer.data(df_geom_11_off + 536 * ccomps * dcomps);

            auto g_z_z_zz_yyz = cbuffer.data(df_geom_11_off + 537 * ccomps * dcomps);

            auto g_z_z_zz_yzz = cbuffer.data(df_geom_11_off + 538 * ccomps * dcomps);

            auto g_z_z_zz_zzz = cbuffer.data(df_geom_11_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxx, g_0_z_z_xxy, g_0_z_z_xxz, g_0_z_z_xyy, g_0_z_z_xyz, g_0_z_z_xzz, g_0_z_z_yyy, g_0_z_z_yyz, g_0_z_z_yzz, g_0_z_z_zzz, g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz, g_z_z_z_xxx, g_z_z_z_xxxz, g_z_z_z_xxy, g_z_z_z_xxyz, g_z_z_z_xxz, g_z_z_z_xxzz, g_z_z_z_xyy, g_z_z_z_xyyz, g_z_z_z_xyz, g_z_z_z_xyzz, g_z_z_z_xzz, g_z_z_z_xzzz, g_z_z_z_yyy, g_z_z_z_yyyz, g_z_z_z_yyz, g_z_z_z_yyzz, g_z_z_z_yzz, g_z_z_z_yzzz, g_z_z_z_zzz, g_z_z_z_zzzz, g_z_z_zz_xxx, g_z_z_zz_xxy, g_z_z_zz_xxz, g_z_z_zz_xyy, g_z_z_zz_xyz, g_z_z_zz_xzz, g_z_z_zz_yyy, g_z_z_zz_yyz, g_z_z_zz_yzz, g_z_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zz_xxx[k] = -g_0_z_z_xxx[k] + g_z_0_z_xxx[k] - g_z_z_z_xxx[k] * ab_z + g_z_z_z_xxxz[k];

                g_z_z_zz_xxy[k] = -g_0_z_z_xxy[k] + g_z_0_z_xxy[k] - g_z_z_z_xxy[k] * ab_z + g_z_z_z_xxyz[k];

                g_z_z_zz_xxz[k] = -g_0_z_z_xxz[k] + g_z_0_z_xxz[k] - g_z_z_z_xxz[k] * ab_z + g_z_z_z_xxzz[k];

                g_z_z_zz_xyy[k] = -g_0_z_z_xyy[k] + g_z_0_z_xyy[k] - g_z_z_z_xyy[k] * ab_z + g_z_z_z_xyyz[k];

                g_z_z_zz_xyz[k] = -g_0_z_z_xyz[k] + g_z_0_z_xyz[k] - g_z_z_z_xyz[k] * ab_z + g_z_z_z_xyzz[k];

                g_z_z_zz_xzz[k] = -g_0_z_z_xzz[k] + g_z_0_z_xzz[k] - g_z_z_z_xzz[k] * ab_z + g_z_z_z_xzzz[k];

                g_z_z_zz_yyy[k] = -g_0_z_z_yyy[k] + g_z_0_z_yyy[k] - g_z_z_z_yyy[k] * ab_z + g_z_z_z_yyyz[k];

                g_z_z_zz_yyz[k] = -g_0_z_z_yyz[k] + g_z_0_z_yyz[k] - g_z_z_z_yyz[k] * ab_z + g_z_z_z_yyzz[k];

                g_z_z_zz_yzz[k] = -g_0_z_z_yzz[k] + g_z_0_z_yzz[k] - g_z_z_z_yzz[k] * ab_z + g_z_z_z_yzzz[k];

                g_z_z_zz_zzz[k] = -g_0_z_z_zzz[k] + g_z_0_z_zzz[k] - g_z_z_z_zzz[k] * ab_z + g_z_z_z_zzzz[k];
            }
        }
    }
}

} // erirec namespace

