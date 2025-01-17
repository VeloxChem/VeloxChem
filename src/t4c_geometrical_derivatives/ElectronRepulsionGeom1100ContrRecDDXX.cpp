#include "ElectronRepulsionGeom1100ContrRecDDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_ddxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_ddxx,
                                            const size_t idx_geom_01_pdxx,
                                            const size_t idx_geom_10_pdxx,
                                            const size_t idx_geom_11_pdxx,
                                            const size_t idx_geom_11_pfxx,
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
            /// Set up components of auxilary buffer : PDSS

            const auto pd_geom_01_off = idx_geom_01_pdxx + i * dcomps + j;

            auto g_0_x_x_xx = cbuffer.data(pd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xy = cbuffer.data(pd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xz = cbuffer.data(pd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_yy = cbuffer.data(pd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_yz = cbuffer.data(pd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_zz = cbuffer.data(pd_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_y_xx = cbuffer.data(pd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_y_xy = cbuffer.data(pd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_y_xz = cbuffer.data(pd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_y_yy = cbuffer.data(pd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_y_yz = cbuffer.data(pd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_y_zz = cbuffer.data(pd_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_z_xx = cbuffer.data(pd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_z_xy = cbuffer.data(pd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_z_xz = cbuffer.data(pd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_z_yy = cbuffer.data(pd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_z_yz = cbuffer.data(pd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_z_zz = cbuffer.data(pd_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_x_xx = cbuffer.data(pd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_x_xy = cbuffer.data(pd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_y_x_xz = cbuffer.data(pd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_y_x_yy = cbuffer.data(pd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_y_x_yz = cbuffer.data(pd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_y_x_zz = cbuffer.data(pd_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_y_y_xx = cbuffer.data(pd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_y_y_xy = cbuffer.data(pd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_y_y_xz = cbuffer.data(pd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_y_y_yy = cbuffer.data(pd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_y_yz = cbuffer.data(pd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_y_zz = cbuffer.data(pd_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_z_xx = cbuffer.data(pd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_z_xy = cbuffer.data(pd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_z_xz = cbuffer.data(pd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_z_yy = cbuffer.data(pd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_z_yz = cbuffer.data(pd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_z_zz = cbuffer.data(pd_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_z_x_xx = cbuffer.data(pd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_z_x_xy = cbuffer.data(pd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_z_x_xz = cbuffer.data(pd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_z_x_yy = cbuffer.data(pd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_z_x_yz = cbuffer.data(pd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_z_x_zz = cbuffer.data(pd_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_z_y_xx = cbuffer.data(pd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_z_y_xy = cbuffer.data(pd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_z_y_xz = cbuffer.data(pd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_z_y_yy = cbuffer.data(pd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_z_y_yz = cbuffer.data(pd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_z_y_zz = cbuffer.data(pd_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_z_z_xx = cbuffer.data(pd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_z_z_xy = cbuffer.data(pd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_z_z_xz = cbuffer.data(pd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_z_z_yy = cbuffer.data(pd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_z_z_yz = cbuffer.data(pd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_z_z_zz = cbuffer.data(pd_geom_01_off + 53 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PDSS

            const auto pd_geom_10_off = idx_geom_10_pdxx + i * dcomps + j;

            auto g_x_0_x_xx = cbuffer.data(pd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xy = cbuffer.data(pd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xz = cbuffer.data(pd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_yy = cbuffer.data(pd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_yz = cbuffer.data(pd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_zz = cbuffer.data(pd_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_y_xx = cbuffer.data(pd_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_y_xy = cbuffer.data(pd_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_y_xz = cbuffer.data(pd_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_y_yy = cbuffer.data(pd_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_y_yz = cbuffer.data(pd_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_y_zz = cbuffer.data(pd_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_z_xx = cbuffer.data(pd_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_z_xy = cbuffer.data(pd_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_z_xz = cbuffer.data(pd_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_z_yy = cbuffer.data(pd_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_z_yz = cbuffer.data(pd_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_z_zz = cbuffer.data(pd_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_x_xx = cbuffer.data(pd_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_x_xy = cbuffer.data(pd_geom_10_off + 19 * ccomps * dcomps);

            auto g_y_0_x_xz = cbuffer.data(pd_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_x_yy = cbuffer.data(pd_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_x_yz = cbuffer.data(pd_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_x_zz = cbuffer.data(pd_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_y_xx = cbuffer.data(pd_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_y_xy = cbuffer.data(pd_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_y_xz = cbuffer.data(pd_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_y_yy = cbuffer.data(pd_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_y_yz = cbuffer.data(pd_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_y_zz = cbuffer.data(pd_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_z_xx = cbuffer.data(pd_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_z_xy = cbuffer.data(pd_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_z_xz = cbuffer.data(pd_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_z_yy = cbuffer.data(pd_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_z_yz = cbuffer.data(pd_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_z_zz = cbuffer.data(pd_geom_10_off + 35 * ccomps * dcomps);

            auto g_z_0_x_xx = cbuffer.data(pd_geom_10_off + 36 * ccomps * dcomps);

            auto g_z_0_x_xy = cbuffer.data(pd_geom_10_off + 37 * ccomps * dcomps);

            auto g_z_0_x_xz = cbuffer.data(pd_geom_10_off + 38 * ccomps * dcomps);

            auto g_z_0_x_yy = cbuffer.data(pd_geom_10_off + 39 * ccomps * dcomps);

            auto g_z_0_x_yz = cbuffer.data(pd_geom_10_off + 40 * ccomps * dcomps);

            auto g_z_0_x_zz = cbuffer.data(pd_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_y_xx = cbuffer.data(pd_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_y_xy = cbuffer.data(pd_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_y_xz = cbuffer.data(pd_geom_10_off + 44 * ccomps * dcomps);

            auto g_z_0_y_yy = cbuffer.data(pd_geom_10_off + 45 * ccomps * dcomps);

            auto g_z_0_y_yz = cbuffer.data(pd_geom_10_off + 46 * ccomps * dcomps);

            auto g_z_0_y_zz = cbuffer.data(pd_geom_10_off + 47 * ccomps * dcomps);

            auto g_z_0_z_xx = cbuffer.data(pd_geom_10_off + 48 * ccomps * dcomps);

            auto g_z_0_z_xy = cbuffer.data(pd_geom_10_off + 49 * ccomps * dcomps);

            auto g_z_0_z_xz = cbuffer.data(pd_geom_10_off + 50 * ccomps * dcomps);

            auto g_z_0_z_yy = cbuffer.data(pd_geom_10_off + 51 * ccomps * dcomps);

            auto g_z_0_z_yz = cbuffer.data(pd_geom_10_off + 52 * ccomps * dcomps);

            auto g_z_0_z_zz = cbuffer.data(pd_geom_10_off + 53 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PDSS

            const auto pd_geom_11_off = idx_geom_11_pdxx + i * dcomps + j;

            auto g_x_x_x_xx = cbuffer.data(pd_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xy = cbuffer.data(pd_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xz = cbuffer.data(pd_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_yy = cbuffer.data(pd_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_yz = cbuffer.data(pd_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_zz = cbuffer.data(pd_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_y_xx = cbuffer.data(pd_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_y_xy = cbuffer.data(pd_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_y_xz = cbuffer.data(pd_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_y_yy = cbuffer.data(pd_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_y_yz = cbuffer.data(pd_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_y_zz = cbuffer.data(pd_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_z_xx = cbuffer.data(pd_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_z_xy = cbuffer.data(pd_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_z_xz = cbuffer.data(pd_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_z_yy = cbuffer.data(pd_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_z_yz = cbuffer.data(pd_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_z_zz = cbuffer.data(pd_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_y_x_xx = cbuffer.data(pd_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_y_x_xy = cbuffer.data(pd_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_y_x_xz = cbuffer.data(pd_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_y_x_yy = cbuffer.data(pd_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_y_x_yz = cbuffer.data(pd_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_y_x_zz = cbuffer.data(pd_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_y_y_xx = cbuffer.data(pd_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_y_y_xy = cbuffer.data(pd_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_y_y_xz = cbuffer.data(pd_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_y_y_yy = cbuffer.data(pd_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_y_y_yz = cbuffer.data(pd_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_y_zz = cbuffer.data(pd_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_y_z_xx = cbuffer.data(pd_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_z_xy = cbuffer.data(pd_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_z_xz = cbuffer.data(pd_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_y_z_yy = cbuffer.data(pd_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_z_yz = cbuffer.data(pd_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_z_zz = cbuffer.data(pd_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_z_x_xx = cbuffer.data(pd_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_z_x_xy = cbuffer.data(pd_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_z_x_xz = cbuffer.data(pd_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_z_x_yy = cbuffer.data(pd_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_z_x_yz = cbuffer.data(pd_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_z_x_zz = cbuffer.data(pd_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_z_y_xx = cbuffer.data(pd_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_z_y_xy = cbuffer.data(pd_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_z_y_xz = cbuffer.data(pd_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_z_y_yy = cbuffer.data(pd_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_z_y_yz = cbuffer.data(pd_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_z_y_zz = cbuffer.data(pd_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_z_z_xx = cbuffer.data(pd_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_z_z_xy = cbuffer.data(pd_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_z_z_xz = cbuffer.data(pd_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_z_z_yy = cbuffer.data(pd_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_z_z_yz = cbuffer.data(pd_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_z_z_zz = cbuffer.data(pd_geom_11_off + 53 * ccomps * dcomps);

            auto g_y_x_x_xx = cbuffer.data(pd_geom_11_off + 54 * ccomps * dcomps);

            auto g_y_x_x_xy = cbuffer.data(pd_geom_11_off + 55 * ccomps * dcomps);

            auto g_y_x_x_xz = cbuffer.data(pd_geom_11_off + 56 * ccomps * dcomps);

            auto g_y_x_x_yy = cbuffer.data(pd_geom_11_off + 57 * ccomps * dcomps);

            auto g_y_x_x_yz = cbuffer.data(pd_geom_11_off + 58 * ccomps * dcomps);

            auto g_y_x_x_zz = cbuffer.data(pd_geom_11_off + 59 * ccomps * dcomps);

            auto g_y_x_y_xx = cbuffer.data(pd_geom_11_off + 60 * ccomps * dcomps);

            auto g_y_x_y_xy = cbuffer.data(pd_geom_11_off + 61 * ccomps * dcomps);

            auto g_y_x_y_xz = cbuffer.data(pd_geom_11_off + 62 * ccomps * dcomps);

            auto g_y_x_y_yy = cbuffer.data(pd_geom_11_off + 63 * ccomps * dcomps);

            auto g_y_x_y_yz = cbuffer.data(pd_geom_11_off + 64 * ccomps * dcomps);

            auto g_y_x_y_zz = cbuffer.data(pd_geom_11_off + 65 * ccomps * dcomps);

            auto g_y_x_z_xx = cbuffer.data(pd_geom_11_off + 66 * ccomps * dcomps);

            auto g_y_x_z_xy = cbuffer.data(pd_geom_11_off + 67 * ccomps * dcomps);

            auto g_y_x_z_xz = cbuffer.data(pd_geom_11_off + 68 * ccomps * dcomps);

            auto g_y_x_z_yy = cbuffer.data(pd_geom_11_off + 69 * ccomps * dcomps);

            auto g_y_x_z_yz = cbuffer.data(pd_geom_11_off + 70 * ccomps * dcomps);

            auto g_y_x_z_zz = cbuffer.data(pd_geom_11_off + 71 * ccomps * dcomps);

            auto g_y_y_x_xx = cbuffer.data(pd_geom_11_off + 72 * ccomps * dcomps);

            auto g_y_y_x_xy = cbuffer.data(pd_geom_11_off + 73 * ccomps * dcomps);

            auto g_y_y_x_xz = cbuffer.data(pd_geom_11_off + 74 * ccomps * dcomps);

            auto g_y_y_x_yy = cbuffer.data(pd_geom_11_off + 75 * ccomps * dcomps);

            auto g_y_y_x_yz = cbuffer.data(pd_geom_11_off + 76 * ccomps * dcomps);

            auto g_y_y_x_zz = cbuffer.data(pd_geom_11_off + 77 * ccomps * dcomps);

            auto g_y_y_y_xx = cbuffer.data(pd_geom_11_off + 78 * ccomps * dcomps);

            auto g_y_y_y_xy = cbuffer.data(pd_geom_11_off + 79 * ccomps * dcomps);

            auto g_y_y_y_xz = cbuffer.data(pd_geom_11_off + 80 * ccomps * dcomps);

            auto g_y_y_y_yy = cbuffer.data(pd_geom_11_off + 81 * ccomps * dcomps);

            auto g_y_y_y_yz = cbuffer.data(pd_geom_11_off + 82 * ccomps * dcomps);

            auto g_y_y_y_zz = cbuffer.data(pd_geom_11_off + 83 * ccomps * dcomps);

            auto g_y_y_z_xx = cbuffer.data(pd_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_y_z_xy = cbuffer.data(pd_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_y_z_xz = cbuffer.data(pd_geom_11_off + 86 * ccomps * dcomps);

            auto g_y_y_z_yy = cbuffer.data(pd_geom_11_off + 87 * ccomps * dcomps);

            auto g_y_y_z_yz = cbuffer.data(pd_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_y_z_zz = cbuffer.data(pd_geom_11_off + 89 * ccomps * dcomps);

            auto g_y_z_x_xx = cbuffer.data(pd_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_z_x_xy = cbuffer.data(pd_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_z_x_xz = cbuffer.data(pd_geom_11_off + 92 * ccomps * dcomps);

            auto g_y_z_x_yy = cbuffer.data(pd_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_z_x_yz = cbuffer.data(pd_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_z_x_zz = cbuffer.data(pd_geom_11_off + 95 * ccomps * dcomps);

            auto g_y_z_y_xx = cbuffer.data(pd_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_z_y_xy = cbuffer.data(pd_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_z_y_xz = cbuffer.data(pd_geom_11_off + 98 * ccomps * dcomps);

            auto g_y_z_y_yy = cbuffer.data(pd_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_z_y_yz = cbuffer.data(pd_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_z_y_zz = cbuffer.data(pd_geom_11_off + 101 * ccomps * dcomps);

            auto g_y_z_z_xx = cbuffer.data(pd_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_z_z_xy = cbuffer.data(pd_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_z_z_xz = cbuffer.data(pd_geom_11_off + 104 * ccomps * dcomps);

            auto g_y_z_z_yy = cbuffer.data(pd_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_z_z_yz = cbuffer.data(pd_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_z_z_zz = cbuffer.data(pd_geom_11_off + 107 * ccomps * dcomps);

            auto g_z_x_x_xx = cbuffer.data(pd_geom_11_off + 108 * ccomps * dcomps);

            auto g_z_x_x_xy = cbuffer.data(pd_geom_11_off + 109 * ccomps * dcomps);

            auto g_z_x_x_xz = cbuffer.data(pd_geom_11_off + 110 * ccomps * dcomps);

            auto g_z_x_x_yy = cbuffer.data(pd_geom_11_off + 111 * ccomps * dcomps);

            auto g_z_x_x_yz = cbuffer.data(pd_geom_11_off + 112 * ccomps * dcomps);

            auto g_z_x_x_zz = cbuffer.data(pd_geom_11_off + 113 * ccomps * dcomps);

            auto g_z_x_y_xx = cbuffer.data(pd_geom_11_off + 114 * ccomps * dcomps);

            auto g_z_x_y_xy = cbuffer.data(pd_geom_11_off + 115 * ccomps * dcomps);

            auto g_z_x_y_xz = cbuffer.data(pd_geom_11_off + 116 * ccomps * dcomps);

            auto g_z_x_y_yy = cbuffer.data(pd_geom_11_off + 117 * ccomps * dcomps);

            auto g_z_x_y_yz = cbuffer.data(pd_geom_11_off + 118 * ccomps * dcomps);

            auto g_z_x_y_zz = cbuffer.data(pd_geom_11_off + 119 * ccomps * dcomps);

            auto g_z_x_z_xx = cbuffer.data(pd_geom_11_off + 120 * ccomps * dcomps);

            auto g_z_x_z_xy = cbuffer.data(pd_geom_11_off + 121 * ccomps * dcomps);

            auto g_z_x_z_xz = cbuffer.data(pd_geom_11_off + 122 * ccomps * dcomps);

            auto g_z_x_z_yy = cbuffer.data(pd_geom_11_off + 123 * ccomps * dcomps);

            auto g_z_x_z_yz = cbuffer.data(pd_geom_11_off + 124 * ccomps * dcomps);

            auto g_z_x_z_zz = cbuffer.data(pd_geom_11_off + 125 * ccomps * dcomps);

            auto g_z_y_x_xx = cbuffer.data(pd_geom_11_off + 126 * ccomps * dcomps);

            auto g_z_y_x_xy = cbuffer.data(pd_geom_11_off + 127 * ccomps * dcomps);

            auto g_z_y_x_xz = cbuffer.data(pd_geom_11_off + 128 * ccomps * dcomps);

            auto g_z_y_x_yy = cbuffer.data(pd_geom_11_off + 129 * ccomps * dcomps);

            auto g_z_y_x_yz = cbuffer.data(pd_geom_11_off + 130 * ccomps * dcomps);

            auto g_z_y_x_zz = cbuffer.data(pd_geom_11_off + 131 * ccomps * dcomps);

            auto g_z_y_y_xx = cbuffer.data(pd_geom_11_off + 132 * ccomps * dcomps);

            auto g_z_y_y_xy = cbuffer.data(pd_geom_11_off + 133 * ccomps * dcomps);

            auto g_z_y_y_xz = cbuffer.data(pd_geom_11_off + 134 * ccomps * dcomps);

            auto g_z_y_y_yy = cbuffer.data(pd_geom_11_off + 135 * ccomps * dcomps);

            auto g_z_y_y_yz = cbuffer.data(pd_geom_11_off + 136 * ccomps * dcomps);

            auto g_z_y_y_zz = cbuffer.data(pd_geom_11_off + 137 * ccomps * dcomps);

            auto g_z_y_z_xx = cbuffer.data(pd_geom_11_off + 138 * ccomps * dcomps);

            auto g_z_y_z_xy = cbuffer.data(pd_geom_11_off + 139 * ccomps * dcomps);

            auto g_z_y_z_xz = cbuffer.data(pd_geom_11_off + 140 * ccomps * dcomps);

            auto g_z_y_z_yy = cbuffer.data(pd_geom_11_off + 141 * ccomps * dcomps);

            auto g_z_y_z_yz = cbuffer.data(pd_geom_11_off + 142 * ccomps * dcomps);

            auto g_z_y_z_zz = cbuffer.data(pd_geom_11_off + 143 * ccomps * dcomps);

            auto g_z_z_x_xx = cbuffer.data(pd_geom_11_off + 144 * ccomps * dcomps);

            auto g_z_z_x_xy = cbuffer.data(pd_geom_11_off + 145 * ccomps * dcomps);

            auto g_z_z_x_xz = cbuffer.data(pd_geom_11_off + 146 * ccomps * dcomps);

            auto g_z_z_x_yy = cbuffer.data(pd_geom_11_off + 147 * ccomps * dcomps);

            auto g_z_z_x_yz = cbuffer.data(pd_geom_11_off + 148 * ccomps * dcomps);

            auto g_z_z_x_zz = cbuffer.data(pd_geom_11_off + 149 * ccomps * dcomps);

            auto g_z_z_y_xx = cbuffer.data(pd_geom_11_off + 150 * ccomps * dcomps);

            auto g_z_z_y_xy = cbuffer.data(pd_geom_11_off + 151 * ccomps * dcomps);

            auto g_z_z_y_xz = cbuffer.data(pd_geom_11_off + 152 * ccomps * dcomps);

            auto g_z_z_y_yy = cbuffer.data(pd_geom_11_off + 153 * ccomps * dcomps);

            auto g_z_z_y_yz = cbuffer.data(pd_geom_11_off + 154 * ccomps * dcomps);

            auto g_z_z_y_zz = cbuffer.data(pd_geom_11_off + 155 * ccomps * dcomps);

            auto g_z_z_z_xx = cbuffer.data(pd_geom_11_off + 156 * ccomps * dcomps);

            auto g_z_z_z_xy = cbuffer.data(pd_geom_11_off + 157 * ccomps * dcomps);

            auto g_z_z_z_xz = cbuffer.data(pd_geom_11_off + 158 * ccomps * dcomps);

            auto g_z_z_z_yy = cbuffer.data(pd_geom_11_off + 159 * ccomps * dcomps);

            auto g_z_z_z_yz = cbuffer.data(pd_geom_11_off + 160 * ccomps * dcomps);

            auto g_z_z_z_zz = cbuffer.data(pd_geom_11_off + 161 * ccomps * dcomps);

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

            auto g_x_x_y_xxy = cbuffer.data(pf_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_y_xyy = cbuffer.data(pf_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_y_xyz = cbuffer.data(pf_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_y_yyy = cbuffer.data(pf_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_y_yyz = cbuffer.data(pf_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_y_yzz = cbuffer.data(pf_geom_11_off + 18 * ccomps * dcomps);

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

            auto g_x_y_z_xxz = cbuffer.data(pf_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_z_xyz = cbuffer.data(pf_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_z_xzz = cbuffer.data(pf_geom_11_off + 55 * ccomps * dcomps);

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

            auto g_x_z_y_xxy = cbuffer.data(pf_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_z_y_xyy = cbuffer.data(pf_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_y_xyz = cbuffer.data(pf_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_z_y_yyy = cbuffer.data(pf_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_y_yyz = cbuffer.data(pf_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_z_y_yzz = cbuffer.data(pf_geom_11_off + 78 * ccomps * dcomps);

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

            auto g_y_x_z_xxz = cbuffer.data(pf_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_x_z_xyz = cbuffer.data(pf_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_x_z_xzz = cbuffer.data(pf_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_x_z_yyz = cbuffer.data(pf_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_x_z_yzz = cbuffer.data(pf_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_x_z_zzz = cbuffer.data(pf_geom_11_off + 119 * ccomps * dcomps);

            auto g_y_y_x_xxx = cbuffer.data(pf_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_y_x_xxy = cbuffer.data(pf_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_y_x_xxz = cbuffer.data(pf_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_y_x_xyy = cbuffer.data(pf_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_y_x_xyz = cbuffer.data(pf_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_y_x_xzz = cbuffer.data(pf_geom_11_off + 125 * ccomps * dcomps);

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

            auto g_y_y_z_yyz = cbuffer.data(pf_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_y_z_yzz = cbuffer.data(pf_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_y_z_zzz = cbuffer.data(pf_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_z_x_xxx = cbuffer.data(pf_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_z_x_xxy = cbuffer.data(pf_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_z_x_xxz = cbuffer.data(pf_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_z_x_xyy = cbuffer.data(pf_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_z_x_xyz = cbuffer.data(pf_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_z_x_xzz = cbuffer.data(pf_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_z_y_xxx = cbuffer.data(pf_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_z_y_xxy = cbuffer.data(pf_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_z_y_xxz = cbuffer.data(pf_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_z_y_xyy = cbuffer.data(pf_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_z_y_xyz = cbuffer.data(pf_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_z_y_xzz = cbuffer.data(pf_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_z_y_yyy = cbuffer.data(pf_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_z_y_yyz = cbuffer.data(pf_geom_11_off + 167 * ccomps * dcomps);

            auto g_y_z_y_yzz = cbuffer.data(pf_geom_11_off + 168 * ccomps * dcomps);

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

            auto g_z_x_y_xxy = cbuffer.data(pf_geom_11_off + 191 * ccomps * dcomps);

            auto g_z_x_y_xyy = cbuffer.data(pf_geom_11_off + 193 * ccomps * dcomps);

            auto g_z_x_y_xyz = cbuffer.data(pf_geom_11_off + 194 * ccomps * dcomps);

            auto g_z_x_y_yyy = cbuffer.data(pf_geom_11_off + 196 * ccomps * dcomps);

            auto g_z_x_y_yyz = cbuffer.data(pf_geom_11_off + 197 * ccomps * dcomps);

            auto g_z_x_y_yzz = cbuffer.data(pf_geom_11_off + 198 * ccomps * dcomps);

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

            auto g_z_y_y_xxx = cbuffer.data(pf_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_y_y_xxy = cbuffer.data(pf_geom_11_off + 221 * ccomps * dcomps);

            auto g_z_y_y_xxz = cbuffer.data(pf_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_y_y_xyy = cbuffer.data(pf_geom_11_off + 223 * ccomps * dcomps);

            auto g_z_y_y_xyz = cbuffer.data(pf_geom_11_off + 224 * ccomps * dcomps);

            auto g_z_y_y_xzz = cbuffer.data(pf_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_y_y_yyy = cbuffer.data(pf_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_y_y_yyz = cbuffer.data(pf_geom_11_off + 227 * ccomps * dcomps);

            auto g_z_y_y_yzz = cbuffer.data(pf_geom_11_off + 228 * ccomps * dcomps);

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

            auto g_z_z_y_xxx = cbuffer.data(pf_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_z_y_xxy = cbuffer.data(pf_geom_11_off + 251 * ccomps * dcomps);

            auto g_z_z_y_xxz = cbuffer.data(pf_geom_11_off + 252 * ccomps * dcomps);

            auto g_z_z_y_xyy = cbuffer.data(pf_geom_11_off + 253 * ccomps * dcomps);

            auto g_z_z_y_xyz = cbuffer.data(pf_geom_11_off + 254 * ccomps * dcomps);

            auto g_z_z_y_xzz = cbuffer.data(pf_geom_11_off + 255 * ccomps * dcomps);

            auto g_z_z_y_yyy = cbuffer.data(pf_geom_11_off + 256 * ccomps * dcomps);

            auto g_z_z_y_yyz = cbuffer.data(pf_geom_11_off + 257 * ccomps * dcomps);

            auto g_z_z_y_yzz = cbuffer.data(pf_geom_11_off + 258 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_ddxx

            const auto dd_geom_11_off = idx_geom_11_ddxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_xx_xx = cbuffer.data(dd_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_xy = cbuffer.data(dd_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_xz = cbuffer.data(dd_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xx_yy = cbuffer.data(dd_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xx_yz = cbuffer.data(dd_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xx_zz = cbuffer.data(dd_geom_11_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, g_x_x_x_xx, g_x_x_x_xxx, g_x_x_x_xxy, g_x_x_x_xxz, g_x_x_x_xy, g_x_x_x_xyy, g_x_x_x_xyz, g_x_x_x_xz, g_x_x_x_xzz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xx_xx[k] = -g_0_x_x_xx[k] + g_x_0_x_xx[k] - g_x_x_x_xx[k] * ab_x + g_x_x_x_xxx[k];

                g_x_x_xx_xy[k] = -g_0_x_x_xy[k] + g_x_0_x_xy[k] - g_x_x_x_xy[k] * ab_x + g_x_x_x_xxy[k];

                g_x_x_xx_xz[k] = -g_0_x_x_xz[k] + g_x_0_x_xz[k] - g_x_x_x_xz[k] * ab_x + g_x_x_x_xxz[k];

                g_x_x_xx_yy[k] = -g_0_x_x_yy[k] + g_x_0_x_yy[k] - g_x_x_x_yy[k] * ab_x + g_x_x_x_xyy[k];

                g_x_x_xx_yz[k] = -g_0_x_x_yz[k] + g_x_0_x_yz[k] - g_x_x_x_yz[k] * ab_x + g_x_x_x_xyz[k];

                g_x_x_xx_zz[k] = -g_0_x_x_zz[k] + g_x_0_x_zz[k] - g_x_x_x_zz[k] * ab_x + g_x_x_x_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_x_xy_xx = cbuffer.data(dd_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xy_xy = cbuffer.data(dd_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xy_xz = cbuffer.data(dd_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xy_yy = cbuffer.data(dd_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xy_yz = cbuffer.data(dd_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xy_zz = cbuffer.data(dd_geom_11_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xx, g_x_x_x_xxy, g_x_x_x_xy, g_x_x_x_xyy, g_x_x_x_xyz, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yyy, g_x_x_x_yyz, g_x_x_x_yz, g_x_x_x_yzz, g_x_x_x_zz, g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xy_xx[k] = -g_x_x_x_xx[k] * ab_y + g_x_x_x_xxy[k];

                g_x_x_xy_xy[k] = -g_x_x_x_xy[k] * ab_y + g_x_x_x_xyy[k];

                g_x_x_xy_xz[k] = -g_x_x_x_xz[k] * ab_y + g_x_x_x_xyz[k];

                g_x_x_xy_yy[k] = -g_x_x_x_yy[k] * ab_y + g_x_x_x_yyy[k];

                g_x_x_xy_yz[k] = -g_x_x_x_yz[k] * ab_y + g_x_x_x_yyz[k];

                g_x_x_xy_zz[k] = -g_x_x_x_zz[k] * ab_y + g_x_x_x_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_x_xz_xx = cbuffer.data(dd_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xz_xy = cbuffer.data(dd_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xz_xz = cbuffer.data(dd_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xz_yy = cbuffer.data(dd_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xz_yz = cbuffer.data(dd_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xz_zz = cbuffer.data(dd_geom_11_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xx, g_x_x_x_xxz, g_x_x_x_xy, g_x_x_x_xyz, g_x_x_x_xz, g_x_x_x_xzz, g_x_x_x_yy, g_x_x_x_yyz, g_x_x_x_yz, g_x_x_x_yzz, g_x_x_x_zz, g_x_x_x_zzz, g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xz_xx[k] = -g_x_x_x_xx[k] * ab_z + g_x_x_x_xxz[k];

                g_x_x_xz_xy[k] = -g_x_x_x_xy[k] * ab_z + g_x_x_x_xyz[k];

                g_x_x_xz_xz[k] = -g_x_x_x_xz[k] * ab_z + g_x_x_x_xzz[k];

                g_x_x_xz_yy[k] = -g_x_x_x_yy[k] * ab_z + g_x_x_x_yyz[k];

                g_x_x_xz_yz[k] = -g_x_x_x_yz[k] * ab_z + g_x_x_x_yzz[k];

                g_x_x_xz_zz[k] = -g_x_x_x_zz[k] * ab_z + g_x_x_x_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_x_x_yy_xx = cbuffer.data(dd_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_yy_xy = cbuffer.data(dd_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_yy_xz = cbuffer.data(dd_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_yy_yy = cbuffer.data(dd_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_yy_yz = cbuffer.data(dd_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_yy_zz = cbuffer.data(dd_geom_11_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_y_xx, g_x_x_y_xxy, g_x_x_y_xy, g_x_x_y_xyy, g_x_x_y_xyz, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yyy, g_x_x_y_yyz, g_x_x_y_yz, g_x_x_y_yzz, g_x_x_y_zz, g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yy_xx[k] = -g_x_x_y_xx[k] * ab_y + g_x_x_y_xxy[k];

                g_x_x_yy_xy[k] = -g_x_x_y_xy[k] * ab_y + g_x_x_y_xyy[k];

                g_x_x_yy_xz[k] = -g_x_x_y_xz[k] * ab_y + g_x_x_y_xyz[k];

                g_x_x_yy_yy[k] = -g_x_x_y_yy[k] * ab_y + g_x_x_y_yyy[k];

                g_x_x_yy_yz[k] = -g_x_x_y_yz[k] * ab_y + g_x_x_y_yyz[k];

                g_x_x_yy_zz[k] = -g_x_x_y_zz[k] * ab_y + g_x_x_y_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_yz_xx = cbuffer.data(dd_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_yz_xy = cbuffer.data(dd_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_yz_xz = cbuffer.data(dd_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_yz_yy = cbuffer.data(dd_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_yz_yz = cbuffer.data(dd_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_yz_zz = cbuffer.data(dd_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, g_x_x_z_xx, g_x_x_z_xxy, g_x_x_z_xy, g_x_x_z_xyy, g_x_x_z_xyz, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yyy, g_x_x_z_yyz, g_x_x_z_yz, g_x_x_z_yzz, g_x_x_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yz_xx[k] = -g_x_x_z_xx[k] * ab_y + g_x_x_z_xxy[k];

                g_x_x_yz_xy[k] = -g_x_x_z_xy[k] * ab_y + g_x_x_z_xyy[k];

                g_x_x_yz_xz[k] = -g_x_x_z_xz[k] * ab_y + g_x_x_z_xyz[k];

                g_x_x_yz_yy[k] = -g_x_x_z_yy[k] * ab_y + g_x_x_z_yyy[k];

                g_x_x_yz_yz[k] = -g_x_x_z_yz[k] * ab_y + g_x_x_z_yyz[k];

                g_x_x_yz_zz[k] = -g_x_x_z_zz[k] * ab_y + g_x_x_z_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_x_x_zz_xx = cbuffer.data(dd_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_zz_xy = cbuffer.data(dd_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_zz_xz = cbuffer.data(dd_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_zz_yy = cbuffer.data(dd_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_zz_yz = cbuffer.data(dd_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_zz_zz = cbuffer.data(dd_geom_11_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_z_xx, g_x_x_z_xxz, g_x_x_z_xy, g_x_x_z_xyz, g_x_x_z_xz, g_x_x_z_xzz, g_x_x_z_yy, g_x_x_z_yyz, g_x_x_z_yz, g_x_x_z_yzz, g_x_x_z_zz, g_x_x_z_zzz, g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zz_xx[k] = -g_x_x_z_xx[k] * ab_z + g_x_x_z_xxz[k];

                g_x_x_zz_xy[k] = -g_x_x_z_xy[k] * ab_z + g_x_x_z_xyz[k];

                g_x_x_zz_xz[k] = -g_x_x_z_xz[k] * ab_z + g_x_x_z_xzz[k];

                g_x_x_zz_yy[k] = -g_x_x_z_yy[k] * ab_z + g_x_x_z_yyz[k];

                g_x_x_zz_yz[k] = -g_x_x_z_yz[k] * ab_z + g_x_x_z_yzz[k];

                g_x_x_zz_zz[k] = -g_x_x_z_zz[k] * ab_z + g_x_x_z_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_x_y_xx_xx = cbuffer.data(dd_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_xx_xy = cbuffer.data(dd_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_xx_xz = cbuffer.data(dd_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_y_xx_yy = cbuffer.data(dd_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_y_xx_yz = cbuffer.data(dd_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_xx_zz = cbuffer.data(dd_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_x_y_x_xx, g_x_y_x_xxx, g_x_y_x_xxy, g_x_y_x_xxz, g_x_y_x_xy, g_x_y_x_xyy, g_x_y_x_xyz, g_x_y_x_xz, g_x_y_x_xzz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xx_xx[k] = -g_0_y_x_xx[k] - g_x_y_x_xx[k] * ab_x + g_x_y_x_xxx[k];

                g_x_y_xx_xy[k] = -g_0_y_x_xy[k] - g_x_y_x_xy[k] * ab_x + g_x_y_x_xxy[k];

                g_x_y_xx_xz[k] = -g_0_y_x_xz[k] - g_x_y_x_xz[k] * ab_x + g_x_y_x_xxz[k];

                g_x_y_xx_yy[k] = -g_0_y_x_yy[k] - g_x_y_x_yy[k] * ab_x + g_x_y_x_xyy[k];

                g_x_y_xx_yz[k] = -g_0_y_x_yz[k] - g_x_y_x_yz[k] * ab_x + g_x_y_x_xyz[k];

                g_x_y_xx_zz[k] = -g_0_y_x_zz[k] - g_x_y_x_zz[k] * ab_x + g_x_y_x_xzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_x_y_xy_xx = cbuffer.data(dd_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_y_xy_xy = cbuffer.data(dd_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_y_xy_xz = cbuffer.data(dd_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_xy_yy = cbuffer.data(dd_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_xy_yz = cbuffer.data(dd_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_xy_zz = cbuffer.data(dd_geom_11_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, g_x_y_y_xx, g_x_y_y_xxx, g_x_y_y_xxy, g_x_y_y_xxz, g_x_y_y_xy, g_x_y_y_xyy, g_x_y_y_xyz, g_x_y_y_xz, g_x_y_y_xzz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xy_xx[k] = -g_0_y_y_xx[k] - g_x_y_y_xx[k] * ab_x + g_x_y_y_xxx[k];

                g_x_y_xy_xy[k] = -g_0_y_y_xy[k] - g_x_y_y_xy[k] * ab_x + g_x_y_y_xxy[k];

                g_x_y_xy_xz[k] = -g_0_y_y_xz[k] - g_x_y_y_xz[k] * ab_x + g_x_y_y_xxz[k];

                g_x_y_xy_yy[k] = -g_0_y_y_yy[k] - g_x_y_y_yy[k] * ab_x + g_x_y_y_xyy[k];

                g_x_y_xy_yz[k] = -g_0_y_y_yz[k] - g_x_y_y_yz[k] * ab_x + g_x_y_y_xyz[k];

                g_x_y_xy_zz[k] = -g_0_y_y_zz[k] - g_x_y_y_zz[k] * ab_x + g_x_y_y_xzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_x_y_xz_xx = cbuffer.data(dd_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_xz_xy = cbuffer.data(dd_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_xz_xz = cbuffer.data(dd_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_xz_yy = cbuffer.data(dd_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_xz_yz = cbuffer.data(dd_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_xz_zz = cbuffer.data(dd_geom_11_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_x_xx, g_x_y_x_xxz, g_x_y_x_xy, g_x_y_x_xyz, g_x_y_x_xz, g_x_y_x_xzz, g_x_y_x_yy, g_x_y_x_yyz, g_x_y_x_yz, g_x_y_x_yzz, g_x_y_x_zz, g_x_y_x_zzz, g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xz_xx[k] = -g_x_y_x_xx[k] * ab_z + g_x_y_x_xxz[k];

                g_x_y_xz_xy[k] = -g_x_y_x_xy[k] * ab_z + g_x_y_x_xyz[k];

                g_x_y_xz_xz[k] = -g_x_y_x_xz[k] * ab_z + g_x_y_x_xzz[k];

                g_x_y_xz_yy[k] = -g_x_y_x_yy[k] * ab_z + g_x_y_x_yyz[k];

                g_x_y_xz_yz[k] = -g_x_y_x_yz[k] * ab_z + g_x_y_x_yzz[k];

                g_x_y_xz_zz[k] = -g_x_y_x_zz[k] * ab_z + g_x_y_x_zzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_x_y_yy_xx = cbuffer.data(dd_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_yy_xy = cbuffer.data(dd_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_yy_xz = cbuffer.data(dd_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_y_yy_yy = cbuffer.data(dd_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_yy_yz = cbuffer.data(dd_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_yy_zz = cbuffer.data(dd_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz, g_x_y_y_xx, g_x_y_y_xxy, g_x_y_y_xy, g_x_y_y_xyy, g_x_y_y_xyz, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yyy, g_x_y_y_yyz, g_x_y_y_yz, g_x_y_y_yzz, g_x_y_y_zz, g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yy_xx[k] = g_x_0_y_xx[k] - g_x_y_y_xx[k] * ab_y + g_x_y_y_xxy[k];

                g_x_y_yy_xy[k] = g_x_0_y_xy[k] - g_x_y_y_xy[k] * ab_y + g_x_y_y_xyy[k];

                g_x_y_yy_xz[k] = g_x_0_y_xz[k] - g_x_y_y_xz[k] * ab_y + g_x_y_y_xyz[k];

                g_x_y_yy_yy[k] = g_x_0_y_yy[k] - g_x_y_y_yy[k] * ab_y + g_x_y_y_yyy[k];

                g_x_y_yy_yz[k] = g_x_0_y_yz[k] - g_x_y_y_yz[k] * ab_y + g_x_y_y_yyz[k];

                g_x_y_yy_zz[k] = g_x_0_y_zz[k] - g_x_y_y_zz[k] * ab_y + g_x_y_y_yzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_x_y_yz_xx = cbuffer.data(dd_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_yz_xy = cbuffer.data(dd_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_yz_xz = cbuffer.data(dd_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_yz_yy = cbuffer.data(dd_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_yz_yz = cbuffer.data(dd_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_yz_zz = cbuffer.data(dd_geom_11_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_y_xx, g_x_y_y_xxz, g_x_y_y_xy, g_x_y_y_xyz, g_x_y_y_xz, g_x_y_y_xzz, g_x_y_y_yy, g_x_y_y_yyz, g_x_y_y_yz, g_x_y_y_yzz, g_x_y_y_zz, g_x_y_y_zzz, g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yz_xx[k] = -g_x_y_y_xx[k] * ab_z + g_x_y_y_xxz[k];

                g_x_y_yz_xy[k] = -g_x_y_y_xy[k] * ab_z + g_x_y_y_xyz[k];

                g_x_y_yz_xz[k] = -g_x_y_y_xz[k] * ab_z + g_x_y_y_xzz[k];

                g_x_y_yz_yy[k] = -g_x_y_y_yy[k] * ab_z + g_x_y_y_yyz[k];

                g_x_y_yz_yz[k] = -g_x_y_y_yz[k] * ab_z + g_x_y_y_yzz[k];

                g_x_y_yz_zz[k] = -g_x_y_y_zz[k] * ab_z + g_x_y_y_zzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_x_y_zz_xx = cbuffer.data(dd_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_zz_xy = cbuffer.data(dd_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_zz_xz = cbuffer.data(dd_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_zz_yy = cbuffer.data(dd_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_zz_yz = cbuffer.data(dd_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_zz_zz = cbuffer.data(dd_geom_11_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_z_xx, g_x_y_z_xxz, g_x_y_z_xy, g_x_y_z_xyz, g_x_y_z_xz, g_x_y_z_xzz, g_x_y_z_yy, g_x_y_z_yyz, g_x_y_z_yz, g_x_y_z_yzz, g_x_y_z_zz, g_x_y_z_zzz, g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zz_xx[k] = -g_x_y_z_xx[k] * ab_z + g_x_y_z_xxz[k];

                g_x_y_zz_xy[k] = -g_x_y_z_xy[k] * ab_z + g_x_y_z_xyz[k];

                g_x_y_zz_xz[k] = -g_x_y_z_xz[k] * ab_z + g_x_y_z_xzz[k];

                g_x_y_zz_yy[k] = -g_x_y_z_yy[k] * ab_z + g_x_y_z_yyz[k];

                g_x_y_zz_yz[k] = -g_x_y_z_yz[k] * ab_z + g_x_y_z_yzz[k];

                g_x_y_zz_zz[k] = -g_x_y_z_zz[k] * ab_z + g_x_y_z_zzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_x_z_xx_xx = cbuffer.data(dd_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_z_xx_xy = cbuffer.data(dd_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_xx_xz = cbuffer.data(dd_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_z_xx_yy = cbuffer.data(dd_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_z_xx_yz = cbuffer.data(dd_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_xx_zz = cbuffer.data(dd_geom_11_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_x_z_x_xx, g_x_z_x_xxx, g_x_z_x_xxy, g_x_z_x_xxz, g_x_z_x_xy, g_x_z_x_xyy, g_x_z_x_xyz, g_x_z_x_xz, g_x_z_x_xzz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xx_xx[k] = -g_0_z_x_xx[k] - g_x_z_x_xx[k] * ab_x + g_x_z_x_xxx[k];

                g_x_z_xx_xy[k] = -g_0_z_x_xy[k] - g_x_z_x_xy[k] * ab_x + g_x_z_x_xxy[k];

                g_x_z_xx_xz[k] = -g_0_z_x_xz[k] - g_x_z_x_xz[k] * ab_x + g_x_z_x_xxz[k];

                g_x_z_xx_yy[k] = -g_0_z_x_yy[k] - g_x_z_x_yy[k] * ab_x + g_x_z_x_xyy[k];

                g_x_z_xx_yz[k] = -g_0_z_x_yz[k] - g_x_z_x_yz[k] * ab_x + g_x_z_x_xyz[k];

                g_x_z_xx_zz[k] = -g_0_z_x_zz[k] - g_x_z_x_zz[k] * ab_x + g_x_z_x_xzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_x_z_xy_xx = cbuffer.data(dd_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_z_xy_xy = cbuffer.data(dd_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_z_xy_xz = cbuffer.data(dd_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_z_xy_yy = cbuffer.data(dd_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_z_xy_yz = cbuffer.data(dd_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_z_xy_zz = cbuffer.data(dd_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_x_xx, g_x_z_x_xxy, g_x_z_x_xy, g_x_z_x_xyy, g_x_z_x_xyz, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yyy, g_x_z_x_yyz, g_x_z_x_yz, g_x_z_x_yzz, g_x_z_x_zz, g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xy_xx[k] = -g_x_z_x_xx[k] * ab_y + g_x_z_x_xxy[k];

                g_x_z_xy_xy[k] = -g_x_z_x_xy[k] * ab_y + g_x_z_x_xyy[k];

                g_x_z_xy_xz[k] = -g_x_z_x_xz[k] * ab_y + g_x_z_x_xyz[k];

                g_x_z_xy_yy[k] = -g_x_z_x_yy[k] * ab_y + g_x_z_x_yyy[k];

                g_x_z_xy_yz[k] = -g_x_z_x_yz[k] * ab_y + g_x_z_x_yyz[k];

                g_x_z_xy_zz[k] = -g_x_z_x_zz[k] * ab_y + g_x_z_x_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_x_z_xz_xx = cbuffer.data(dd_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_z_xz_xy = cbuffer.data(dd_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_z_xz_xz = cbuffer.data(dd_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_z_xz_yy = cbuffer.data(dd_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_z_xz_yz = cbuffer.data(dd_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_z_xz_zz = cbuffer.data(dd_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz, g_x_z_z_xx, g_x_z_z_xxx, g_x_z_z_xxy, g_x_z_z_xxz, g_x_z_z_xy, g_x_z_z_xyy, g_x_z_z_xyz, g_x_z_z_xz, g_x_z_z_xzz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xz_xx[k] = -g_0_z_z_xx[k] - g_x_z_z_xx[k] * ab_x + g_x_z_z_xxx[k];

                g_x_z_xz_xy[k] = -g_0_z_z_xy[k] - g_x_z_z_xy[k] * ab_x + g_x_z_z_xxy[k];

                g_x_z_xz_xz[k] = -g_0_z_z_xz[k] - g_x_z_z_xz[k] * ab_x + g_x_z_z_xxz[k];

                g_x_z_xz_yy[k] = -g_0_z_z_yy[k] - g_x_z_z_yy[k] * ab_x + g_x_z_z_xyy[k];

                g_x_z_xz_yz[k] = -g_0_z_z_yz[k] - g_x_z_z_yz[k] * ab_x + g_x_z_z_xyz[k];

                g_x_z_xz_zz[k] = -g_0_z_z_zz[k] - g_x_z_z_zz[k] * ab_x + g_x_z_z_xzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_x_z_yy_xx = cbuffer.data(dd_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_z_yy_xy = cbuffer.data(dd_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_z_yy_xz = cbuffer.data(dd_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_z_yy_yy = cbuffer.data(dd_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_z_yy_yz = cbuffer.data(dd_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_z_yy_zz = cbuffer.data(dd_geom_11_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_y_xx, g_x_z_y_xxy, g_x_z_y_xy, g_x_z_y_xyy, g_x_z_y_xyz, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yyy, g_x_z_y_yyz, g_x_z_y_yz, g_x_z_y_yzz, g_x_z_y_zz, g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yy_xx[k] = -g_x_z_y_xx[k] * ab_y + g_x_z_y_xxy[k];

                g_x_z_yy_xy[k] = -g_x_z_y_xy[k] * ab_y + g_x_z_y_xyy[k];

                g_x_z_yy_xz[k] = -g_x_z_y_xz[k] * ab_y + g_x_z_y_xyz[k];

                g_x_z_yy_yy[k] = -g_x_z_y_yy[k] * ab_y + g_x_z_y_yyy[k];

                g_x_z_yy_yz[k] = -g_x_z_y_yz[k] * ab_y + g_x_z_y_yyz[k];

                g_x_z_yy_zz[k] = -g_x_z_y_zz[k] * ab_y + g_x_z_y_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_x_z_yz_xx = cbuffer.data(dd_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_z_yz_xy = cbuffer.data(dd_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_z_yz_xz = cbuffer.data(dd_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_z_yz_yy = cbuffer.data(dd_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_z_yz_yz = cbuffer.data(dd_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_z_yz_zz = cbuffer.data(dd_geom_11_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz, g_x_z_z_xx, g_x_z_z_xxy, g_x_z_z_xy, g_x_z_z_xyy, g_x_z_z_xyz, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yyy, g_x_z_z_yyz, g_x_z_z_yz, g_x_z_z_yzz, g_x_z_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yz_xx[k] = -g_x_z_z_xx[k] * ab_y + g_x_z_z_xxy[k];

                g_x_z_yz_xy[k] = -g_x_z_z_xy[k] * ab_y + g_x_z_z_xyy[k];

                g_x_z_yz_xz[k] = -g_x_z_z_xz[k] * ab_y + g_x_z_z_xyz[k];

                g_x_z_yz_yy[k] = -g_x_z_z_yy[k] * ab_y + g_x_z_z_yyy[k];

                g_x_z_yz_yz[k] = -g_x_z_z_yz[k] * ab_y + g_x_z_z_yyz[k];

                g_x_z_yz_zz[k] = -g_x_z_z_zz[k] * ab_y + g_x_z_z_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_x_z_zz_xx = cbuffer.data(dd_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_z_zz_xy = cbuffer.data(dd_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_z_zz_xz = cbuffer.data(dd_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_z_zz_yy = cbuffer.data(dd_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_z_zz_yz = cbuffer.data(dd_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_z_zz_zz = cbuffer.data(dd_geom_11_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz, g_x_z_z_xx, g_x_z_z_xxz, g_x_z_z_xy, g_x_z_z_xyz, g_x_z_z_xz, g_x_z_z_xzz, g_x_z_z_yy, g_x_z_z_yyz, g_x_z_z_yz, g_x_z_z_yzz, g_x_z_z_zz, g_x_z_z_zzz, g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zz_xx[k] = g_x_0_z_xx[k] - g_x_z_z_xx[k] * ab_z + g_x_z_z_xxz[k];

                g_x_z_zz_xy[k] = g_x_0_z_xy[k] - g_x_z_z_xy[k] * ab_z + g_x_z_z_xyz[k];

                g_x_z_zz_xz[k] = g_x_0_z_xz[k] - g_x_z_z_xz[k] * ab_z + g_x_z_z_xzz[k];

                g_x_z_zz_yy[k] = g_x_0_z_yy[k] - g_x_z_z_yy[k] * ab_z + g_x_z_z_yyz[k];

                g_x_z_zz_yz[k] = g_x_0_z_yz[k] - g_x_z_z_yz[k] * ab_z + g_x_z_z_yzz[k];

                g_x_z_zz_zz[k] = g_x_0_z_zz[k] - g_x_z_z_zz[k] * ab_z + g_x_z_z_zzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_y_x_xx_xx = cbuffer.data(dd_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_x_xx_xy = cbuffer.data(dd_geom_11_off + 109 * ccomps * dcomps);

            auto g_y_x_xx_xz = cbuffer.data(dd_geom_11_off + 110 * ccomps * dcomps);

            auto g_y_x_xx_yy = cbuffer.data(dd_geom_11_off + 111 * ccomps * dcomps);

            auto g_y_x_xx_yz = cbuffer.data(dd_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_x_xx_zz = cbuffer.data(dd_geom_11_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, g_y_x_x_xx, g_y_x_x_xxx, g_y_x_x_xxy, g_y_x_x_xxz, g_y_x_x_xy, g_y_x_x_xyy, g_y_x_x_xyz, g_y_x_x_xz, g_y_x_x_xzz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xx_xx[k] = g_y_0_x_xx[k] - g_y_x_x_xx[k] * ab_x + g_y_x_x_xxx[k];

                g_y_x_xx_xy[k] = g_y_0_x_xy[k] - g_y_x_x_xy[k] * ab_x + g_y_x_x_xxy[k];

                g_y_x_xx_xz[k] = g_y_0_x_xz[k] - g_y_x_x_xz[k] * ab_x + g_y_x_x_xxz[k];

                g_y_x_xx_yy[k] = g_y_0_x_yy[k] - g_y_x_x_yy[k] * ab_x + g_y_x_x_xyy[k];

                g_y_x_xx_yz[k] = g_y_0_x_yz[k] - g_y_x_x_yz[k] * ab_x + g_y_x_x_xyz[k];

                g_y_x_xx_zz[k] = g_y_0_x_zz[k] - g_y_x_x_zz[k] * ab_x + g_y_x_x_xzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_y_x_xy_xx = cbuffer.data(dd_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_x_xy_xy = cbuffer.data(dd_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_x_xy_xz = cbuffer.data(dd_geom_11_off + 116 * ccomps * dcomps);

            auto g_y_x_xy_yy = cbuffer.data(dd_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_x_xy_yz = cbuffer.data(dd_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_x_xy_zz = cbuffer.data(dd_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz, g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz, g_y_x_y_xx, g_y_x_y_xxx, g_y_x_y_xxy, g_y_x_y_xxz, g_y_x_y_xy, g_y_x_y_xyy, g_y_x_y_xyz, g_y_x_y_xz, g_y_x_y_xzz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xy_xx[k] = g_y_0_y_xx[k] - g_y_x_y_xx[k] * ab_x + g_y_x_y_xxx[k];

                g_y_x_xy_xy[k] = g_y_0_y_xy[k] - g_y_x_y_xy[k] * ab_x + g_y_x_y_xxy[k];

                g_y_x_xy_xz[k] = g_y_0_y_xz[k] - g_y_x_y_xz[k] * ab_x + g_y_x_y_xxz[k];

                g_y_x_xy_yy[k] = g_y_0_y_yy[k] - g_y_x_y_yy[k] * ab_x + g_y_x_y_xyy[k];

                g_y_x_xy_yz[k] = g_y_0_y_yz[k] - g_y_x_y_yz[k] * ab_x + g_y_x_y_xyz[k];

                g_y_x_xy_zz[k] = g_y_0_y_zz[k] - g_y_x_y_zz[k] * ab_x + g_y_x_y_xzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_y_x_xz_xx = cbuffer.data(dd_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_x_xz_xy = cbuffer.data(dd_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_x_xz_xz = cbuffer.data(dd_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_x_xz_yy = cbuffer.data(dd_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_x_xz_yz = cbuffer.data(dd_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_x_xz_zz = cbuffer.data(dd_geom_11_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_x_xx, g_y_x_x_xxz, g_y_x_x_xy, g_y_x_x_xyz, g_y_x_x_xz, g_y_x_x_xzz, g_y_x_x_yy, g_y_x_x_yyz, g_y_x_x_yz, g_y_x_x_yzz, g_y_x_x_zz, g_y_x_x_zzz, g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xz_xx[k] = -g_y_x_x_xx[k] * ab_z + g_y_x_x_xxz[k];

                g_y_x_xz_xy[k] = -g_y_x_x_xy[k] * ab_z + g_y_x_x_xyz[k];

                g_y_x_xz_xz[k] = -g_y_x_x_xz[k] * ab_z + g_y_x_x_xzz[k];

                g_y_x_xz_yy[k] = -g_y_x_x_yy[k] * ab_z + g_y_x_x_yyz[k];

                g_y_x_xz_yz[k] = -g_y_x_x_yz[k] * ab_z + g_y_x_x_yzz[k];

                g_y_x_xz_zz[k] = -g_y_x_x_zz[k] * ab_z + g_y_x_x_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_y_x_yy_xx = cbuffer.data(dd_geom_11_off + 126 * ccomps * dcomps);

            auto g_y_x_yy_xy = cbuffer.data(dd_geom_11_off + 127 * ccomps * dcomps);

            auto g_y_x_yy_xz = cbuffer.data(dd_geom_11_off + 128 * ccomps * dcomps);

            auto g_y_x_yy_yy = cbuffer.data(dd_geom_11_off + 129 * ccomps * dcomps);

            auto g_y_x_yy_yz = cbuffer.data(dd_geom_11_off + 130 * ccomps * dcomps);

            auto g_y_x_yy_zz = cbuffer.data(dd_geom_11_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_y_x_y_xx, g_y_x_y_xxy, g_y_x_y_xy, g_y_x_y_xyy, g_y_x_y_xyz, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yyy, g_y_x_y_yyz, g_y_x_y_yz, g_y_x_y_yzz, g_y_x_y_zz, g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yy_xx[k] = -g_0_x_y_xx[k] - g_y_x_y_xx[k] * ab_y + g_y_x_y_xxy[k];

                g_y_x_yy_xy[k] = -g_0_x_y_xy[k] - g_y_x_y_xy[k] * ab_y + g_y_x_y_xyy[k];

                g_y_x_yy_xz[k] = -g_0_x_y_xz[k] - g_y_x_y_xz[k] * ab_y + g_y_x_y_xyz[k];

                g_y_x_yy_yy[k] = -g_0_x_y_yy[k] - g_y_x_y_yy[k] * ab_y + g_y_x_y_yyy[k];

                g_y_x_yy_yz[k] = -g_0_x_y_yz[k] - g_y_x_y_yz[k] * ab_y + g_y_x_y_yyz[k];

                g_y_x_yy_zz[k] = -g_0_x_y_zz[k] - g_y_x_y_zz[k] * ab_y + g_y_x_y_yzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_y_x_yz_xx = cbuffer.data(dd_geom_11_off + 132 * ccomps * dcomps);

            auto g_y_x_yz_xy = cbuffer.data(dd_geom_11_off + 133 * ccomps * dcomps);

            auto g_y_x_yz_xz = cbuffer.data(dd_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_x_yz_yy = cbuffer.data(dd_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_x_yz_yz = cbuffer.data(dd_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_x_yz_zz = cbuffer.data(dd_geom_11_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_y_xx, g_y_x_y_xxz, g_y_x_y_xy, g_y_x_y_xyz, g_y_x_y_xz, g_y_x_y_xzz, g_y_x_y_yy, g_y_x_y_yyz, g_y_x_y_yz, g_y_x_y_yzz, g_y_x_y_zz, g_y_x_y_zzz, g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yz_xx[k] = -g_y_x_y_xx[k] * ab_z + g_y_x_y_xxz[k];

                g_y_x_yz_xy[k] = -g_y_x_y_xy[k] * ab_z + g_y_x_y_xyz[k];

                g_y_x_yz_xz[k] = -g_y_x_y_xz[k] * ab_z + g_y_x_y_xzz[k];

                g_y_x_yz_yy[k] = -g_y_x_y_yy[k] * ab_z + g_y_x_y_yyz[k];

                g_y_x_yz_yz[k] = -g_y_x_y_yz[k] * ab_z + g_y_x_y_yzz[k];

                g_y_x_yz_zz[k] = -g_y_x_y_zz[k] * ab_z + g_y_x_y_zzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_y_x_zz_xx = cbuffer.data(dd_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_x_zz_xy = cbuffer.data(dd_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_x_zz_xz = cbuffer.data(dd_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_x_zz_yy = cbuffer.data(dd_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_x_zz_yz = cbuffer.data(dd_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_x_zz_zz = cbuffer.data(dd_geom_11_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_z_xx, g_y_x_z_xxz, g_y_x_z_xy, g_y_x_z_xyz, g_y_x_z_xz, g_y_x_z_xzz, g_y_x_z_yy, g_y_x_z_yyz, g_y_x_z_yz, g_y_x_z_yzz, g_y_x_z_zz, g_y_x_z_zzz, g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zz_xx[k] = -g_y_x_z_xx[k] * ab_z + g_y_x_z_xxz[k];

                g_y_x_zz_xy[k] = -g_y_x_z_xy[k] * ab_z + g_y_x_z_xyz[k];

                g_y_x_zz_xz[k] = -g_y_x_z_xz[k] * ab_z + g_y_x_z_xzz[k];

                g_y_x_zz_yy[k] = -g_y_x_z_yy[k] * ab_z + g_y_x_z_yyz[k];

                g_y_x_zz_yz[k] = -g_y_x_z_yz[k] * ab_z + g_y_x_z_yzz[k];

                g_y_x_zz_zz[k] = -g_y_x_z_zz[k] * ab_z + g_y_x_z_zzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_y_y_xx_xx = cbuffer.data(dd_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_y_xx_xy = cbuffer.data(dd_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_y_xx_xz = cbuffer.data(dd_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_y_xx_yy = cbuffer.data(dd_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_y_xx_yz = cbuffer.data(dd_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_y_xx_zz = cbuffer.data(dd_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_x_xx, g_y_y_x_xxx, g_y_y_x_xxy, g_y_y_x_xxz, g_y_y_x_xy, g_y_y_x_xyy, g_y_y_x_xyz, g_y_y_x_xz, g_y_y_x_xzz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xx_xx[k] = -g_y_y_x_xx[k] * ab_x + g_y_y_x_xxx[k];

                g_y_y_xx_xy[k] = -g_y_y_x_xy[k] * ab_x + g_y_y_x_xxy[k];

                g_y_y_xx_xz[k] = -g_y_y_x_xz[k] * ab_x + g_y_y_x_xxz[k];

                g_y_y_xx_yy[k] = -g_y_y_x_yy[k] * ab_x + g_y_y_x_xyy[k];

                g_y_y_xx_yz[k] = -g_y_y_x_yz[k] * ab_x + g_y_y_x_xyz[k];

                g_y_y_xx_zz[k] = -g_y_y_x_zz[k] * ab_x + g_y_y_x_xzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_y_y_xy_xx = cbuffer.data(dd_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_y_xy_xy = cbuffer.data(dd_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_y_xy_xz = cbuffer.data(dd_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_y_xy_yy = cbuffer.data(dd_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_y_xy_yz = cbuffer.data(dd_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_y_xy_zz = cbuffer.data(dd_geom_11_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz, g_y_y_y_xx, g_y_y_y_xxx, g_y_y_y_xxy, g_y_y_y_xxz, g_y_y_y_xy, g_y_y_y_xyy, g_y_y_y_xyz, g_y_y_y_xz, g_y_y_y_xzz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xy_xx[k] = -g_y_y_y_xx[k] * ab_x + g_y_y_y_xxx[k];

                g_y_y_xy_xy[k] = -g_y_y_y_xy[k] * ab_x + g_y_y_y_xxy[k];

                g_y_y_xy_xz[k] = -g_y_y_y_xz[k] * ab_x + g_y_y_y_xxz[k];

                g_y_y_xy_yy[k] = -g_y_y_y_yy[k] * ab_x + g_y_y_y_xyy[k];

                g_y_y_xy_yz[k] = -g_y_y_y_yz[k] * ab_x + g_y_y_y_xyz[k];

                g_y_y_xy_zz[k] = -g_y_y_y_zz[k] * ab_x + g_y_y_y_xzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_y_y_xz_xx = cbuffer.data(dd_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_y_xz_xy = cbuffer.data(dd_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_y_xz_xz = cbuffer.data(dd_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_y_xz_yy = cbuffer.data(dd_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_y_xz_yz = cbuffer.data(dd_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_y_xz_zz = cbuffer.data(dd_geom_11_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz, g_y_y_z_xx, g_y_y_z_xxx, g_y_y_z_xxy, g_y_y_z_xxz, g_y_y_z_xy, g_y_y_z_xyy, g_y_y_z_xyz, g_y_y_z_xz, g_y_y_z_xzz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xz_xx[k] = -g_y_y_z_xx[k] * ab_x + g_y_y_z_xxx[k];

                g_y_y_xz_xy[k] = -g_y_y_z_xy[k] * ab_x + g_y_y_z_xxy[k];

                g_y_y_xz_xz[k] = -g_y_y_z_xz[k] * ab_x + g_y_y_z_xxz[k];

                g_y_y_xz_yy[k] = -g_y_y_z_yy[k] * ab_x + g_y_y_z_xyy[k];

                g_y_y_xz_yz[k] = -g_y_y_z_yz[k] * ab_x + g_y_y_z_xyz[k];

                g_y_y_xz_zz[k] = -g_y_y_z_zz[k] * ab_x + g_y_y_z_xzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_y_y_yy_xx = cbuffer.data(dd_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_y_yy_xy = cbuffer.data(dd_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_y_yy_xz = cbuffer.data(dd_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_y_yy_yy = cbuffer.data(dd_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_y_yy_yz = cbuffer.data(dd_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_y_yy_zz = cbuffer.data(dd_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz, g_y_y_y_xx, g_y_y_y_xxy, g_y_y_y_xy, g_y_y_y_xyy, g_y_y_y_xyz, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yyy, g_y_y_y_yyz, g_y_y_y_yz, g_y_y_y_yzz, g_y_y_y_zz, g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yy_xx[k] = -g_0_y_y_xx[k] + g_y_0_y_xx[k] - g_y_y_y_xx[k] * ab_y + g_y_y_y_xxy[k];

                g_y_y_yy_xy[k] = -g_0_y_y_xy[k] + g_y_0_y_xy[k] - g_y_y_y_xy[k] * ab_y + g_y_y_y_xyy[k];

                g_y_y_yy_xz[k] = -g_0_y_y_xz[k] + g_y_0_y_xz[k] - g_y_y_y_xz[k] * ab_y + g_y_y_y_xyz[k];

                g_y_y_yy_yy[k] = -g_0_y_y_yy[k] + g_y_0_y_yy[k] - g_y_y_y_yy[k] * ab_y + g_y_y_y_yyy[k];

                g_y_y_yy_yz[k] = -g_0_y_y_yz[k] + g_y_0_y_yz[k] - g_y_y_y_yz[k] * ab_y + g_y_y_y_yyz[k];

                g_y_y_yy_zz[k] = -g_0_y_y_zz[k] + g_y_0_y_zz[k] - g_y_y_y_zz[k] * ab_y + g_y_y_y_yzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_y_y_yz_xx = cbuffer.data(dd_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_y_yz_xy = cbuffer.data(dd_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_y_yz_xz = cbuffer.data(dd_geom_11_off + 170 * ccomps * dcomps);

            auto g_y_y_yz_yy = cbuffer.data(dd_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_y_yz_yz = cbuffer.data(dd_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_y_yz_zz = cbuffer.data(dd_geom_11_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_y_xx, g_y_y_y_xxz, g_y_y_y_xy, g_y_y_y_xyz, g_y_y_y_xz, g_y_y_y_xzz, g_y_y_y_yy, g_y_y_y_yyz, g_y_y_y_yz, g_y_y_y_yzz, g_y_y_y_zz, g_y_y_y_zzz, g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yz_xx[k] = -g_y_y_y_xx[k] * ab_z + g_y_y_y_xxz[k];

                g_y_y_yz_xy[k] = -g_y_y_y_xy[k] * ab_z + g_y_y_y_xyz[k];

                g_y_y_yz_xz[k] = -g_y_y_y_xz[k] * ab_z + g_y_y_y_xzz[k];

                g_y_y_yz_yy[k] = -g_y_y_y_yy[k] * ab_z + g_y_y_y_yyz[k];

                g_y_y_yz_yz[k] = -g_y_y_y_yz[k] * ab_z + g_y_y_y_yzz[k];

                g_y_y_yz_zz[k] = -g_y_y_y_zz[k] * ab_z + g_y_y_y_zzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_y_y_zz_xx = cbuffer.data(dd_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_y_zz_xy = cbuffer.data(dd_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_y_zz_xz = cbuffer.data(dd_geom_11_off + 176 * ccomps * dcomps);

            auto g_y_y_zz_yy = cbuffer.data(dd_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_y_zz_yz = cbuffer.data(dd_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_y_zz_zz = cbuffer.data(dd_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_z_xx, g_y_y_z_xxz, g_y_y_z_xy, g_y_y_z_xyz, g_y_y_z_xz, g_y_y_z_xzz, g_y_y_z_yy, g_y_y_z_yyz, g_y_y_z_yz, g_y_y_z_yzz, g_y_y_z_zz, g_y_y_z_zzz, g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zz_xx[k] = -g_y_y_z_xx[k] * ab_z + g_y_y_z_xxz[k];

                g_y_y_zz_xy[k] = -g_y_y_z_xy[k] * ab_z + g_y_y_z_xyz[k];

                g_y_y_zz_xz[k] = -g_y_y_z_xz[k] * ab_z + g_y_y_z_xzz[k];

                g_y_y_zz_yy[k] = -g_y_y_z_yy[k] * ab_z + g_y_y_z_yyz[k];

                g_y_y_zz_yz[k] = -g_y_y_z_yz[k] * ab_z + g_y_y_z_yzz[k];

                g_y_y_zz_zz[k] = -g_y_y_z_zz[k] * ab_z + g_y_y_z_zzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_y_z_xx_xx = cbuffer.data(dd_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_z_xx_xy = cbuffer.data(dd_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_z_xx_xz = cbuffer.data(dd_geom_11_off + 182 * ccomps * dcomps);

            auto g_y_z_xx_yy = cbuffer.data(dd_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_z_xx_yz = cbuffer.data(dd_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_z_xx_zz = cbuffer.data(dd_geom_11_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_x_xx, g_y_z_x_xxx, g_y_z_x_xxy, g_y_z_x_xxz, g_y_z_x_xy, g_y_z_x_xyy, g_y_z_x_xyz, g_y_z_x_xz, g_y_z_x_xzz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz, g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xx_xx[k] = -g_y_z_x_xx[k] * ab_x + g_y_z_x_xxx[k];

                g_y_z_xx_xy[k] = -g_y_z_x_xy[k] * ab_x + g_y_z_x_xxy[k];

                g_y_z_xx_xz[k] = -g_y_z_x_xz[k] * ab_x + g_y_z_x_xxz[k];

                g_y_z_xx_yy[k] = -g_y_z_x_yy[k] * ab_x + g_y_z_x_xyy[k];

                g_y_z_xx_yz[k] = -g_y_z_x_yz[k] * ab_x + g_y_z_x_xyz[k];

                g_y_z_xx_zz[k] = -g_y_z_x_zz[k] * ab_x + g_y_z_x_xzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_y_z_xy_xx = cbuffer.data(dd_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_z_xy_xy = cbuffer.data(dd_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_z_xy_xz = cbuffer.data(dd_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_z_xy_yy = cbuffer.data(dd_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_z_xy_yz = cbuffer.data(dd_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_z_xy_zz = cbuffer.data(dd_geom_11_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz, g_y_z_y_xx, g_y_z_y_xxx, g_y_z_y_xxy, g_y_z_y_xxz, g_y_z_y_xy, g_y_z_y_xyy, g_y_z_y_xyz, g_y_z_y_xz, g_y_z_y_xzz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xy_xx[k] = -g_y_z_y_xx[k] * ab_x + g_y_z_y_xxx[k];

                g_y_z_xy_xy[k] = -g_y_z_y_xy[k] * ab_x + g_y_z_y_xxy[k];

                g_y_z_xy_xz[k] = -g_y_z_y_xz[k] * ab_x + g_y_z_y_xxz[k];

                g_y_z_xy_yy[k] = -g_y_z_y_yy[k] * ab_x + g_y_z_y_xyy[k];

                g_y_z_xy_yz[k] = -g_y_z_y_yz[k] * ab_x + g_y_z_y_xyz[k];

                g_y_z_xy_zz[k] = -g_y_z_y_zz[k] * ab_x + g_y_z_y_xzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_y_z_xz_xx = cbuffer.data(dd_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_z_xz_xy = cbuffer.data(dd_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_z_xz_xz = cbuffer.data(dd_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_z_xz_yy = cbuffer.data(dd_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_z_xz_yz = cbuffer.data(dd_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_z_xz_zz = cbuffer.data(dd_geom_11_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz, g_y_z_z_xx, g_y_z_z_xxx, g_y_z_z_xxy, g_y_z_z_xxz, g_y_z_z_xy, g_y_z_z_xyy, g_y_z_z_xyz, g_y_z_z_xz, g_y_z_z_xzz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xz_xx[k] = -g_y_z_z_xx[k] * ab_x + g_y_z_z_xxx[k];

                g_y_z_xz_xy[k] = -g_y_z_z_xy[k] * ab_x + g_y_z_z_xxy[k];

                g_y_z_xz_xz[k] = -g_y_z_z_xz[k] * ab_x + g_y_z_z_xxz[k];

                g_y_z_xz_yy[k] = -g_y_z_z_yy[k] * ab_x + g_y_z_z_xyy[k];

                g_y_z_xz_yz[k] = -g_y_z_z_yz[k] * ab_x + g_y_z_z_xyz[k];

                g_y_z_xz_zz[k] = -g_y_z_z_zz[k] * ab_x + g_y_z_z_xzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_y_z_yy_xx = cbuffer.data(dd_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_z_yy_xy = cbuffer.data(dd_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_z_yy_xz = cbuffer.data(dd_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_z_yy_yy = cbuffer.data(dd_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_z_yy_yz = cbuffer.data(dd_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_z_yy_zz = cbuffer.data(dd_geom_11_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_y_z_y_xx, g_y_z_y_xxy, g_y_z_y_xy, g_y_z_y_xyy, g_y_z_y_xyz, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yyy, g_y_z_y_yyz, g_y_z_y_yz, g_y_z_y_yzz, g_y_z_y_zz, g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yy_xx[k] = -g_0_z_y_xx[k] - g_y_z_y_xx[k] * ab_y + g_y_z_y_xxy[k];

                g_y_z_yy_xy[k] = -g_0_z_y_xy[k] - g_y_z_y_xy[k] * ab_y + g_y_z_y_xyy[k];

                g_y_z_yy_xz[k] = -g_0_z_y_xz[k] - g_y_z_y_xz[k] * ab_y + g_y_z_y_xyz[k];

                g_y_z_yy_yy[k] = -g_0_z_y_yy[k] - g_y_z_y_yy[k] * ab_y + g_y_z_y_yyy[k];

                g_y_z_yy_yz[k] = -g_0_z_y_yz[k] - g_y_z_y_yz[k] * ab_y + g_y_z_y_yyz[k];

                g_y_z_yy_zz[k] = -g_0_z_y_zz[k] - g_y_z_y_zz[k] * ab_y + g_y_z_y_yzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_y_z_yz_xx = cbuffer.data(dd_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_z_yz_xy = cbuffer.data(dd_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_z_yz_xz = cbuffer.data(dd_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_z_yz_yy = cbuffer.data(dd_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_z_yz_yz = cbuffer.data(dd_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_z_yz_zz = cbuffer.data(dd_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz, g_y_z_z_xx, g_y_z_z_xxy, g_y_z_z_xy, g_y_z_z_xyy, g_y_z_z_xyz, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yyy, g_y_z_z_yyz, g_y_z_z_yz, g_y_z_z_yzz, g_y_z_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yz_xx[k] = -g_0_z_z_xx[k] - g_y_z_z_xx[k] * ab_y + g_y_z_z_xxy[k];

                g_y_z_yz_xy[k] = -g_0_z_z_xy[k] - g_y_z_z_xy[k] * ab_y + g_y_z_z_xyy[k];

                g_y_z_yz_xz[k] = -g_0_z_z_xz[k] - g_y_z_z_xz[k] * ab_y + g_y_z_z_xyz[k];

                g_y_z_yz_yy[k] = -g_0_z_z_yy[k] - g_y_z_z_yy[k] * ab_y + g_y_z_z_yyy[k];

                g_y_z_yz_yz[k] = -g_0_z_z_yz[k] - g_y_z_z_yz[k] * ab_y + g_y_z_z_yyz[k];

                g_y_z_yz_zz[k] = -g_0_z_z_zz[k] - g_y_z_z_zz[k] * ab_y + g_y_z_z_yzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_y_z_zz_xx = cbuffer.data(dd_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_z_zz_xy = cbuffer.data(dd_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_z_zz_xz = cbuffer.data(dd_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_z_zz_yy = cbuffer.data(dd_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_z_zz_yz = cbuffer.data(dd_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_z_zz_zz = cbuffer.data(dd_geom_11_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz, g_y_z_z_xx, g_y_z_z_xxz, g_y_z_z_xy, g_y_z_z_xyz, g_y_z_z_xz, g_y_z_z_xzz, g_y_z_z_yy, g_y_z_z_yyz, g_y_z_z_yz, g_y_z_z_yzz, g_y_z_z_zz, g_y_z_z_zzz, g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zz_xx[k] = g_y_0_z_xx[k] - g_y_z_z_xx[k] * ab_z + g_y_z_z_xxz[k];

                g_y_z_zz_xy[k] = g_y_0_z_xy[k] - g_y_z_z_xy[k] * ab_z + g_y_z_z_xyz[k];

                g_y_z_zz_xz[k] = g_y_0_z_xz[k] - g_y_z_z_xz[k] * ab_z + g_y_z_z_xzz[k];

                g_y_z_zz_yy[k] = g_y_0_z_yy[k] - g_y_z_z_yy[k] * ab_z + g_y_z_z_yyz[k];

                g_y_z_zz_yz[k] = g_y_0_z_yz[k] - g_y_z_z_yz[k] * ab_z + g_y_z_z_yzz[k];

                g_y_z_zz_zz[k] = g_y_0_z_zz[k] - g_y_z_z_zz[k] * ab_z + g_y_z_z_zzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_z_x_xx_xx = cbuffer.data(dd_geom_11_off + 216 * ccomps * dcomps);

            auto g_z_x_xx_xy = cbuffer.data(dd_geom_11_off + 217 * ccomps * dcomps);

            auto g_z_x_xx_xz = cbuffer.data(dd_geom_11_off + 218 * ccomps * dcomps);

            auto g_z_x_xx_yy = cbuffer.data(dd_geom_11_off + 219 * ccomps * dcomps);

            auto g_z_x_xx_yz = cbuffer.data(dd_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_x_xx_zz = cbuffer.data(dd_geom_11_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, g_z_x_x_xx, g_z_x_x_xxx, g_z_x_x_xxy, g_z_x_x_xxz, g_z_x_x_xy, g_z_x_x_xyy, g_z_x_x_xyz, g_z_x_x_xz, g_z_x_x_xzz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz, g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xx_xx[k] = g_z_0_x_xx[k] - g_z_x_x_xx[k] * ab_x + g_z_x_x_xxx[k];

                g_z_x_xx_xy[k] = g_z_0_x_xy[k] - g_z_x_x_xy[k] * ab_x + g_z_x_x_xxy[k];

                g_z_x_xx_xz[k] = g_z_0_x_xz[k] - g_z_x_x_xz[k] * ab_x + g_z_x_x_xxz[k];

                g_z_x_xx_yy[k] = g_z_0_x_yy[k] - g_z_x_x_yy[k] * ab_x + g_z_x_x_xyy[k];

                g_z_x_xx_yz[k] = g_z_0_x_yz[k] - g_z_x_x_yz[k] * ab_x + g_z_x_x_xyz[k];

                g_z_x_xx_zz[k] = g_z_0_x_zz[k] - g_z_x_x_zz[k] * ab_x + g_z_x_x_xzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_z_x_xy_xx = cbuffer.data(dd_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_x_xy_xy = cbuffer.data(dd_geom_11_off + 223 * ccomps * dcomps);

            auto g_z_x_xy_xz = cbuffer.data(dd_geom_11_off + 224 * ccomps * dcomps);

            auto g_z_x_xy_yy = cbuffer.data(dd_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_x_xy_yz = cbuffer.data(dd_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_x_xy_zz = cbuffer.data(dd_geom_11_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_x_xx, g_z_x_x_xxy, g_z_x_x_xy, g_z_x_x_xyy, g_z_x_x_xyz, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yyy, g_z_x_x_yyz, g_z_x_x_yz, g_z_x_x_yzz, g_z_x_x_zz, g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xy_xx[k] = -g_z_x_x_xx[k] * ab_y + g_z_x_x_xxy[k];

                g_z_x_xy_xy[k] = -g_z_x_x_xy[k] * ab_y + g_z_x_x_xyy[k];

                g_z_x_xy_xz[k] = -g_z_x_x_xz[k] * ab_y + g_z_x_x_xyz[k];

                g_z_x_xy_yy[k] = -g_z_x_x_yy[k] * ab_y + g_z_x_x_yyy[k];

                g_z_x_xy_yz[k] = -g_z_x_x_yz[k] * ab_y + g_z_x_x_yyz[k];

                g_z_x_xy_zz[k] = -g_z_x_x_zz[k] * ab_y + g_z_x_x_yzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_z_x_xz_xx = cbuffer.data(dd_geom_11_off + 228 * ccomps * dcomps);

            auto g_z_x_xz_xy = cbuffer.data(dd_geom_11_off + 229 * ccomps * dcomps);

            auto g_z_x_xz_xz = cbuffer.data(dd_geom_11_off + 230 * ccomps * dcomps);

            auto g_z_x_xz_yy = cbuffer.data(dd_geom_11_off + 231 * ccomps * dcomps);

            auto g_z_x_xz_yz = cbuffer.data(dd_geom_11_off + 232 * ccomps * dcomps);

            auto g_z_x_xz_zz = cbuffer.data(dd_geom_11_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz, g_z_x_z_xx, g_z_x_z_xxx, g_z_x_z_xxy, g_z_x_z_xxz, g_z_x_z_xy, g_z_x_z_xyy, g_z_x_z_xyz, g_z_x_z_xz, g_z_x_z_xzz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xz_xx[k] = g_z_0_z_xx[k] - g_z_x_z_xx[k] * ab_x + g_z_x_z_xxx[k];

                g_z_x_xz_xy[k] = g_z_0_z_xy[k] - g_z_x_z_xy[k] * ab_x + g_z_x_z_xxy[k];

                g_z_x_xz_xz[k] = g_z_0_z_xz[k] - g_z_x_z_xz[k] * ab_x + g_z_x_z_xxz[k];

                g_z_x_xz_yy[k] = g_z_0_z_yy[k] - g_z_x_z_yy[k] * ab_x + g_z_x_z_xyy[k];

                g_z_x_xz_yz[k] = g_z_0_z_yz[k] - g_z_x_z_yz[k] * ab_x + g_z_x_z_xyz[k];

                g_z_x_xz_zz[k] = g_z_0_z_zz[k] - g_z_x_z_zz[k] * ab_x + g_z_x_z_xzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_z_x_yy_xx = cbuffer.data(dd_geom_11_off + 234 * ccomps * dcomps);

            auto g_z_x_yy_xy = cbuffer.data(dd_geom_11_off + 235 * ccomps * dcomps);

            auto g_z_x_yy_xz = cbuffer.data(dd_geom_11_off + 236 * ccomps * dcomps);

            auto g_z_x_yy_yy = cbuffer.data(dd_geom_11_off + 237 * ccomps * dcomps);

            auto g_z_x_yy_yz = cbuffer.data(dd_geom_11_off + 238 * ccomps * dcomps);

            auto g_z_x_yy_zz = cbuffer.data(dd_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_y_xx, g_z_x_y_xxy, g_z_x_y_xy, g_z_x_y_xyy, g_z_x_y_xyz, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yyy, g_z_x_y_yyz, g_z_x_y_yz, g_z_x_y_yzz, g_z_x_y_zz, g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yy_xx[k] = -g_z_x_y_xx[k] * ab_y + g_z_x_y_xxy[k];

                g_z_x_yy_xy[k] = -g_z_x_y_xy[k] * ab_y + g_z_x_y_xyy[k];

                g_z_x_yy_xz[k] = -g_z_x_y_xz[k] * ab_y + g_z_x_y_xyz[k];

                g_z_x_yy_yy[k] = -g_z_x_y_yy[k] * ab_y + g_z_x_y_yyy[k];

                g_z_x_yy_yz[k] = -g_z_x_y_yz[k] * ab_y + g_z_x_y_yyz[k];

                g_z_x_yy_zz[k] = -g_z_x_y_zz[k] * ab_y + g_z_x_y_yzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_z_x_yz_xx = cbuffer.data(dd_geom_11_off + 240 * ccomps * dcomps);

            auto g_z_x_yz_xy = cbuffer.data(dd_geom_11_off + 241 * ccomps * dcomps);

            auto g_z_x_yz_xz = cbuffer.data(dd_geom_11_off + 242 * ccomps * dcomps);

            auto g_z_x_yz_yy = cbuffer.data(dd_geom_11_off + 243 * ccomps * dcomps);

            auto g_z_x_yz_yz = cbuffer.data(dd_geom_11_off + 244 * ccomps * dcomps);

            auto g_z_x_yz_zz = cbuffer.data(dd_geom_11_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz, g_z_x_z_xx, g_z_x_z_xxy, g_z_x_z_xy, g_z_x_z_xyy, g_z_x_z_xyz, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yyy, g_z_x_z_yyz, g_z_x_z_yz, g_z_x_z_yzz, g_z_x_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yz_xx[k] = -g_z_x_z_xx[k] * ab_y + g_z_x_z_xxy[k];

                g_z_x_yz_xy[k] = -g_z_x_z_xy[k] * ab_y + g_z_x_z_xyy[k];

                g_z_x_yz_xz[k] = -g_z_x_z_xz[k] * ab_y + g_z_x_z_xyz[k];

                g_z_x_yz_yy[k] = -g_z_x_z_yy[k] * ab_y + g_z_x_z_yyy[k];

                g_z_x_yz_yz[k] = -g_z_x_z_yz[k] * ab_y + g_z_x_z_yyz[k];

                g_z_x_yz_zz[k] = -g_z_x_z_zz[k] * ab_y + g_z_x_z_yzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_z_x_zz_xx = cbuffer.data(dd_geom_11_off + 246 * ccomps * dcomps);

            auto g_z_x_zz_xy = cbuffer.data(dd_geom_11_off + 247 * ccomps * dcomps);

            auto g_z_x_zz_xz = cbuffer.data(dd_geom_11_off + 248 * ccomps * dcomps);

            auto g_z_x_zz_yy = cbuffer.data(dd_geom_11_off + 249 * ccomps * dcomps);

            auto g_z_x_zz_yz = cbuffer.data(dd_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_x_zz_zz = cbuffer.data(dd_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_z_x_z_xx, g_z_x_z_xxz, g_z_x_z_xy, g_z_x_z_xyz, g_z_x_z_xz, g_z_x_z_xzz, g_z_x_z_yy, g_z_x_z_yyz, g_z_x_z_yz, g_z_x_z_yzz, g_z_x_z_zz, g_z_x_z_zzz, g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zz_xx[k] = -g_0_x_z_xx[k] - g_z_x_z_xx[k] * ab_z + g_z_x_z_xxz[k];

                g_z_x_zz_xy[k] = -g_0_x_z_xy[k] - g_z_x_z_xy[k] * ab_z + g_z_x_z_xyz[k];

                g_z_x_zz_xz[k] = -g_0_x_z_xz[k] - g_z_x_z_xz[k] * ab_z + g_z_x_z_xzz[k];

                g_z_x_zz_yy[k] = -g_0_x_z_yy[k] - g_z_x_z_yy[k] * ab_z + g_z_x_z_yyz[k];

                g_z_x_zz_yz[k] = -g_0_x_z_yz[k] - g_z_x_z_yz[k] * ab_z + g_z_x_z_yzz[k];

                g_z_x_zz_zz[k] = -g_0_x_z_zz[k] - g_z_x_z_zz[k] * ab_z + g_z_x_z_zzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_z_y_xx_xx = cbuffer.data(dd_geom_11_off + 252 * ccomps * dcomps);

            auto g_z_y_xx_xy = cbuffer.data(dd_geom_11_off + 253 * ccomps * dcomps);

            auto g_z_y_xx_xz = cbuffer.data(dd_geom_11_off + 254 * ccomps * dcomps);

            auto g_z_y_xx_yy = cbuffer.data(dd_geom_11_off + 255 * ccomps * dcomps);

            auto g_z_y_xx_yz = cbuffer.data(dd_geom_11_off + 256 * ccomps * dcomps);

            auto g_z_y_xx_zz = cbuffer.data(dd_geom_11_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_x_xx, g_z_y_x_xxx, g_z_y_x_xxy, g_z_y_x_xxz, g_z_y_x_xy, g_z_y_x_xyy, g_z_y_x_xyz, g_z_y_x_xz, g_z_y_x_xzz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xx_xx[k] = -g_z_y_x_xx[k] * ab_x + g_z_y_x_xxx[k];

                g_z_y_xx_xy[k] = -g_z_y_x_xy[k] * ab_x + g_z_y_x_xxy[k];

                g_z_y_xx_xz[k] = -g_z_y_x_xz[k] * ab_x + g_z_y_x_xxz[k];

                g_z_y_xx_yy[k] = -g_z_y_x_yy[k] * ab_x + g_z_y_x_xyy[k];

                g_z_y_xx_yz[k] = -g_z_y_x_yz[k] * ab_x + g_z_y_x_xyz[k];

                g_z_y_xx_zz[k] = -g_z_y_x_zz[k] * ab_x + g_z_y_x_xzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_z_y_xy_xx = cbuffer.data(dd_geom_11_off + 258 * ccomps * dcomps);

            auto g_z_y_xy_xy = cbuffer.data(dd_geom_11_off + 259 * ccomps * dcomps);

            auto g_z_y_xy_xz = cbuffer.data(dd_geom_11_off + 260 * ccomps * dcomps);

            auto g_z_y_xy_yy = cbuffer.data(dd_geom_11_off + 261 * ccomps * dcomps);

            auto g_z_y_xy_yz = cbuffer.data(dd_geom_11_off + 262 * ccomps * dcomps);

            auto g_z_y_xy_zz = cbuffer.data(dd_geom_11_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz, g_z_y_y_xx, g_z_y_y_xxx, g_z_y_y_xxy, g_z_y_y_xxz, g_z_y_y_xy, g_z_y_y_xyy, g_z_y_y_xyz, g_z_y_y_xz, g_z_y_y_xzz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xy_xx[k] = -g_z_y_y_xx[k] * ab_x + g_z_y_y_xxx[k];

                g_z_y_xy_xy[k] = -g_z_y_y_xy[k] * ab_x + g_z_y_y_xxy[k];

                g_z_y_xy_xz[k] = -g_z_y_y_xz[k] * ab_x + g_z_y_y_xxz[k];

                g_z_y_xy_yy[k] = -g_z_y_y_yy[k] * ab_x + g_z_y_y_xyy[k];

                g_z_y_xy_yz[k] = -g_z_y_y_yz[k] * ab_x + g_z_y_y_xyz[k];

                g_z_y_xy_zz[k] = -g_z_y_y_zz[k] * ab_x + g_z_y_y_xzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_z_y_xz_xx = cbuffer.data(dd_geom_11_off + 264 * ccomps * dcomps);

            auto g_z_y_xz_xy = cbuffer.data(dd_geom_11_off + 265 * ccomps * dcomps);

            auto g_z_y_xz_xz = cbuffer.data(dd_geom_11_off + 266 * ccomps * dcomps);

            auto g_z_y_xz_yy = cbuffer.data(dd_geom_11_off + 267 * ccomps * dcomps);

            auto g_z_y_xz_yz = cbuffer.data(dd_geom_11_off + 268 * ccomps * dcomps);

            auto g_z_y_xz_zz = cbuffer.data(dd_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz, g_z_y_z_xx, g_z_y_z_xxx, g_z_y_z_xxy, g_z_y_z_xxz, g_z_y_z_xy, g_z_y_z_xyy, g_z_y_z_xyz, g_z_y_z_xz, g_z_y_z_xzz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xz_xx[k] = -g_z_y_z_xx[k] * ab_x + g_z_y_z_xxx[k];

                g_z_y_xz_xy[k] = -g_z_y_z_xy[k] * ab_x + g_z_y_z_xxy[k];

                g_z_y_xz_xz[k] = -g_z_y_z_xz[k] * ab_x + g_z_y_z_xxz[k];

                g_z_y_xz_yy[k] = -g_z_y_z_yy[k] * ab_x + g_z_y_z_xyy[k];

                g_z_y_xz_yz[k] = -g_z_y_z_yz[k] * ab_x + g_z_y_z_xyz[k];

                g_z_y_xz_zz[k] = -g_z_y_z_zz[k] * ab_x + g_z_y_z_xzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_z_y_yy_xx = cbuffer.data(dd_geom_11_off + 270 * ccomps * dcomps);

            auto g_z_y_yy_xy = cbuffer.data(dd_geom_11_off + 271 * ccomps * dcomps);

            auto g_z_y_yy_xz = cbuffer.data(dd_geom_11_off + 272 * ccomps * dcomps);

            auto g_z_y_yy_yy = cbuffer.data(dd_geom_11_off + 273 * ccomps * dcomps);

            auto g_z_y_yy_yz = cbuffer.data(dd_geom_11_off + 274 * ccomps * dcomps);

            auto g_z_y_yy_zz = cbuffer.data(dd_geom_11_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz, g_z_y_y_xx, g_z_y_y_xxy, g_z_y_y_xy, g_z_y_y_xyy, g_z_y_y_xyz, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yyy, g_z_y_y_yyz, g_z_y_y_yz, g_z_y_y_yzz, g_z_y_y_zz, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yy_xx[k] = g_z_0_y_xx[k] - g_z_y_y_xx[k] * ab_y + g_z_y_y_xxy[k];

                g_z_y_yy_xy[k] = g_z_0_y_xy[k] - g_z_y_y_xy[k] * ab_y + g_z_y_y_xyy[k];

                g_z_y_yy_xz[k] = g_z_0_y_xz[k] - g_z_y_y_xz[k] * ab_y + g_z_y_y_xyz[k];

                g_z_y_yy_yy[k] = g_z_0_y_yy[k] - g_z_y_y_yy[k] * ab_y + g_z_y_y_yyy[k];

                g_z_y_yy_yz[k] = g_z_0_y_yz[k] - g_z_y_y_yz[k] * ab_y + g_z_y_y_yyz[k];

                g_z_y_yy_zz[k] = g_z_0_y_zz[k] - g_z_y_y_zz[k] * ab_y + g_z_y_y_yzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_z_y_yz_xx = cbuffer.data(dd_geom_11_off + 276 * ccomps * dcomps);

            auto g_z_y_yz_xy = cbuffer.data(dd_geom_11_off + 277 * ccomps * dcomps);

            auto g_z_y_yz_xz = cbuffer.data(dd_geom_11_off + 278 * ccomps * dcomps);

            auto g_z_y_yz_yy = cbuffer.data(dd_geom_11_off + 279 * ccomps * dcomps);

            auto g_z_y_yz_yz = cbuffer.data(dd_geom_11_off + 280 * ccomps * dcomps);

            auto g_z_y_yz_zz = cbuffer.data(dd_geom_11_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz, g_z_y_z_xx, g_z_y_z_xxy, g_z_y_z_xy, g_z_y_z_xyy, g_z_y_z_xyz, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yyy, g_z_y_z_yyz, g_z_y_z_yz, g_z_y_z_yzz, g_z_y_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yz_xx[k] = g_z_0_z_xx[k] - g_z_y_z_xx[k] * ab_y + g_z_y_z_xxy[k];

                g_z_y_yz_xy[k] = g_z_0_z_xy[k] - g_z_y_z_xy[k] * ab_y + g_z_y_z_xyy[k];

                g_z_y_yz_xz[k] = g_z_0_z_xz[k] - g_z_y_z_xz[k] * ab_y + g_z_y_z_xyz[k];

                g_z_y_yz_yy[k] = g_z_0_z_yy[k] - g_z_y_z_yy[k] * ab_y + g_z_y_z_yyy[k];

                g_z_y_yz_yz[k] = g_z_0_z_yz[k] - g_z_y_z_yz[k] * ab_y + g_z_y_z_yyz[k];

                g_z_y_yz_zz[k] = g_z_0_z_zz[k] - g_z_y_z_zz[k] * ab_y + g_z_y_z_yzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_z_y_zz_xx = cbuffer.data(dd_geom_11_off + 282 * ccomps * dcomps);

            auto g_z_y_zz_xy = cbuffer.data(dd_geom_11_off + 283 * ccomps * dcomps);

            auto g_z_y_zz_xz = cbuffer.data(dd_geom_11_off + 284 * ccomps * dcomps);

            auto g_z_y_zz_yy = cbuffer.data(dd_geom_11_off + 285 * ccomps * dcomps);

            auto g_z_y_zz_yz = cbuffer.data(dd_geom_11_off + 286 * ccomps * dcomps);

            auto g_z_y_zz_zz = cbuffer.data(dd_geom_11_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_z_y_z_xx, g_z_y_z_xxz, g_z_y_z_xy, g_z_y_z_xyz, g_z_y_z_xz, g_z_y_z_xzz, g_z_y_z_yy, g_z_y_z_yyz, g_z_y_z_yz, g_z_y_z_yzz, g_z_y_z_zz, g_z_y_z_zzz, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zz_xx[k] = -g_0_y_z_xx[k] - g_z_y_z_xx[k] * ab_z + g_z_y_z_xxz[k];

                g_z_y_zz_xy[k] = -g_0_y_z_xy[k] - g_z_y_z_xy[k] * ab_z + g_z_y_z_xyz[k];

                g_z_y_zz_xz[k] = -g_0_y_z_xz[k] - g_z_y_z_xz[k] * ab_z + g_z_y_z_xzz[k];

                g_z_y_zz_yy[k] = -g_0_y_z_yy[k] - g_z_y_z_yy[k] * ab_z + g_z_y_z_yyz[k];

                g_z_y_zz_yz[k] = -g_0_y_z_yz[k] - g_z_y_z_yz[k] * ab_z + g_z_y_z_yzz[k];

                g_z_y_zz_zz[k] = -g_0_y_z_zz[k] - g_z_y_z_zz[k] * ab_z + g_z_y_z_zzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_z_z_xx_xx = cbuffer.data(dd_geom_11_off + 288 * ccomps * dcomps);

            auto g_z_z_xx_xy = cbuffer.data(dd_geom_11_off + 289 * ccomps * dcomps);

            auto g_z_z_xx_xz = cbuffer.data(dd_geom_11_off + 290 * ccomps * dcomps);

            auto g_z_z_xx_yy = cbuffer.data(dd_geom_11_off + 291 * ccomps * dcomps);

            auto g_z_z_xx_yz = cbuffer.data(dd_geom_11_off + 292 * ccomps * dcomps);

            auto g_z_z_xx_zz = cbuffer.data(dd_geom_11_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_x_xx, g_z_z_x_xxx, g_z_z_x_xxy, g_z_z_x_xxz, g_z_z_x_xy, g_z_z_x_xyy, g_z_z_x_xyz, g_z_z_x_xz, g_z_z_x_xzz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xx_xx[k] = -g_z_z_x_xx[k] * ab_x + g_z_z_x_xxx[k];

                g_z_z_xx_xy[k] = -g_z_z_x_xy[k] * ab_x + g_z_z_x_xxy[k];

                g_z_z_xx_xz[k] = -g_z_z_x_xz[k] * ab_x + g_z_z_x_xxz[k];

                g_z_z_xx_yy[k] = -g_z_z_x_yy[k] * ab_x + g_z_z_x_xyy[k];

                g_z_z_xx_yz[k] = -g_z_z_x_yz[k] * ab_x + g_z_z_x_xyz[k];

                g_z_z_xx_zz[k] = -g_z_z_x_zz[k] * ab_x + g_z_z_x_xzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_z_z_xy_xx = cbuffer.data(dd_geom_11_off + 294 * ccomps * dcomps);

            auto g_z_z_xy_xy = cbuffer.data(dd_geom_11_off + 295 * ccomps * dcomps);

            auto g_z_z_xy_xz = cbuffer.data(dd_geom_11_off + 296 * ccomps * dcomps);

            auto g_z_z_xy_yy = cbuffer.data(dd_geom_11_off + 297 * ccomps * dcomps);

            auto g_z_z_xy_yz = cbuffer.data(dd_geom_11_off + 298 * ccomps * dcomps);

            auto g_z_z_xy_zz = cbuffer.data(dd_geom_11_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz, g_z_z_y_xx, g_z_z_y_xxx, g_z_z_y_xxy, g_z_z_y_xxz, g_z_z_y_xy, g_z_z_y_xyy, g_z_z_y_xyz, g_z_z_y_xz, g_z_z_y_xzz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xy_xx[k] = -g_z_z_y_xx[k] * ab_x + g_z_z_y_xxx[k];

                g_z_z_xy_xy[k] = -g_z_z_y_xy[k] * ab_x + g_z_z_y_xxy[k];

                g_z_z_xy_xz[k] = -g_z_z_y_xz[k] * ab_x + g_z_z_y_xxz[k];

                g_z_z_xy_yy[k] = -g_z_z_y_yy[k] * ab_x + g_z_z_y_xyy[k];

                g_z_z_xy_yz[k] = -g_z_z_y_yz[k] * ab_x + g_z_z_y_xyz[k];

                g_z_z_xy_zz[k] = -g_z_z_y_zz[k] * ab_x + g_z_z_y_xzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_z_z_xz_xx = cbuffer.data(dd_geom_11_off + 300 * ccomps * dcomps);

            auto g_z_z_xz_xy = cbuffer.data(dd_geom_11_off + 301 * ccomps * dcomps);

            auto g_z_z_xz_xz = cbuffer.data(dd_geom_11_off + 302 * ccomps * dcomps);

            auto g_z_z_xz_yy = cbuffer.data(dd_geom_11_off + 303 * ccomps * dcomps);

            auto g_z_z_xz_yz = cbuffer.data(dd_geom_11_off + 304 * ccomps * dcomps);

            auto g_z_z_xz_zz = cbuffer.data(dd_geom_11_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz, g_z_z_z_xx, g_z_z_z_xxx, g_z_z_z_xxy, g_z_z_z_xxz, g_z_z_z_xy, g_z_z_z_xyy, g_z_z_z_xyz, g_z_z_z_xz, g_z_z_z_xzz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xz_xx[k] = -g_z_z_z_xx[k] * ab_x + g_z_z_z_xxx[k];

                g_z_z_xz_xy[k] = -g_z_z_z_xy[k] * ab_x + g_z_z_z_xxy[k];

                g_z_z_xz_xz[k] = -g_z_z_z_xz[k] * ab_x + g_z_z_z_xxz[k];

                g_z_z_xz_yy[k] = -g_z_z_z_yy[k] * ab_x + g_z_z_z_xyy[k];

                g_z_z_xz_yz[k] = -g_z_z_z_yz[k] * ab_x + g_z_z_z_xyz[k];

                g_z_z_xz_zz[k] = -g_z_z_z_zz[k] * ab_x + g_z_z_z_xzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_z_z_yy_xx = cbuffer.data(dd_geom_11_off + 306 * ccomps * dcomps);

            auto g_z_z_yy_xy = cbuffer.data(dd_geom_11_off + 307 * ccomps * dcomps);

            auto g_z_z_yy_xz = cbuffer.data(dd_geom_11_off + 308 * ccomps * dcomps);

            auto g_z_z_yy_yy = cbuffer.data(dd_geom_11_off + 309 * ccomps * dcomps);

            auto g_z_z_yy_yz = cbuffer.data(dd_geom_11_off + 310 * ccomps * dcomps);

            auto g_z_z_yy_zz = cbuffer.data(dd_geom_11_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_y_xx, g_z_z_y_xxy, g_z_z_y_xy, g_z_z_y_xyy, g_z_z_y_xyz, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yyy, g_z_z_y_yyz, g_z_z_y_yz, g_z_z_y_yzz, g_z_z_y_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yy_xx[k] = -g_z_z_y_xx[k] * ab_y + g_z_z_y_xxy[k];

                g_z_z_yy_xy[k] = -g_z_z_y_xy[k] * ab_y + g_z_z_y_xyy[k];

                g_z_z_yy_xz[k] = -g_z_z_y_xz[k] * ab_y + g_z_z_y_xyz[k];

                g_z_z_yy_yy[k] = -g_z_z_y_yy[k] * ab_y + g_z_z_y_yyy[k];

                g_z_z_yy_yz[k] = -g_z_z_y_yz[k] * ab_y + g_z_z_y_yyz[k];

                g_z_z_yy_zz[k] = -g_z_z_y_zz[k] * ab_y + g_z_z_y_yzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_z_z_yz_xx = cbuffer.data(dd_geom_11_off + 312 * ccomps * dcomps);

            auto g_z_z_yz_xy = cbuffer.data(dd_geom_11_off + 313 * ccomps * dcomps);

            auto g_z_z_yz_xz = cbuffer.data(dd_geom_11_off + 314 * ccomps * dcomps);

            auto g_z_z_yz_yy = cbuffer.data(dd_geom_11_off + 315 * ccomps * dcomps);

            auto g_z_z_yz_yz = cbuffer.data(dd_geom_11_off + 316 * ccomps * dcomps);

            auto g_z_z_yz_zz = cbuffer.data(dd_geom_11_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz, g_z_z_z_xx, g_z_z_z_xxy, g_z_z_z_xy, g_z_z_z_xyy, g_z_z_z_xyz, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yyy, g_z_z_z_yyz, g_z_z_z_yz, g_z_z_z_yzz, g_z_z_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yz_xx[k] = -g_z_z_z_xx[k] * ab_y + g_z_z_z_xxy[k];

                g_z_z_yz_xy[k] = -g_z_z_z_xy[k] * ab_y + g_z_z_z_xyy[k];

                g_z_z_yz_xz[k] = -g_z_z_z_xz[k] * ab_y + g_z_z_z_xyz[k];

                g_z_z_yz_yy[k] = -g_z_z_z_yy[k] * ab_y + g_z_z_z_yyy[k];

                g_z_z_yz_yz[k] = -g_z_z_z_yz[k] * ab_y + g_z_z_z_yyz[k];

                g_z_z_yz_zz[k] = -g_z_z_z_zz[k] * ab_y + g_z_z_z_yzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_z_z_zz_xx = cbuffer.data(dd_geom_11_off + 318 * ccomps * dcomps);

            auto g_z_z_zz_xy = cbuffer.data(dd_geom_11_off + 319 * ccomps * dcomps);

            auto g_z_z_zz_xz = cbuffer.data(dd_geom_11_off + 320 * ccomps * dcomps);

            auto g_z_z_zz_yy = cbuffer.data(dd_geom_11_off + 321 * ccomps * dcomps);

            auto g_z_z_zz_yz = cbuffer.data(dd_geom_11_off + 322 * ccomps * dcomps);

            auto g_z_z_zz_zz = cbuffer.data(dd_geom_11_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, g_z_z_z_xx, g_z_z_z_xxz, g_z_z_z_xy, g_z_z_z_xyz, g_z_z_z_xz, g_z_z_z_xzz, g_z_z_z_yy, g_z_z_z_yyz, g_z_z_z_yz, g_z_z_z_yzz, g_z_z_z_zz, g_z_z_z_zzz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zz_xx[k] = -g_0_z_z_xx[k] + g_z_0_z_xx[k] - g_z_z_z_xx[k] * ab_z + g_z_z_z_xxz[k];

                g_z_z_zz_xy[k] = -g_0_z_z_xy[k] + g_z_0_z_xy[k] - g_z_z_z_xy[k] * ab_z + g_z_z_z_xyz[k];

                g_z_z_zz_xz[k] = -g_0_z_z_xz[k] + g_z_0_z_xz[k] - g_z_z_z_xz[k] * ab_z + g_z_z_z_xzz[k];

                g_z_z_zz_yy[k] = -g_0_z_z_yy[k] + g_z_0_z_yy[k] - g_z_z_z_yy[k] * ab_z + g_z_z_z_yyz[k];

                g_z_z_zz_yz[k] = -g_0_z_z_yz[k] + g_z_0_z_yz[k] - g_z_z_z_yz[k] * ab_z + g_z_z_z_yzz[k];

                g_z_z_zz_zz[k] = -g_0_z_z_zz[k] + g_z_0_z_zz[k] - g_z_z_z_zz[k] * ab_z + g_z_z_z_zzz[k];
            }
        }
    }
}

} // erirec namespace

