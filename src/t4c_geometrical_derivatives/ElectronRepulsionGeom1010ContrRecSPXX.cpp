#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_spxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_spxx,
                                              const size_t idx_geom_0010_spxx,
                                              const size_t idx_geom_0010_sdxx,
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
            /// Set up components of auxilary buffer : SPSS

            const auto sp_geom_0010_off = idx_geom_0010_spxx + i * dcomps + j;

            auto g_0_0_x_0_0_x = cbuffer.data(sp_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_0_y = cbuffer.data(sp_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_0_z = cbuffer.data(sp_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_y_0_0_x = cbuffer.data(sp_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_y_0_0_y = cbuffer.data(sp_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_y_0_0_z = cbuffer.data(sp_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_z_0_0_x = cbuffer.data(sp_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_z_0_0_y = cbuffer.data(sp_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_z_0_0_z = cbuffer.data(sp_geom_0010_off + 8 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SDSS

            const auto sd_geom_0010_off = idx_geom_0010_sdxx + i * dcomps + j;

            auto g_0_0_x_0_0_xx = cbuffer.data(sd_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_0_xy = cbuffer.data(sd_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_0_xz = cbuffer.data(sd_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_0_yy = cbuffer.data(sd_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_0_yz = cbuffer.data(sd_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_0_zz = cbuffer.data(sd_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_y_0_0_xx = cbuffer.data(sd_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_y_0_0_xy = cbuffer.data(sd_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_y_0_0_xz = cbuffer.data(sd_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_y_0_0_yy = cbuffer.data(sd_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_y_0_0_yz = cbuffer.data(sd_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_y_0_0_zz = cbuffer.data(sd_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_z_0_0_xx = cbuffer.data(sd_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_z_0_0_xy = cbuffer.data(sd_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_z_0_0_xz = cbuffer.data(sd_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_z_0_0_yy = cbuffer.data(sd_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_z_0_0_yz = cbuffer.data(sd_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_z_0_0_zz = cbuffer.data(sd_geom_0010_off + 17 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_spxx

            const auto sp_geom_1010_off = idx_geom_1010_spxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_0_x = cbuffer.data(sp_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_0_y = cbuffer.data(sp_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_0_z = cbuffer.data(sp_geom_1010_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_x, g_0_0_x_0_0_xx, g_0_0_x_0_0_xy, g_0_0_x_0_0_xz, g_0_0_x_0_0_y, g_0_0_x_0_0_z, g_x_0_x_0_0_x, g_x_0_x_0_0_y, g_x_0_x_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_0_x[k] = -g_0_0_x_0_0_x[k] * ab_x + g_0_0_x_0_0_xx[k];

                g_x_0_x_0_0_y[k] = -g_0_0_x_0_0_y[k] * ab_x + g_0_0_x_0_0_xy[k];

                g_x_0_x_0_0_z[k] = -g_0_0_x_0_0_z[k] * ab_x + g_0_0_x_0_0_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_0_x = cbuffer.data(sp_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_y_0_0_y = cbuffer.data(sp_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_y_0_0_z = cbuffer.data(sp_geom_1010_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_x, g_0_0_y_0_0_xx, g_0_0_y_0_0_xy, g_0_0_y_0_0_xz, g_0_0_y_0_0_y, g_0_0_y_0_0_z, g_x_0_y_0_0_x, g_x_0_y_0_0_y, g_x_0_y_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_0_x[k] = -g_0_0_y_0_0_x[k] * ab_x + g_0_0_y_0_0_xx[k];

                g_x_0_y_0_0_y[k] = -g_0_0_y_0_0_y[k] * ab_x + g_0_0_y_0_0_xy[k];

                g_x_0_y_0_0_z[k] = -g_0_0_y_0_0_z[k] * ab_x + g_0_0_y_0_0_xz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_0_x = cbuffer.data(sp_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_z_0_0_y = cbuffer.data(sp_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_z_0_0_z = cbuffer.data(sp_geom_1010_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_x, g_0_0_z_0_0_xx, g_0_0_z_0_0_xy, g_0_0_z_0_0_xz, g_0_0_z_0_0_y, g_0_0_z_0_0_z, g_x_0_z_0_0_x, g_x_0_z_0_0_y, g_x_0_z_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_0_x[k] = -g_0_0_z_0_0_x[k] * ab_x + g_0_0_z_0_0_xx[k];

                g_x_0_z_0_0_y[k] = -g_0_0_z_0_0_y[k] * ab_x + g_0_0_z_0_0_xy[k];

                g_x_0_z_0_0_z[k] = -g_0_0_z_0_0_z[k] * ab_x + g_0_0_z_0_0_xz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_0_x = cbuffer.data(sp_geom_1010_off + 9 * ccomps * dcomps);

            auto g_y_0_x_0_0_y = cbuffer.data(sp_geom_1010_off + 10 * ccomps * dcomps);

            auto g_y_0_x_0_0_z = cbuffer.data(sp_geom_1010_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_x, g_0_0_x_0_0_xy, g_0_0_x_0_0_y, g_0_0_x_0_0_yy, g_0_0_x_0_0_yz, g_0_0_x_0_0_z, g_y_0_x_0_0_x, g_y_0_x_0_0_y, g_y_0_x_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_0_x[k] = -g_0_0_x_0_0_x[k] * ab_y + g_0_0_x_0_0_xy[k];

                g_y_0_x_0_0_y[k] = -g_0_0_x_0_0_y[k] * ab_y + g_0_0_x_0_0_yy[k];

                g_y_0_x_0_0_z[k] = -g_0_0_x_0_0_z[k] * ab_y + g_0_0_x_0_0_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_0_x = cbuffer.data(sp_geom_1010_off + 12 * ccomps * dcomps);

            auto g_y_0_y_0_0_y = cbuffer.data(sp_geom_1010_off + 13 * ccomps * dcomps);

            auto g_y_0_y_0_0_z = cbuffer.data(sp_geom_1010_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_x, g_0_0_y_0_0_xy, g_0_0_y_0_0_y, g_0_0_y_0_0_yy, g_0_0_y_0_0_yz, g_0_0_y_0_0_z, g_y_0_y_0_0_x, g_y_0_y_0_0_y, g_y_0_y_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_0_x[k] = -g_0_0_y_0_0_x[k] * ab_y + g_0_0_y_0_0_xy[k];

                g_y_0_y_0_0_y[k] = -g_0_0_y_0_0_y[k] * ab_y + g_0_0_y_0_0_yy[k];

                g_y_0_y_0_0_z[k] = -g_0_0_y_0_0_z[k] * ab_y + g_0_0_y_0_0_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_0_x = cbuffer.data(sp_geom_1010_off + 15 * ccomps * dcomps);

            auto g_y_0_z_0_0_y = cbuffer.data(sp_geom_1010_off + 16 * ccomps * dcomps);

            auto g_y_0_z_0_0_z = cbuffer.data(sp_geom_1010_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_x, g_0_0_z_0_0_xy, g_0_0_z_0_0_y, g_0_0_z_0_0_yy, g_0_0_z_0_0_yz, g_0_0_z_0_0_z, g_y_0_z_0_0_x, g_y_0_z_0_0_y, g_y_0_z_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_0_x[k] = -g_0_0_z_0_0_x[k] * ab_y + g_0_0_z_0_0_xy[k];

                g_y_0_z_0_0_y[k] = -g_0_0_z_0_0_y[k] * ab_y + g_0_0_z_0_0_yy[k];

                g_y_0_z_0_0_z[k] = -g_0_0_z_0_0_z[k] * ab_y + g_0_0_z_0_0_yz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_0_x = cbuffer.data(sp_geom_1010_off + 18 * ccomps * dcomps);

            auto g_z_0_x_0_0_y = cbuffer.data(sp_geom_1010_off + 19 * ccomps * dcomps);

            auto g_z_0_x_0_0_z = cbuffer.data(sp_geom_1010_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_x, g_0_0_x_0_0_xz, g_0_0_x_0_0_y, g_0_0_x_0_0_yz, g_0_0_x_0_0_z, g_0_0_x_0_0_zz, g_z_0_x_0_0_x, g_z_0_x_0_0_y, g_z_0_x_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_0_x[k] = -g_0_0_x_0_0_x[k] * ab_z + g_0_0_x_0_0_xz[k];

                g_z_0_x_0_0_y[k] = -g_0_0_x_0_0_y[k] * ab_z + g_0_0_x_0_0_yz[k];

                g_z_0_x_0_0_z[k] = -g_0_0_x_0_0_z[k] * ab_z + g_0_0_x_0_0_zz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_0_x = cbuffer.data(sp_geom_1010_off + 21 * ccomps * dcomps);

            auto g_z_0_y_0_0_y = cbuffer.data(sp_geom_1010_off + 22 * ccomps * dcomps);

            auto g_z_0_y_0_0_z = cbuffer.data(sp_geom_1010_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_x, g_0_0_y_0_0_xz, g_0_0_y_0_0_y, g_0_0_y_0_0_yz, g_0_0_y_0_0_z, g_0_0_y_0_0_zz, g_z_0_y_0_0_x, g_z_0_y_0_0_y, g_z_0_y_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_0_x[k] = -g_0_0_y_0_0_x[k] * ab_z + g_0_0_y_0_0_xz[k];

                g_z_0_y_0_0_y[k] = -g_0_0_y_0_0_y[k] * ab_z + g_0_0_y_0_0_yz[k];

                g_z_0_y_0_0_z[k] = -g_0_0_y_0_0_z[k] * ab_z + g_0_0_y_0_0_zz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_0_x = cbuffer.data(sp_geom_1010_off + 24 * ccomps * dcomps);

            auto g_z_0_z_0_0_y = cbuffer.data(sp_geom_1010_off + 25 * ccomps * dcomps);

            auto g_z_0_z_0_0_z = cbuffer.data(sp_geom_1010_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_x, g_0_0_z_0_0_xz, g_0_0_z_0_0_y, g_0_0_z_0_0_yz, g_0_0_z_0_0_z, g_0_0_z_0_0_zz, g_z_0_z_0_0_x, g_z_0_z_0_0_y, g_z_0_z_0_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_0_x[k] = -g_0_0_z_0_0_x[k] * ab_z + g_0_0_z_0_0_xz[k];

                g_z_0_z_0_0_y[k] = -g_0_0_z_0_0_y[k] * ab_z + g_0_0_z_0_0_yz[k];

                g_z_0_z_0_0_z[k] = -g_0_0_z_0_0_z[k] * ab_z + g_0_0_z_0_0_zz[k];
            }
        }
    }
}

} // erirec namespace

