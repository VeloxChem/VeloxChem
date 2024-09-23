#include "ElectronRepulsionGeomContrRecPSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto
comp_bra_geom_hrr_electron_repulsion_psxx(CSimdArray<double>&   cbuffer,
                                               const size_t          idx_geom_psxx,
                                               const size_t          idx_geom_ssxx,
                                               const size_t          idx_geom_spxx,
                                               const size_t          idx_ssxx,
                                               const TPoint<double>& r_ab,
                                               const int             c_angmom,
                                               const int             d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{
        c_angmom,
    });

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{
        d_angmom,
    });

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SSSS

            const auto ss_geom_off = idx_geom_ssxx + i * dcomps + j;

            auto g_x_0_0 = cbuffer.data(ss_geom_off);

            auto g_y_0_0 = cbuffer.data(ss_geom_off + ccomps * dcomps);

            auto g_z_0_0 = cbuffer.data(ss_geom_off + 2 * ccomps * dcomps);
            
            /// Set up components of auxilary buffer : SDSS

            auto sp_geom_off = idx_geom_spxx + i * dcomps + j;

            auto g_x_0_x = cbuffer.data(sp_geom_off + 0 * ccomps * dcomps);

            auto g_x_0_y = cbuffer.data(sp_geom_off + 1 * ccomps * dcomps);

            auto g_x_0_z = cbuffer.data(sp_geom_off + 2 * ccomps * dcomps);
            
            sp_geom_off += 3 * ccomps * dcomps;
            
            auto g_y_0_x = cbuffer.data(sp_geom_off + 0 * ccomps * dcomps);

            auto g_y_0_y = cbuffer.data(sp_geom_off + 1 * ccomps * dcomps);

            auto g_y_0_z = cbuffer.data(sp_geom_off + 2 * ccomps * dcomps);
            
            sp_geom_off += 3 * ccomps * dcomps;
            
            auto g_z_0_x = cbuffer.data(sp_geom_off + 0 * ccomps * dcomps);

            auto g_z_0_y = cbuffer.data(sp_geom_off + 1 * ccomps * dcomps);

            auto g_z_0_z = cbuffer.data(sp_geom_off + 2 * ccomps * dcomps);
            
            /// Set up components of auxilary buffer : SSSS

            auto g_0_0 = cbuffer.data(idx_ssxx + i * dcomps + j);

            /// set up bra offset for contr_buffer_psxx

            auto ps_geom_off = idx_geom_psxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_x_0 = cbuffer.data(ps_geom_off + 0 * ccomps * dcomps);

            auto g_x_y_0 = cbuffer.data(ps_geom_off + 1 * ccomps * dcomps);

            auto g_x_z_0 = cbuffer.data(ps_geom_off + 2 * ccomps * dcomps);

#pragma omp simd aligned(g_0_0, g_x_x_0, g_x_0_x, g_x_y_0, g_x_0_y, g_x_z_0, g_x_0_z, g_x_0_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_0[k] = -g_x_0_0[k] * ab_x + g_x_0_x[k] - g_0_0[k];

                g_x_y_0[k] = -g_x_0_0[k] * ab_y + g_x_0_y[k];

                g_x_z_0[k] = -g_x_0_0[k] * ab_z + g_x_0_z[k];
            }

            // updated offset
            
            ps_geom_off += 3 * ccomps * dcomps;
            
            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_y_x_0 = cbuffer.data(ps_geom_off + 0 * ccomps * dcomps);

            auto g_y_y_0 = cbuffer.data(ps_geom_off + 1 * ccomps * dcomps);

            auto g_y_z_0 = cbuffer.data(ps_geom_off + 2 * ccomps * dcomps);

#pragma omp simd aligned(g_0_0, g_y_x_0, g_y_0_x, g_y_y_0, g_y_0_y, g_y_z_0, g_y_0_z, g_y_0_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_0[k] = -g_y_0_0[k] * ab_x + g_y_0_x[k];

                g_y_y_0[k] = -g_y_0_0[k] * ab_y + g_y_0_y[k] - g_0_0[k];

                g_y_z_0[k] = -g_y_0_0[k] * ab_z + g_y_0_z[k];
            }
            
            // updated offset
            
            ps_geom_off += 3 * ccomps * dcomps;
            
            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_z_x_0 = cbuffer.data(ps_geom_off + 0 * ccomps * dcomps);

            auto g_z_y_0 = cbuffer.data(ps_geom_off + 1 * ccomps * dcomps);

            auto g_z_z_0 = cbuffer.data(ps_geom_off + 2 * ccomps * dcomps);

#pragma omp simd aligned(g_0_0, g_z_x_0, g_z_0_x, g_z_y_0, g_z_0_y, g_z_z_0, g_z_0_z, g_z_0_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_0[k] = -g_z_0_0[k] * ab_x + g_z_0_x[k];

                g_z_y_0[k] = -g_z_0_0[k] * ab_y + g_z_0_y[k];

                g_z_z_0[k] = -g_z_0_0[k] * ab_z + g_z_0_z[k] - g_0_0[k];
            }
        }
    }
}

}  // namespace erirec
