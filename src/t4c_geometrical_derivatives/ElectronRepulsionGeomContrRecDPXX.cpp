#include "ElectronRepulsionGeomContrRecDPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto comp_bra_geom_hrr_electron_repulsion_dpxx(CSimdArray<double>&   cbuffer,
                                               const size_t          idx_geom_dpxx,
                                               const size_t          idx_geom_ppxx,
                                               const size_t          idx_geom_pdxx,
                                               const size_t          idx_ppxx,
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
            /// Set up components of auxilary buffer : PDSS

            auto pd_geom_off = idx_geom_pdxx + i * dcomps + j;

            auto g_x_x_xx = cbuffer.data(pd_geom_off + 0 * ccomps * dcomps);

            auto g_x_x_xy = cbuffer.data(pd_geom_off + 1 * ccomps * dcomps);

            auto g_x_x_xz = cbuffer.data(pd_geom_off + 2 * ccomps * dcomps);

            auto g_x_x_yy = cbuffer.data(pd_geom_off + 3 * ccomps * dcomps);

            auto g_x_x_yz = cbuffer.data(pd_geom_off + 4 * ccomps * dcomps);

            auto g_x_x_zz = cbuffer.data(pd_geom_off + 5 * ccomps * dcomps);

            auto g_x_y_xx = cbuffer.data(pd_geom_off + 6 * ccomps * dcomps);

            auto g_x_y_xy = cbuffer.data(pd_geom_off + 7 * ccomps * dcomps);

            auto g_x_y_xz = cbuffer.data(pd_geom_off + 8 * ccomps * dcomps);

            auto g_x_y_yy = cbuffer.data(pd_geom_off + 9 * ccomps * dcomps);

            auto g_x_y_yz = cbuffer.data(pd_geom_off + 10 * ccomps * dcomps);

            auto g_x_y_zz = cbuffer.data(pd_geom_off + 11 * ccomps * dcomps);

            auto g_x_z_xx = cbuffer.data(pd_geom_off + 12 * ccomps * dcomps);

            auto g_x_z_xy = cbuffer.data(pd_geom_off + 13 * ccomps * dcomps);

            auto g_x_z_xz = cbuffer.data(pd_geom_off + 14 * ccomps * dcomps);

            auto g_x_z_yy = cbuffer.data(pd_geom_off + 15 * ccomps * dcomps);

            auto g_x_z_yz = cbuffer.data(pd_geom_off + 16 * ccomps * dcomps);

            auto g_x_z_zz = cbuffer.data(pd_geom_off + 17 * ccomps * dcomps);
            
            pd_geom_off += 18 * ccomps * dcomps;
            
            auto g_y_x_xx = cbuffer.data(pd_geom_off + 0 * ccomps * dcomps);

            auto g_y_x_xy = cbuffer.data(pd_geom_off + 1 * ccomps * dcomps);

            auto g_y_x_xz = cbuffer.data(pd_geom_off + 2 * ccomps * dcomps);

            auto g_y_x_yy = cbuffer.data(pd_geom_off + 3 * ccomps * dcomps);

            auto g_y_x_yz = cbuffer.data(pd_geom_off + 4 * ccomps * dcomps);

            auto g_y_x_zz = cbuffer.data(pd_geom_off + 5 * ccomps * dcomps);

            auto g_y_y_xx = cbuffer.data(pd_geom_off + 6 * ccomps * dcomps);

            auto g_y_y_xy = cbuffer.data(pd_geom_off + 7 * ccomps * dcomps);

            auto g_y_y_xz = cbuffer.data(pd_geom_off + 8 * ccomps * dcomps);

            auto g_y_y_yy = cbuffer.data(pd_geom_off + 9 * ccomps * dcomps);

            auto g_y_y_yz = cbuffer.data(pd_geom_off + 10 * ccomps * dcomps);

            auto g_y_y_zz = cbuffer.data(pd_geom_off + 11 * ccomps * dcomps);

            auto g_y_z_xx = cbuffer.data(pd_geom_off + 12 * ccomps * dcomps);

            auto g_y_z_xy = cbuffer.data(pd_geom_off + 13 * ccomps * dcomps);

            auto g_y_z_xz = cbuffer.data(pd_geom_off + 14 * ccomps * dcomps);

            auto g_y_z_yy = cbuffer.data(pd_geom_off + 15 * ccomps * dcomps);

            auto g_y_z_yz = cbuffer.data(pd_geom_off + 16 * ccomps * dcomps);

            auto g_y_z_zz = cbuffer.data(pd_geom_off + 17 * ccomps * dcomps);
            
            pd_geom_off += 18 * ccomps * dcomps;
            
            auto g_z_x_xx = cbuffer.data(pd_geom_off + 0 * ccomps * dcomps);

            auto g_z_x_xy = cbuffer.data(pd_geom_off + 1 * ccomps * dcomps);

            auto g_z_x_xz = cbuffer.data(pd_geom_off + 2 * ccomps * dcomps);

            auto g_z_x_yy = cbuffer.data(pd_geom_off + 3 * ccomps * dcomps);

            auto g_z_x_yz = cbuffer.data(pd_geom_off + 4 * ccomps * dcomps);

            auto g_z_x_zz = cbuffer.data(pd_geom_off + 5 * ccomps * dcomps);

            auto g_z_y_xx = cbuffer.data(pd_geom_off + 6 * ccomps * dcomps);

            auto g_z_y_xy = cbuffer.data(pd_geom_off + 7 * ccomps * dcomps);

            auto g_z_y_xz = cbuffer.data(pd_geom_off + 8 * ccomps * dcomps);

            auto g_z_y_yy = cbuffer.data(pd_geom_off + 9 * ccomps * dcomps);

            auto g_z_y_yz = cbuffer.data(pd_geom_off + 10 * ccomps * dcomps);

            auto g_z_y_zz = cbuffer.data(pd_geom_off + 11 * ccomps * dcomps);

            auto g_z_z_xx = cbuffer.data(pd_geom_off + 12 * ccomps * dcomps);

            auto g_z_z_xy = cbuffer.data(pd_geom_off + 13 * ccomps * dcomps);

            auto g_z_z_xz = cbuffer.data(pd_geom_off + 14 * ccomps * dcomps);

            auto g_z_z_yy = cbuffer.data(pd_geom_off + 15 * ccomps * dcomps);

            auto g_z_z_yz = cbuffer.data(pd_geom_off + 16 * ccomps * dcomps);

            auto g_z_z_zz = cbuffer.data(pd_geom_off + 17 * ccomps * dcomps);
            
            /// Set up components of auxilary buffer : PPSS

            auto pp_geom_off = idx_geom_ppxx + i * dcomps + j;

            auto g_x_x_x = cbuffer.data(pp_geom_off + 0 * ccomps * dcomps);

            auto g_x_x_y = cbuffer.data(pp_geom_off + 1 * ccomps * dcomps);

            auto g_x_x_z = cbuffer.data(pp_geom_off + 2 * ccomps * dcomps);
            
            auto g_x_y_x = cbuffer.data(pp_geom_off + 3 * ccomps * dcomps);

            auto g_x_y_y = cbuffer.data(pp_geom_off + 4 * ccomps * dcomps);

            auto g_x_y_z = cbuffer.data(pp_geom_off + 5 * ccomps * dcomps);
            
            auto g_x_z_x = cbuffer.data(pp_geom_off + 6 * ccomps * dcomps);

            auto g_x_z_y = cbuffer.data(pp_geom_off + 7 * ccomps * dcomps);

            auto g_x_z_z = cbuffer.data(pp_geom_off + 8 * ccomps * dcomps);
            
            pp_geom_off += 9 * ccomps * dcomps;
            
            auto g_y_x_x = cbuffer.data(pp_geom_off + 0 * ccomps * dcomps);

            auto g_y_x_y = cbuffer.data(pp_geom_off + 1 * ccomps * dcomps);

            auto g_y_x_z = cbuffer.data(pp_geom_off + 2 * ccomps * dcomps);
            
            auto g_y_y_x = cbuffer.data(pp_geom_off + 3 * ccomps * dcomps);

            auto g_y_y_y = cbuffer.data(pp_geom_off + 4 * ccomps * dcomps);

            auto g_y_y_z = cbuffer.data(pp_geom_off + 5 * ccomps * dcomps);
            
            auto g_y_z_x = cbuffer.data(pp_geom_off + 6 * ccomps * dcomps);

            auto g_y_z_y = cbuffer.data(pp_geom_off + 7 * ccomps * dcomps);

            auto g_y_z_z = cbuffer.data(pp_geom_off + 8 * ccomps * dcomps);
            
            pp_geom_off += 9 * ccomps * dcomps;
            
            auto g_z_x_x = cbuffer.data(pp_geom_off + 0 * ccomps * dcomps);
            
            auto g_z_x_y = cbuffer.data(pp_geom_off + 1 * ccomps * dcomps);

            auto g_z_x_z = cbuffer.data(pp_geom_off + 2 * ccomps * dcomps);
            
            auto g_z_y_x = cbuffer.data(pp_geom_off + 3 * ccomps * dcomps);

            auto g_z_y_y = cbuffer.data(pp_geom_off + 4 * ccomps * dcomps);

            auto g_z_y_z = cbuffer.data(pp_geom_off + 5 * ccomps * dcomps);
            
            auto g_z_z_x = cbuffer.data(pp_geom_off + 6 * ccomps * dcomps);

            auto g_z_z_y = cbuffer.data(pp_geom_off + 7 * ccomps * dcomps);

            auto g_z_z_z = cbuffer.data(pp_geom_off + 8 * ccomps * dcomps);
            
            /// Set up components of auxilary buffer : PPSS

            auto pp_off = idx_ppxx + i * dcomps + j;

            auto g_x_x = cbuffer.data(pp_off + 0 * ccomps * dcomps);

            auto g_x_y = cbuffer.data(pp_off + 1 * ccomps * dcomps);

            auto g_x_z = cbuffer.data(pp_off + 2 * ccomps * dcomps);
            
            auto g_y_x = cbuffer.data(pp_off + 3 * ccomps * dcomps);

            auto g_y_y = cbuffer.data(pp_off + 4 * ccomps * dcomps);

            auto g_y_z = cbuffer.data(pp_off + 5 * ccomps * dcomps);
            
            auto g_z_x = cbuffer.data(pp_off + 6 * ccomps * dcomps);

            auto g_z_y = cbuffer.data(pp_off + 7 * ccomps * dcomps);

            auto g_z_z = cbuffer.data(pp_off + 8 * ccomps * dcomps);
            
            /// set up bra offset for contr_buffer_ddxx

            auto dp_geom_off = idx_geom_dpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_xx_x = cbuffer.data(dp_geom_off + 0 * ccomps * dcomps);

            auto g_x_xx_y = cbuffer.data(dp_geom_off + 1 * ccomps * dcomps);

            auto g_x_xx_z = cbuffer.data(dp_geom_off + 2 * ccomps * dcomps);

#pragma omp simd aligned(g_x_xx_x, g_x_xx_y, g_x_xx_z, g_x_x_xx, g_x_x_xy, g_x_x_xz, g_x_x, g_x_y, g_x_z : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_xx_x[k] = -g_x_x_x[k] * ab_x + g_x_x_xx[k] - g_x_x[k];
                
                g_x_xx_y[k] = -g_x_x_y[k] * ab_x + g_x_x_xy[k] - g_x_y[k];
                
                g_x_xx_z[k] = -g_x_x_z[k] * ab_x + g_x_x_xz[k] - g_x_z[k];
            }
            
            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_xy_x = cbuffer.data(dp_geom_off + 3 * ccomps * dcomps);

            auto g_x_xy_y = cbuffer.data(dp_geom_off + 4 * ccomps * dcomps);

            auto g_x_xy_z = cbuffer.data(dp_geom_off + 5 * ccomps * dcomps);

#pragma omp simd aligned(g_x_xy_x, g_x_xy_y, g_x_xy_z, g_x_x_xy, g_x_x_yy, g_x_x_yz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_xy_x[k] = -g_x_x_x[k] * ab_y + g_x_x_xy[k];
                
                g_x_xy_y[k] = -g_x_x_y[k] * ab_y + g_x_x_yy[k];
                
                g_x_xy_z[k] = -g_x_x_z[k] * ab_y + g_x_x_yz[k];
            }
            
            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_xz_x = cbuffer.data(dp_geom_off + 6 * ccomps * dcomps);

            auto g_x_xz_y = cbuffer.data(dp_geom_off + 7 * ccomps * dcomps);

            auto g_x_xz_z = cbuffer.data(dp_geom_off + 8 * ccomps * dcomps);

#pragma omp simd aligned(g_x_xz_x, g_x_xz_y, g_x_xz_z, g_x_x_xz, g_x_x_yz, g_x_x_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_xz_x[k] = -g_x_x_x[k] * ab_z + g_x_x_xz[k];
                
                g_x_xz_y[k] = -g_x_x_y[k] * ab_z + g_x_x_yz[k];
                
                g_x_xz_z[k] = -g_x_x_z[k] * ab_z + g_x_x_zz[k];
            }
            
            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_yy_x = cbuffer.data(dp_geom_off + 9 * ccomps * dcomps);

            auto g_x_yy_y = cbuffer.data(dp_geom_off + 10 * ccomps * dcomps);

            auto g_x_yy_z = cbuffer.data(dp_geom_off + 11 * ccomps * dcomps);

#pragma omp simd aligned(g_x_yy_x, g_x_yy_y, g_x_yy_z, g_x_y_xy, g_x_y_yy, g_x_y_yz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_yy_x[k] = -g_x_y_x[k] * ab_y + g_x_y_xy[k];
                
                g_x_yy_y[k] = -g_x_y_y[k] * ab_y + g_x_y_yy[k];
                
                g_x_yy_z[k] = -g_x_y_z[k] * ab_y + g_x_y_yz[k];
            }
            
            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_yz_x = cbuffer.data(dp_geom_off + 12 * ccomps * dcomps);

            auto g_x_yz_y = cbuffer.data(dp_geom_off + 13 * ccomps * dcomps);

            auto g_x_yz_z = cbuffer.data(dp_geom_off + 14 * ccomps * dcomps);

#pragma omp simd aligned(g_x_yz_x, g_x_yz_y, g_x_yz_z, g_x_y_xz, g_x_y_yz, g_x_y_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_yz_x[k] = -g_x_y_x[k] * ab_z + g_x_y_xz[k];
                
                g_x_yz_y[k] = -g_x_y_y[k] * ab_z + g_x_y_yz[k];
                
                g_x_yz_z[k] = -g_x_y_z[k] * ab_z + g_x_y_zz[k];
            }
            
            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_zz_x = cbuffer.data(dp_geom_off + 15 * ccomps * dcomps);

            auto g_x_zz_y = cbuffer.data(dp_geom_off + 16 * ccomps * dcomps);

            auto g_x_zz_z = cbuffer.data(dp_geom_off + 17 * ccomps * dcomps);

#pragma omp simd aligned(g_x_zz_x, g_x_zz_y, g_x_zz_z, g_x_z_xz, g_x_z_yz, g_x_z_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_zz_x[k] = -g_x_z_x[k] * ab_z + g_x_z_xz[k];
                
                g_x_zz_y[k] = -g_x_z_y[k] * ab_z + g_x_z_yz[k];
                
                g_x_zz_z[k] = -g_x_z_z[k] * ab_z + g_x_z_zz[k];
            }
            
            // updated offset
            
            dp_geom_off += 18 * ccomps * dcomps;
            
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_y_xx_x = cbuffer.data(dp_geom_off + 0 * ccomps * dcomps);

            auto g_y_xx_y = cbuffer.data(dp_geom_off + 1 * ccomps * dcomps);

            auto g_y_xx_z = cbuffer.data(dp_geom_off + 2 * ccomps * dcomps);

#pragma omp simd aligned(g_y_xx_x, g_y_xx_y, g_y_xx_z, g_y_x_xx, g_y_x_xy, g_y_x_xz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_xx_x[k] = -g_y_x_x[k] * ab_x + g_y_x_xx[k];
                
                g_y_xx_y[k] = -g_y_x_y[k] * ab_x + g_y_x_xy[k];
                
                g_y_xx_z[k] = -g_y_x_z[k] * ab_x + g_y_x_xz[k];
            }
            
            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_y_xy_x = cbuffer.data(dp_geom_off + 3 * ccomps * dcomps);

            auto g_y_xy_y = cbuffer.data(dp_geom_off + 4 * ccomps * dcomps);

            auto g_y_xy_z = cbuffer.data(dp_geom_off + 5 * ccomps * dcomps);

#pragma omp simd aligned(g_y_xy_x, g_y_xy_y, g_y_xy_z, g_y_x_xy, g_y_x_yy, g_y_x_yz, g_x_x, g_x_y, g_x_z : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_xy_x[k] = -g_y_x_x[k] * ab_y + g_y_x_xy[k] - g_x_x[k];
                
                g_y_xy_y[k] = -g_y_x_y[k] * ab_y + g_y_x_yy[k] - g_x_y[k];
                
                g_y_xy_z[k] = -g_y_x_z[k] * ab_y + g_y_x_yz[k] - g_x_z[k];
            }
            
            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_y_xz_x = cbuffer.data(dp_geom_off + 6 * ccomps * dcomps);

            auto g_y_xz_y = cbuffer.data(dp_geom_off + 7 * ccomps * dcomps);

            auto g_y_xz_z = cbuffer.data(dp_geom_off + 8 * ccomps * dcomps);

#pragma omp simd aligned(g_y_xz_x, g_y_xz_y, g_y_xz_z, g_y_x_xz, g_y_x_yz, g_y_x_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_xz_x[k] = -g_y_x_x[k] * ab_z + g_y_x_xz[k];
                
                g_y_xz_y[k] = -g_y_x_y[k] * ab_z + g_y_x_yz[k];
                
                g_y_xz_z[k] = -g_y_x_z[k] * ab_z + g_y_x_zz[k];
            }
            
            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_y_yy_x = cbuffer.data(dp_geom_off + 9 * ccomps * dcomps);

            auto g_y_yy_y = cbuffer.data(dp_geom_off + 10 * ccomps * dcomps);

            auto g_y_yy_z = cbuffer.data(dp_geom_off + 11 * ccomps * dcomps);

#pragma omp simd aligned(g_y_yy_x, g_y_yy_y, g_y_yy_z, g_y_y_xy, g_y_y_yy, g_y_y_yz, g_y_x, g_y_y, g_y_z : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_yy_x[k] = -g_y_y_x[k] * ab_y + g_y_y_xy[k] - g_y_x[k];
                
                g_y_yy_y[k] = -g_y_y_y[k] * ab_y + g_y_y_yy[k] - g_y_y[k];
                
                g_y_yy_z[k] = -g_y_y_z[k] * ab_y + g_y_y_yz[k] - g_y_z[k];
            }
            
            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_y_yz_x = cbuffer.data(dp_geom_off + 12 * ccomps * dcomps);

            auto g_y_yz_y = cbuffer.data(dp_geom_off + 13 * ccomps * dcomps);

            auto g_y_yz_z = cbuffer.data(dp_geom_off + 14 * ccomps * dcomps);

#pragma omp simd aligned(g_y_yz_x, g_y_yz_y, g_y_yz_z, g_y_y_xz, g_y_y_yz, g_y_y_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_yz_x[k] = -g_y_y_x[k] * ab_z + g_y_y_xz[k];
                
                g_y_yz_y[k] = -g_y_y_y[k] * ab_z + g_y_y_yz[k];
                
                g_y_yz_z[k] = -g_y_y_z[k] * ab_z + g_y_y_zz[k];
            }
            
            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_y_zz_x = cbuffer.data(dp_geom_off + 15 * ccomps * dcomps);

            auto g_y_zz_y = cbuffer.data(dp_geom_off + 16 * ccomps * dcomps);

            auto g_y_zz_z = cbuffer.data(dp_geom_off + 17 * ccomps * dcomps);

#pragma omp simd aligned(g_y_zz_x, g_y_zz_y, g_y_zz_z, g_y_z_xz, g_y_z_yz, g_y_z_zz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_zz_x[k] = -g_y_z_x[k] * ab_z + g_y_z_xz[k];
                
                g_y_zz_y[k] = -g_y_z_y[k] * ab_z + g_y_z_yz[k];
                
                g_y_zz_z[k] = -g_y_z_z[k] * ab_z + g_y_z_zz[k];
            }
            
            // updated offset
            
            dp_geom_off += 18 * ccomps * dcomps;
            
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_z_xx_x = cbuffer.data(dp_geom_off + 0 * ccomps * dcomps);

            auto g_z_xx_y = cbuffer.data(dp_geom_off + 1 * ccomps * dcomps);

            auto g_z_xx_z = cbuffer.data(dp_geom_off + 2 * ccomps * dcomps);

#pragma omp simd aligned(g_z_xx_x, g_z_xx_y, g_z_xx_z, g_z_x_xx, g_z_x_xy, g_z_x_xz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_xx_x[k] = -g_z_x_x[k] * ab_x + g_z_x_xx[k];
                
                g_z_xx_y[k] = -g_z_x_y[k] * ab_x + g_z_x_xy[k];
                
                g_z_xx_z[k] = -g_z_x_z[k] * ab_x + g_z_x_xz[k];
            }
            
            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_z_xy_x = cbuffer.data(dp_geom_off + 3 * ccomps * dcomps);

            auto g_z_xy_y = cbuffer.data(dp_geom_off + 4 * ccomps * dcomps);

            auto g_z_xy_z = cbuffer.data(dp_geom_off + 5 * ccomps * dcomps);

#pragma omp simd aligned(g_z_xy_x, g_z_xy_y, g_z_xy_z, g_z_x_xy, g_z_x_yy, g_z_x_yz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_xy_x[k] = -g_z_x_x[k] * ab_y + g_z_x_xy[k];
                
                g_z_xy_y[k] = -g_z_x_y[k] * ab_y + g_z_x_yy[k];
                
                g_z_xy_z[k] = -g_z_x_z[k] * ab_y + g_z_x_yz[k];
            }
            
            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_z_xz_x = cbuffer.data(dp_geom_off + 6 * ccomps * dcomps);

            auto g_z_xz_y = cbuffer.data(dp_geom_off + 7 * ccomps * dcomps);

            auto g_z_xz_z = cbuffer.data(dp_geom_off + 8 * ccomps * dcomps);

#pragma omp simd aligned(g_z_xz_x, g_z_xz_y, g_z_xz_z, g_z_x_xz, g_z_x_yz, g_z_x_zz, g_x_x, g_x_y, g_x_z : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_xz_x[k] = -g_z_x_x[k] * ab_z + g_z_x_xz[k] - g_x_x[k];
                
                g_z_xz_y[k] = -g_z_x_y[k] * ab_z + g_z_x_yz[k] - g_x_y[k];
                
                g_z_xz_z[k] = -g_z_x_z[k] * ab_z + g_z_x_zz[k] - g_x_z[k];
            }
            
            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_z_yy_x = cbuffer.data(dp_geom_off + 9 * ccomps * dcomps);

            auto g_z_yy_y = cbuffer.data(dp_geom_off + 10 * ccomps * dcomps);

            auto g_z_yy_z = cbuffer.data(dp_geom_off + 11 * ccomps * dcomps);

#pragma omp simd aligned(g_z_yy_x, g_z_yy_y, g_z_yy_z, g_z_y_xy, g_z_y_yy, g_z_y_yz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_yy_x[k] = -g_z_y_x[k] * ab_y + g_z_y_xy[k];
                
                g_z_yy_y[k] = -g_z_y_y[k] * ab_y + g_z_y_yy[k];
                
                g_z_yy_z[k] = -g_z_y_z[k] * ab_y + g_z_y_yz[k];
            }
            
            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_z_yz_x = cbuffer.data(dp_geom_off + 12 * ccomps * dcomps);

            auto g_z_yz_y = cbuffer.data(dp_geom_off + 13 * ccomps * dcomps);

            auto g_z_yz_z = cbuffer.data(dp_geom_off + 14 * ccomps * dcomps);

#pragma omp simd aligned(g_z_yz_x, g_z_yz_y, g_z_yz_z, g_z_y_xz, g_z_y_yz, g_z_y_zz, g_y_x, g_y_y, g_y_z : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_yz_x[k] = -g_z_y_x[k] * ab_z + g_z_y_xz[k] - g_y_x[k];
                
                g_z_yz_y[k] = -g_z_y_y[k] * ab_z + g_z_y_yz[k] - g_y_y[k];
                
                g_z_yz_z[k] = -g_z_y_z[k] * ab_z + g_z_y_zz[k] - g_y_z[k];
            }
            
            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_z_zz_x = cbuffer.data(dp_geom_off + 15 * ccomps * dcomps);

            auto g_z_zz_y = cbuffer.data(dp_geom_off + 16 * ccomps * dcomps);

            auto g_z_zz_z = cbuffer.data(dp_geom_off + 17 * ccomps * dcomps);

#pragma omp simd aligned(g_z_zz_x, g_z_zz_y, g_z_zz_z, g_z_z_xz, g_z_z_yz, g_z_z_zz, g_z_x, g_z_y, g_z_z : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_zz_x[k] = -g_z_z_x[k] * ab_z + g_z_z_xz[k] - g_z_x[k];
                
                g_z_zz_y[k] = -g_z_z_y[k] * ab_z + g_z_z_yz[k] - g_z_y[k];
                
                g_z_zz_z[k] = -g_z_z_z[k] * ab_z + g_z_z_zz[k] - g_z_z[k];
            }
        }
    }
}

}  // namespace erirec
