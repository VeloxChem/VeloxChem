#include "ThreeCenterElectronRepulsionPrimRecSSS.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_sss(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_eri_0_sss,
                                 CSimdArray<double>&       factors,
                                 const size_t              idx_ovl,
                                 const CSimdArray<double>& bf_data,
                                 const size_t              idx_bvals) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    /// Set up components of targeted buffer : prim_buffer_0_sss

    auto g_0_0_0_0 = pbuffer.data(idx_eri_0_sss);

    /// Set up overlap factors

    auto fovl_acd = factors.data(idx_ovl);

    /// Set up Boys function values

    auto bf_values = bf_data.data(idx_bvals);

#pragma omp simd aligned(g_0_0_0_0, fovl_acd, bf_values : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_0_0_0[i] = fovl_acd[i] * bf_values[i];
    }
    
}

} // t3ceri namespace

