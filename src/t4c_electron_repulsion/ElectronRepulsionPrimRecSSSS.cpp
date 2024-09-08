#include "ElectronRepulsionPrimRecSSSS.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_ssss(CSimdArray<double>&       pbuffer,
                                  const size_t              idx_eri_0_ssss,
                                  CSimdArray<double>&       factors,
                                  const size_t              idx_ovl,
                                  const CSimdArray<double>& bf_data,
                                  const size_t              idx_bvals) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    /// Set up components of targeted buffer : prim_buffer_0_ssss

    auto g_0_0_0_0_0 = pbuffer.data(idx_eri_0_ssss);

    /// Set up overlap factors

    auto fovl_abcd = factors.data(idx_ovl);

    /// Set up Boys function values

    auto bf_values = bf_data.data(idx_bvals);

#pragma omp simd aligned(g_0_0_0_0_0, fovl_abcd, bf_values : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_0_0_0_0[i] = fovl_abcd[i] * bf_values[i];
    }
}

}  // namespace erirec
