#include "ElectronRepulsionPrimRecSSSS.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_ssss(CSimdArray<double>& prim_buffer_0_ssss,
                                  const double* fovl_abcd,
                                  const double* bf_values) -> void
{
    const auto ndims = prim_buffer_0_ssss.number_of_columns();

    /// Set up components of targeted buffer : prim_buffer_0_ssss

    auto g_0_0_0_0_0 = prim_buffer_0_ssss[0];

    #pragma omp simd aligned(g_0_0_0_0_0, fovl_abcd, bf_values : 64)
    for (int i = 0; i < ndims; i++)
    {
        g_0_0_0_0_0[i] = fovl_abcd[i] * bf_values[i];
    }
}

} // erirec namespace

